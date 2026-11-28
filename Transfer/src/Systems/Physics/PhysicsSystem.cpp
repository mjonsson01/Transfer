// File: Transfer/src/Systems/Physics/PhysicsSystem.cpp

#include "PhysicsSystem.hpp"
#include "DynamoEngine/Math/DeterministicMath.hpp"
#include "DynamoEngine/Math/DeterministicRandom.hpp"
#include "Systems/Physics/PhysicsConstants.hpp"

// What each random number is for, so two decisions about the same body never get the same number (see
// DynamoEngine/Math/DeterministicRandom.hpp). Changing these changes how every shatter plays out.
constexpr uint64_t RANDOM_PURPOSE_FRAGMENT_SPEED = 1;
constexpr uint64_t RANDOM_PURPOSE_FRAGMENT_SIZE = 2;

// The key all of one body's random choices this tick are hashed from: the universe's seed, the tick, and the body
// (its spawn ID and exact position, so two bodies breaking on the same tick never share their random numbers).
static uint64_t bodyRandomKey(const UniverseState& universe_state, const GravitationalBody& body)
{
    uint64_t key = DynamoEngine::combineHash(universe_state.random_seed, universe_state.tick_count);
    key = DynamoEngine::combineHash(key, static_cast<uint64_t>(static_cast<int64_t>(body.macroIdentifier)));
    key = DynamoEngine::combineHash(key, DynamoEngine::doubleBits(body.position.x_val));
    key = DynamoEngine::combineHash(key, DynamoEngine::doubleBits(body.position.y_val));
    return key;
}

PhysicsSystem::PhysicsSystem() : m_profiler(PROFILER_WINDOW_FRAMES)
{
    // The F5 overlay's rows, in this order
    m_profile.planet_collisions = m_profiler.addSection("Coll: planets");
    m_profile.planet_debris_collisions = m_profiler.addSection("Coll: planet-debris");
    m_profile.debris_sort = m_profiler.addSection("Sort debris by cell");
    // (Rows whose code doesn't run are left out of the overlay: "Grav: debris-debris" only shows while that gravity is
    // switched on.)
    m_profile.debris_grid = m_profiler.addSection("Coll: debris grid");
    m_profile.debris_pairs = m_profiler.addSection("Coll: debris pairs");
    m_profile.ship_collisions = m_profiler.addSection("Coll: ship");
    m_profile.verlet = m_profiler.addSection("Verlet (2 halves)");
    m_profile.planet_gravity = m_profiler.addSection("Grav: planets");
    m_profile.debris_planet_gravity = m_profiler.addSection("Grav: planet-debris");
    m_profile.debris_debris_gravity = m_profiler.addSection("Grav: debris-debris");
    m_profile.ship_gravity = m_profiler.addSection("Grav: ship");
    m_profile.walls_and_box = m_profiler.addSection("Walls & box");
    m_profile.cleanup = m_profiler.addSection("Cleanup");
    m_profile.whole_tick = m_profiler.addSection("Whole tick");
}

PhysicsSystem::~PhysicsSystem() {}

void PhysicsSystem::advanceUniverseState(UniverseState& universe_state, const ShipControls& ship_controls)
{
    // Only the real tick feeds the F5 profile: the trajectory predictions run the same functions on copies
    m_profiling_real_tick = true;
    {
        DynamoEngine::ProfileScope tick_timer(m_profiler, m_profile.whole_tick);
        // Debris sort, collisions, the ship's controls, then motion (Verlet + forces + walls), then cleanup.
        universe_state.tick_count += 1; // part of every random number hashed this tick (bodyRandomKey)

        // Neighbors next to each other in memory, so the debris pair pass reads nearby memory. FIRST, before anything
        // keeps particle indices for later in the tick (rememberCollisionBoxSides does). Every
        // DEBRIS_SORT_INTERVAL_TICKS ticks: the order stays nearly right in between.
        m_ticks_since_debris_sort += 1;
        if (m_ticks_since_debris_sort >= DEBRIS_SORT_INTERVAL_TICKS)
        {
            DynamoEngine::ProfileScope sort_timer(m_profiler, m_profile.debris_sort);
            particleGrid.sortByCell(universe_state.particles);
            m_ticks_since_debris_sort = 0;
        }

        // Collisions push overlapping bodies apart and could shove one across the collision box's edge, so remember
        // which side everything was on before them, and put anything shoved across back on its side right after
        rememberCollisionBoxSides(universe_state);
        handleCollisionEvents(universe_state);
        bounceOffCollisionBox(universe_state, BoxSides::BeforeCollisions);

        applyShipInputPhysics(universe_state, ship_controls);

        integrateUniversePhysics(universe_state);

        cleanupGravBodies(universe_state);
    } // the whole tick is recorded here, before profiling is switched off again
    m_profiling_real_tick = false;
    // One profiler "frame" per tick: the F5 numbers are per tick, however many ticks a frame runs to catch up
    m_profiler.endFrame();
}

void PhysicsSystem::integrateUniversePhysics(UniverseState& universe_state)
{
    // Update forces (grav, later will add electromagnetic)
    {
        DynamoEngine::ProfileScope timer(m_profiler, m_profile.verlet, m_profiling_real_tick);
        integrateVelocityVerletFirstHalfStep(universe_state);
    }
    recomputeAllForces(universe_state);
    {
        DynamoEngine::ProfileScope timer(m_profiler, m_profile.verlet, m_profiling_real_tick); // adds to the same row
        integrateVelocityVerletSecondHalfStep(universe_state);
    }
    bounceOffWorldWalls(universe_state); // walls are part of motion, so the trajectory predictions bounce too
    bounceOffCollisionBox(universe_state, BoxSides::StartOfStep);
}

// One axis of a wall bounce: keeps the body's edge inside [min_edge, max_edge] and turns it back if it was heading out
static void bounceInsideWalls(double& position, double& velocity, double radius, double min_edge, double max_edge)
{
    if (max_edge - min_edge <= 2.0 * radius)
    {
        position = (min_edge + max_edge) / 2.0; // wider than the world: park it in the middle
        velocity = 0.0;
        return;
    }
    if (position - radius < min_edge)
    {
        position = min_edge + radius;
        if (velocity < 0.0)
        {
            velocity = -velocity * WALL_RESTITUTION;
        }
    }
    else if (position + radius > max_edge)
    {
        position = max_edge - radius;
        if (velocity > 0.0)
        {
            velocity = -velocity * WALL_RESTITUTION;
        }
    }
}

void PhysicsSystem::bounceOffWorldWalls(UniverseState& universe_state)
{
    DynamoEngine::ProfileScope timer(m_profiler, m_profile.walls_and_box, m_profiling_real_tick);
    if (!universe_state.reflective_walls)
    {
        return;
    }
    const DynamoEngine::Vector2D& low = universe_state.world_min_corner;
    const DynamoEngine::Vector2D& high = universe_state.world_max_corner;
    for (std::vector<GravitationalBody>* bodies : {&universe_state.macro_bodies, &universe_state.particles})
    {
        for (GravitationalBody& body : *bodies)
        {
            bounceInsideWalls(body.position.x_val, body.velocity.x_val, body.radius, low.x_val, high.x_val);
            bounceInsideWalls(body.position.y_val, body.velocity.y_val, body.radius, low.y_val, high.y_val);
        }
    }
    if (universe_state.has_ship)
    {
        universe_state.ship.bounceInsideWalls(low, high, WALL_RESTITUTION);
    }
}

// A body against the collision box's edges, from whichever side it was on at the start of the tick (`previous`)
static void bounceOffBox(DynamoEngine::Vector2D& position, DynamoEngine::Vector2D& velocity,
                         const DynamoEngine::Vector2D& previous, double radius, const DynamoEngine::Vector2D& box_min,
                         const DynamoEngine::Vector2D& box_max)
{
    const bool was_inside = previous.x_val >= box_min.x_val && previous.x_val <= box_max.x_val &&
                            previous.y_val >= box_min.y_val && previous.y_val <= box_max.y_val;
    if (was_inside)
    {
        bounceInsideWalls(position.x_val, velocity.x_val, radius, box_min.x_val, box_max.x_val);
        bounceInsideWalls(position.y_val, velocity.y_val, radius, box_min.y_val, box_max.y_val);
        return;
    }

    // From outside: if it now overlaps the box, push it back out through the side it came from
    const bool overlaps = position.x_val + radius > box_min.x_val && position.x_val - radius < box_max.x_val &&
                          position.y_val + radius > box_min.y_val && position.y_val - radius < box_max.y_val;
    if (!overlaps)
    {
        return;
    }
    if (previous.x_val < box_min.x_val)
    {
        position.x_val = box_min.x_val - radius;
        if (velocity.x_val > 0.0)
            velocity.x_val = -velocity.x_val * WALL_RESTITUTION;
    }
    else if (previous.x_val > box_max.x_val)
    {
        position.x_val = box_max.x_val + radius;
        if (velocity.x_val < 0.0)
            velocity.x_val = -velocity.x_val * WALL_RESTITUTION;
    }
    else if (previous.y_val < box_min.y_val)
    {
        position.y_val = box_min.y_val - radius;
        if (velocity.y_val > 0.0)
            velocity.y_val = -velocity.y_val * WALL_RESTITUTION;
    }
    else
    {
        position.y_val = box_max.y_val + radius;
        if (velocity.y_val < 0.0)
            velocity.y_val = -velocity.y_val * WALL_RESTITUTION;
    }
}

void PhysicsSystem::rememberCollisionBoxSides(const UniverseState& universe_state)
{
    DynamoEngine::ProfileScope timer(m_profiler, m_profile.walls_and_box, m_profiling_real_tick);
    m_box_sides_macro_bodies.clear();
    m_box_sides_particles.clear();
    if (!universe_state.has_collision_box)
    {
        return;
    }
    for (const GravitationalBody& body : universe_state.macro_bodies)
    {
        m_box_sides_macro_bodies.push_back(body.position);
    }
    for (const GravitationalBody& body : universe_state.particles)
    {
        m_box_sides_particles.push_back(body.position);
    }
}

// Each body against the box, judged from `remembered` (BeforeCollisions) or from its previousPosition (StartOfStep).
// Bodies added since the sides were remembered (new debris) are judged from where they are.
static void bounceListOffBox(std::vector<GravitationalBody>& bodies,
                             const std::vector<DynamoEngine::Vector2D>& remembered, bool use_remembered,
                             const DynamoEngine::Vector2D& box_min, const DynamoEngine::Vector2D& box_max)
{
    for (size_t i = 0; i < bodies.size(); ++i)
    {
        GravitationalBody& body = bodies[i];
        DynamoEngine::Vector2D side_reference = body.previousPosition;
        if (use_remembered)
        {
            side_reference = (i < remembered.size()) ? remembered[i] : body.position;
        }
        bounceOffBox(body.position, body.velocity, side_reference, body.radius, box_min, box_max);
    }
}

void PhysicsSystem::bounceOffCollisionBox(UniverseState& universe_state, BoxSides sides)
{
    DynamoEngine::ProfileScope timer(m_profiler, m_profile.walls_and_box, m_profiling_real_tick);
    if (!universe_state.has_collision_box)
    {
        return;
    }
    const DynamoEngine::Vector2D& box_min = universe_state.collision_box_min_corner;
    const DynamoEngine::Vector2D& box_max = universe_state.collision_box_max_corner;
    const bool use_remembered = (sides == BoxSides::BeforeCollisions);
    bounceListOffBox(universe_state.macro_bodies, m_box_sides_macro_bodies, use_remembered, box_min, box_max);
    bounceListOffBox(universe_state.particles, m_box_sides_particles, use_remembered, box_min, box_max);
    if (universe_state.has_ship)
    {
        universe_state.ship.bounceOffBox(box_min, box_max, WALL_RESTITUTION);
    }
}

void PhysicsSystem::cleanUp()
{
    // Any necessary cleanup code for the physics system
}

void PhysicsSystem::updateGravBodyInstantiations(UniverseState& universe_state, const PlayerInput& player_input,
                                                 const SpawnSettings& spawn_settings)
{
    const PlayerCommands& commands = player_input.commands;
    if (commands.spawn)
    {
        if (commands.spawn_aim.is_cluster)
        {
            createParticleCluster(universe_state, commands.spawn_aim, spawn_settings);
        }
        else
        {
            createMacroBody(universe_state, commands.spawn_aim, spawn_settings);
        }
    }
    // Undo the newest spawn (plain Delete)
    if (commands.undo_newest_spawn)
    {
        removeNewestSpawn(universe_state);
    }
    // Wipe every body (Shift + Delete)
    if (commands.clear_all_bodies)
    {
        universe_state.macro_bodies.clear();
        universe_state.particles.clear();
    }
}
void PhysicsSystem::removeNewestSpawn(UniverseState& universe_state)
{
    std::vector<GravitationalBody>& macro_bodies = universe_state.macro_bodies;
    std::vector<GravitationalBody>& particles = universe_state.particles;

    // Every spawn gets the next ID and its debris inherits it, so the newest spawn still in the scene is simply the
    // highest ID any body still carries. (A spawn that merged into another body no longer exists, so it's skipped.)
    int newest_id = -1;
    for (const GravitationalBody& body : macro_bodies)
    {
        newest_id = std::max(newest_id, body.macroIdentifier);
    }
    for (const GravitationalBody& particle : particles)
    {
        newest_id = std::max(newest_id, particle.macroIdentifier);
    }
    if (newest_id < 0)
    {
        return; // nothing the player spawned is left
    }

    // Remove every body carrying that ID: the planet or cluster itself, and any debris it broke into.
    // Erased right away (not just marked), because this also has to work while time is stopped.
    std::erase_if(macro_bodies,
                  [newest_id](const GravitationalBody& body) { return body.macroIdentifier == newest_id; });
    std::erase_if(particles, [newest_id](const GravitationalBody& body) { return body.macroIdentifier == newest_id; });
}

// --------- COLLISION HANDLING --------- //

static inline CollisionInfo getCollisionInfo(const GravitationalBody& a, const GravitationalBody& b)
{
    DynamoEngine::Vector2D r_vector = b.position - a.position;

    // Quick rejection first, without the square root and division below: most pairs the grid hands out don't touch.
    // Compares squared distances, with a margin of a billionth so that pairs right at the edge still go to the exact
    // test below: the answer is always the same as the exact test's.
    double reach = a.radius + b.radius;
    if (r_vector.squareMagnitude() > reach * reach * (1.0 + 1e-9))
    {
        CollisionInfo no_contact = {};
        no_contact.shouldCollide = false;
        return no_contact;
    }

    double distance = r_vector.magnitude();
    DynamoEngine::Vector2D unit_normal_vector =
        (distance > 1e-8) ? (r_vector / distance) : DynamoEngine::Vector2D(1.0, 0.0);

    DynamoEngine::Vector2D relative_velocity_vector = b.velocity - a.velocity;
    double normal_speed = relative_velocity_vector.dot(unit_normal_vector);
    double abs_normal_speed = std::abs(normal_speed);
    bool should_collide = (distance < b.radius + a.radius);
    double approaching_speed = -normal_speed; // positive if bodies are moving towards each other
    bool should_blow_up = (approaching_speed >= MIN_SHATTER_SPEED);
    return {distance,         unit_normal_vector, relative_velocity_vector, normal_speed,
            abs_normal_speed, should_collide,     should_blow_up,           &a};
}

static inline GravitationalBodyPair pickMassPair(GravitationalBody& a, GravitationalBody& b)
{
    if (abs(a.mass) == abs(b.mass))
        return {&a, &b, 1.0, true};
    if (abs(a.mass) > abs(b.mass))
        return {&a, &b, abs(a.mass / b.mass), false};
    return {&b, &a, abs(b.mass / a.mass), false};
}

void PhysicsSystem::handleCollisionEvents(UniverseState& universe_state)
{
    // Shatters during these passes put their fragments in m_pending_fragments, not in particles
    handleMacroMacroCollisions(universe_state);
    handleMacroParticleCollisions(universe_state);
    // The ship before the debris pairs: debris the ship shoves into other debris is separated in this same tick, so
    // the single debris pass covers it (until 2026-11-26 a second debris pass ran after the ship instead)
    handleShipCollisions(universe_state);
    handleParticleParticleCollisions(universe_state);

    // All loops over particles are finished: now it's safe to grow the vector
    std::vector<GravitationalBody>& particles = universe_state.particles;
    particles.insert(particles.end(), m_pending_fragments.begin(), m_pending_fragments.end());
    m_pending_fragments.clear(); // clear() keeps the memory, so next tick's shatters don't allocate again
}

size_t PhysicsSystem::liveParticleCount(const UniverseState& universe_state) const
{
    return universe_state.particles.size() + m_pending_fragments.size();
}

void PhysicsSystem::handleMacroMacroCollisions(UniverseState& universe_state)
{
    DynamoEngine::ProfileScope timer(m_profiler, m_profile.planet_collisions, m_profiling_real_tick);
    std::vector<GravitationalBody>& macro_body_list = universe_state.macro_bodies;
    size_t num_macro_bodies = macro_body_list.size();

    for (size_t i = 0; i < num_macro_bodies; ++i)
    {
        GravitationalBody& first_body = macro_body_list[i];
        if (!first_body.isCollidable || first_body.isMacroGhost || first_body.isMarkedForDeletion)
        {
            continue;
        }

        // j = i + 1: each unordered pair is visited exactly once (previously visited twice)
        for (size_t j = i + 1; j < num_macro_bodies; ++j)
        {
            GravitationalBody& second_body = macro_body_list[j];
            if (!second_body.isCollidable || second_body.isMacroGhost || second_body.isMarkedForDeletion)
            {
                continue;
            }

            CollisionInfo collision_info = getCollisionInfo(first_body, second_body);
            if (!collision_info.shouldCollide)
            {
                continue;
            }

            GravitationalBodyPair grav_body_pair = pickMassPair(first_body, second_body);
            handleDynamicCollision(grav_body_pair, collision_info, universe_state);
            // first_body may have just shattered or been absorbed: it's gone, so it can't hit anything else this tick
            if (first_body.isMarkedForDeletion)
            {
                break;
            }
        }
    }
}

void PhysicsSystem::handleMacroParticleCollisions(UniverseState& universe_state)
{
    DynamoEngine::ProfileScope timer(m_profiler, m_profile.planet_debris_collisions, m_profiling_real_tick);
    std::vector<GravitationalBody>& particles = universe_state.particles;
    std::vector<GravitationalBody>& macro_bodies = universe_state.macro_bodies;

    for (auto& particle : particles)
    {
        // A macro ghost ignores planets and the ship but still hits debris. So a macro-ghost PARTICLE skips this whole
        // pass (it's all planets), while a macro-ghost PLANET still takes part (debris is what it collides with)
        if (!particle.isCollidable || particle.isMacroGhost || particle.isMarkedForDeletion)
        {
            continue;
        }

        for (auto& macro_body : macro_bodies)
        {
            if (!macro_body.isCollidable || macro_body.isMarkedForDeletion)
            {
                continue;
            }
            CollisionInfo collision_info = getCollisionInfo(particle, macro_body);
            if (!collision_info.shouldCollide)
            {
                continue;
            }

            GravitationalBodyPair grav_body_pair = pickMassPair(particle, macro_body);
            handleDynamicCollision(grav_body_pair, collision_info, universe_state);
            // The particle may have just been absorbed: it's gone, so it can't hit anything else this tick
            if (particle.isMarkedForDeletion)
            {
                break;
            }
        }
    }
}

void PhysicsSystem::handleParticleParticleCollisions(UniverseState& universe_state)
{
    std::vector<GravitationalBody>& particles = universe_state.particles;

    {
        DynamoEngine::ProfileScope timer(m_profiler, m_profile.debris_grid, m_profiling_real_tick);
        particleGrid.build(particles);
    }
    DynamoEngine::ProfileScope timer(m_profiler, m_profile.debris_pairs, m_profiling_real_tick); // the rest

    // Alternate the sweep direction every tick. The pairs are resolved one after another and each push changes what
    // the next pair sees, so a sweep that always runs the same way (left to right: the debris is sorted by position)
    // pushes everything a little in that direction every tick. Piles pack tighter, and a cloud can slowly pick up a
    // spin. Sweeping back the other way every other tick cancels that bias out.
    m_sweep_debris_backward = !m_sweep_debris_backward;

    std::vector<size_t> candidates;
    const size_t num_particles = particles.size();
    for (size_t step = 0; step < num_particles; ++step)
    {
        size_t i = m_sweep_debris_backward ? (num_particles - 1 - step) : step;
        resolveContactsOfParticle(i, m_sweep_debris_backward, universe_state, candidates);
    }
}

void PhysicsSystem::resolveContactsOfParticle(size_t i, bool backward, UniverseState& universe_state,
                                              std::vector<size_t>& candidates)
{
    std::vector<GravitationalBody>& particles = universe_state.particles;
    if (!particles[i].isCollidable || particles[i].isMarkedForDeletion)
    {
        return;
    }

    particleGrid.queryCandidates(i, candidates);
    const size_t candidate_count = candidates.size();
    for (size_t step = 0; step < candidate_count; ++step)
    {
        // Backward sweeps take the candidates in reverse too, so the whole order of pairs is reversed
        size_t j = backward ? candidates[candidate_count - 1 - step] : candidates[step];
        if (!particles[j].isCollidable || particles[j].isMarkedForDeletion)
        {
            continue;
        }

        CollisionInfo collision_info = getCollisionInfo(particles[i], particles[j]);
        if (!collision_info.shouldCollide)
        {
            continue;
        }

        GravitationalBodyPair grav_body_pair = pickMassPair(particles[i], particles[j]);
        handleDynamicCollision(grav_body_pair, collision_info, universe_state);
        // particles[i] may have just been absorbed: it's gone, so it can't hit anything else this tick
        if (particles[i].isMarkedForDeletion)
        {
            return;
        }
    }
}

void PhysicsSystem::handleShipCollisions(UniverseState& universe_state)
{
    DynamoEngine::ProfileScope timer(m_profiler, m_profile.ship_collisions, m_profiling_real_tick);
    if (!universe_state.has_ship)
    {
        return;
    }
    Starship& ship = universe_state.ship;
    double impact_energy_this_tick = 0.0;

    // Planets
    for (GravitationalBody& body : universe_state.macro_bodies)
    {
        if (!body.isCollidable || body.isMacroGhost || body.isMarkedForDeletion)
        {
            continue;
        }
        DynamoEngine::CircleContact contact =
            DynamoEngine::circleVsConvexPolygon(body.position, body.radius, ship.collisionPolygon());
        if (contact.touching)
        {
            impact_energy_this_tick += resolveShipContact(ship, body, contact);
        }
    }

    // Debris: a cheap distance check against the ship's bounding radius first, so only the few particles near the
    // ship get the full polygon test
    double ship_reach = ship.boundingRadius();
    for (GravitationalBody& particle : universe_state.particles)
    {
        // The ship counts as a big body, like a planet: macro ghosts pass through it (as macro-ghost planets do above)
        if (!particle.isCollidable || particle.isMacroGhost || particle.isMarkedForDeletion)
        {
            continue;
        }
        double reach = ship_reach + particle.radius;
        if ((particle.position - ship.center()).squareMagnitude() > reach * reach)
        {
            continue;
        }
        DynamoEngine::CircleContact contact =
            DynamoEngine::circleVsConvexPolygon(particle.position, particle.radius, ship.collisionPolygon());
        if (contact.touching)
        {
            impact_energy_this_tick += resolveShipContact(ship, particle, contact);
        }
    }

    ship.setLastImpactEnergy(impact_energy_this_tick);
}

double PhysicsSystem::resolveShipContact(Starship& ship, GravitationalBody& body,
                                         const DynamoEngine::CircleContact& contact)
{
    // Inverse masses say how easily each side is moved. A force-static body counts as infinitely heavy (0).
    double ship_inverse_mass = 1.0 / ship.mass();
    double body_inverse_mass = body.isForceStatic ? 0.0 : body.invMass;
    double total_inverse_mass = ship_inverse_mass + body_inverse_mass;
    if (total_inverse_mass <= 0.0)
    {
        return 0.0; // only possible with negative masses (see REWORK); nothing sensible to do
    }

    // 1. Push them apart so they just touch, each moving in proportion to how light it is.
    //    contact.normal points from the ship toward the body, so the ship moves against it and the body along it.
    ship.moveBy(contact.normal * (-contact.depth * ship_inverse_mass / total_inverse_mass));
    body.position += contact.normal * (contact.depth * body_inverse_mass / total_inverse_mass);

    // 2. Bounce, but only if they're still moving toward each other along the normal (same rule as the bodies')
    double normal_speed = (body.velocity - ship.velocity()).dot(contact.normal); // < 0: closing
    if (normal_speed >= 0.0)
    {
        return 0.0;
    }

    // The impulse that turns the closing speed around and keeps SHIP_RESTITUTION of it, equal and opposite on both
    // sides, so momentum is conserved
    double impulse_size = -(1.0 + SHIP_RESTITUTION) * normal_speed / total_inverse_mass;
    DynamoEngine::Vector2D impulse = contact.normal * impulse_size;
    body.velocity += impulse * body_inverse_mass;
    ship.applyImpulse(impulse * -1.0);

    // The kinetic energy the bounce absorbed: 1/2 * reduced mass * closing speed^2 * (1 - e^2).
    // The reduced mass 1 / (1/m_ship + 1/m_body) is about the ship's own mass against a planet, less against debris.
    double reduced_mass = 1.0 / total_inverse_mass;
    return 0.5 * reduced_mass * normal_speed * normal_speed * (1.0 - SHIP_RESTITUTION * SHIP_RESTITUTION);
}
void PhysicsSystem::handleDynamicCollision(GravitationalBodyPair& grav_body_pair, const CollisionInfo& collisionInfo,
                                           UniverseState& universe_state)
{
    GravitationalBody& heavier = *grav_body_pair.heavierBody;
    GravitationalBody& lighter = *grav_body_pair.lighterBody;

    if (heavier.isBounce || lighter.isBounce)
    {
        bounceUsingCollisionInfo(lighter, heavier, collisionInfo);
        return;
    }

    // Each body decides for itself whether a hit this hard breaks it: the lighter one if it's shatterable, the heavier
    // one only if it's shatterable AND the masses are comparable (a much heavier body shrugs the hit off).
    // Deciding per body means it can't matter which one is called "lighter" (for equal masses: whichever came first).
    bool lighter_shatters = lighter.isShatterable;
    bool heavier_shatters = heavier.isShatterable && grav_body_pair.ratio <= MUTUAL_SHATTER_MASS_RATIO_THRESHOLD;
    if (collisionInfo.shouldBlowUp && (lighter_shatters || heavier_shatters))
    {
        DynamoEngine::Vector2D toward_lighter = (lighter.position - heavier.position).normalize();
        DynamoEngine::Vector2D impact_point = heavier.position + toward_lighter * heavier.radius;

        if (heavier_shatters)
        {
            substituteWithParticlesFromImpact(heavier, m_pending_fragments, DEFAULT_FRAGMENT_COUNT, impact_point,
                                              bodyRandomKey(universe_state, heavier));
        }
        if (lighter_shatters)
        {
            substituteWithParticlesFromImpact(lighter, m_pending_fragments, DEFAULT_FRAGMENT_COUNT, impact_point,
                                              bodyRandomKey(universe_state, lighter));
        }
        return;
    }

    if (collisionInfo.absNormalSpeed >= MAX_ACCRETION_COLLISION_SPEED)
    {
        // Too fast to cleanly merge, not fast enough to shatter: always bounce.
        bounceUsingCollisionInfo(lighter, heavier, collisionInfo);
        return;
    }

    // Gentle contact: merge if there's enough size disparity to look right, else bounce.
    double accretion_ratio_threshold;
    if (heavier.isMacro && lighter.isMacro)
    {
        accretion_ratio_threshold = MIN_BODY_BODY_ACCRETION_THRESHOLD_RATIO;
    }
    else if (heavier.isMacro || lighter.isMacro)
    {
        accretion_ratio_threshold = MIN_BODY_PARTICLE_ACCRETION_THRESHOLD_RATIO;
    }
    else
    {
        accretion_ratio_threshold = MIN_PARTICLE_PARTICLE_ACCRETION_THRESHOLD_RATIO;
    }

    // The ratio compares the SIZES of the masses (pickMassPair uses abs), so opposite signs follow the same rule: small
    // negative debris eats a planet down (handleAccretion subtracts) until the planet is under 10x its size, then
    // bounces off it like any gentle contact.
    if (lighter.isAccretable && grav_body_pair.ratio >= accretion_ratio_threshold)
    {
        // Macro bodies crumble into fragments that then accrete individually; particles merge directly.
        // This is a very gentle collision compared to our normal explosion collision and requires that major mass
        // disparity so we decrease fragment density
        bool can_crumble = lighter.isMacro && lighter.isShatterable;
        if (can_crumble)
        {
            DynamoEngine::Vector2D toward_lighter = (lighter.position - heavier.position).normalize();
            DynamoEngine::Vector2D contact_point = heavier.position + toward_lighter * heavier.radius;
            substituteWithParticlesFromImpact(lighter, m_pending_fragments, DEFAULT_FRAGMENT_COUNT / 3, contact_point,
                                              bodyRandomKey(universe_state, lighter));
        }
        else if (heavier.isParticle)
        {
            // Particles can't accrete, so we just handle an elastic collision.
            bounceUsingCollisionInfo(lighter, heavier, collisionInfo);
        }
        else
        {
            handleAccretion(grav_body_pair);
        }
    }
    else
    {
        bounceUsingCollisionInfo(lighter, heavier, collisionInfo);
    }
}
// The unit vector pointing from `from` to `to`. Two bodies exactly on top of each other have no such direction, so a
// "random" one is used instead: any direction pushes them apart, and a different one per pair keeps a stack of
// exactly-overlapping pairs (e.g. two identical clusters spawned on the same spot) from all sliding the same way.
// Hashed from the two positions' exact bits, so the same pair always gets the same direction.
static DynamoEngine::Vector2D directionFromTo(const DynamoEngine::Vector2D& from, const DynamoEngine::Vector2D& to)
{
    DynamoEngine::Vector2D offset = to - from;
    double distance = offset.magnitude();
    if (firstWithinEpsilonOfSecond(distance, 0.0))
    {
        uint64_t hash = DynamoEngine::combineHash(DynamoEngine::doubleBits(from.x_val), DynamoEngine::doubleBits(from.y_val));
        hash = DynamoEngine::combineHash(hash, DynamoEngine::doubleBits(to.x_val));
        hash = DynamoEngine::combineHash(hash, DynamoEngine::doubleBits(to.y_val));
        double angle = DynamoEngine::hashToRange(hash, 0.0, DynamoEngine::TWO_PI);
        return {DynamoEngine::deterministicCos(angle), DynamoEngine::deterministicSin(angle)};
    }
    return offset / distance;
}

void PhysicsSystem::bounceUsingCollisionInfo(GravitationalBody& lighter, GravitationalBody& heavier,
                                             const CollisionInfo& collision_info)
{
    DynamoEngine::Vector2D lighter_to_heavier;
    if (firstWithinEpsilonOfSecond(collision_info.distance, 0.0))
    {
        // Exactly on top of each other: there is no direction to reuse, so directionFromTo picks a random one
        lighter_to_heavier = directionFromTo(lighter.position, heavier.position);
    }
    else if (collision_info.bodyA == &lighter)
    {
        lighter_to_heavier = collision_info.unitNormalVector; // it points from body A (the lighter) to the heavier
    }
    else
    {
        lighter_to_heavier = collision_info.unitNormalVector * -1.0; // body A is the heavier one: turn it around
    }
    handleElasticCollisions(lighter, heavier, collision_info.distance, lighter_to_heavier);
}

void PhysicsSystem::handleElasticCollisions(GravitationalBody& smallerBody, GravitationalBody& largerBody,
                                            double distance, const DynamoEngine::Vector2D& smaller_to_larger)
{
    if (smallerBody.isForceStatic && largerBody.isForceStatic)
    {
        const DynamoEngine::Vector2D& n = smaller_to_larger;

        double closing_speed = (smallerBody.velocity - largerBody.velocity).dot(n); // > 0: the gap is shrinking
        if (closing_speed > 0.0)
        {
            // Equal masses: the closing speed is turned around (times the loss factor), half the change on each
            DynamoEngine::Vector2D velocity_change = n * ((1.0 + ELASTIC_LOSS_FACTOR) * closing_speed / 2.0);
            smallerBody.velocity -= velocity_change;
            largerBody.velocity += velocity_change;
        }

        double penetration = (smallerBody.radius + largerBody.radius) - distance;
        if (penetration > 0.0)
        {
            smallerBody.position -= n * (penetration / 2.0);
            largerBody.position += n * (penetration / 2.0);
        }
        return;
    }
    if (smallerBody.isForceStatic != largerBody.isForceStatic)
    {
        // One static body: an immovable wall (infinite mass), so only the other body (dyn) changes
        GravitationalBody& dyn = smallerBody.isForceStatic ? largerBody : smallerBody;
        GravitationalBody& stat = smallerBody.isForceStatic ? smallerBody : largerBody;

        // From the static body to the other one
        DynamoEngine::Vector2D n = smallerBody.isForceStatic ? smaller_to_larger : smaller_to_larger * -1.0;

        // The closing speed RELATIVE to the static body: a static body can be drifting, so the other body's own
        // velocity isn't enough (a resting body hit by a drifting static one would otherwise feel nothing and get
        // passed through)
        double v_n = (dyn.velocity - stat.velocity).dot(n);
        if (v_n < 0.0)
        {
            dyn.velocity -= n * (1.0 + ELASTIC_LOSS_FACTOR) * v_n;
        }

        // Push the other body all the way out: the static one can't move, and a drifting static body would otherwise
        // keep sliding deeper into it every tick
        double penetration = (dyn.radius + stat.radius) - distance;
        if (penetration > 0.0)
        {
            dyn.position += n * penetration;
        }
        return;
    }

    if (firstWithinEpsilonOfSecond(smallerBody.mass, 0.0) || firstWithinEpsilonOfSecond(largerBody.mass, 0.0))
    {
        return;
    }

    const DynamoEngine::Vector2D& normal_vector = smaller_to_larger;
    double v_smaller_n = smallerBody.velocity.dot(normal_vector);
    double v_larger_n = largerBody.velocity.dot(normal_vector);
    // A bounce uses the SIZE of each mass (std::abs), so a negative-mass body bounces like a normal one of the same
    // size (its gravity still has the negative sign). With the signed masses this formula divides by
    // m_smaller + m_larger, which is 0 for -m hitting +m, and for unequal opposite-sign masses it multiplies the
    // speed on every contact (-2m hitting +m: about 3x per bounce), so a pile reached infinite speed within a second.
    double m_smaller = std::abs(smallerBody.mass);
    double m_larger = std::abs(largerBody.mass);

    // Only bounce bodies that are moving TOWARD each other. Overlapping bodies that are already separating (e.g.
    // fragments born overlapping) must keep separating -- bouncing them would flip them back together every tick, so
    // they'd stick.
    double closing_speed = v_smaller_n - v_larger_n; // > 0: the gap along the normal is shrinking
    if (closing_speed > 0.0)
    {
        double v_smaller_n_new =
            (v_smaller_n * (m_smaller - m_larger) + 2 * m_larger * v_larger_n) / (m_smaller + m_larger);
        double v_larger_n_new =
            (v_larger_n * (m_larger - m_smaller) + 2 * m_smaller * v_smaller_n) / (m_smaller + m_larger);

        smallerBody.velocity += normal_vector * (v_smaller_n_new - v_smaller_n) * ELASTIC_LOSS_FACTOR;
        largerBody.velocity += normal_vector * (v_larger_n_new - v_larger_n) * ELASTIC_LOSS_FACTOR;
    }

    double penetration = (smallerBody.radius + largerBody.radius) - distance;
    if (penetration > 0.0)
    {
        constexpr double percent = 0.8;
        constexpr double slop = 0.01;
        double correction_magnitude = std::max(penetration - slop, 0.0) * percent;
        DynamoEngine::Vector2D correction = normal_vector * correction_magnitude;

        if (smallerBody.isForceStatic && !largerBody.isForceStatic)
        {
            largerBody.position += correction;
            largerBody.previousPosition = largerBody.position;
        }
        else if (!smallerBody.isForceStatic && largerBody.isForceStatic)
        {
            smallerBody.position -= correction;
            smallerBody.previousPosition = smallerBody.position;
        }
        else
        {
            // Share the push by inverse mass (the lighter body moves more): ratios, so they work at any scale.
            // That only makes sense when both masses have the same sign; with one negative mass (see REWORK) the
            // shares can blow up, so split the push evenly instead.
            double smaller_share = 0.5;
            double larger_share = 0.5;
            if (smallerBody.invMass * largerBody.invMass > 0.0) // same sign
            {
                double total_inv_mass = smallerBody.invMass + largerBody.invMass;
                smaller_share = smallerBody.invMass / total_inv_mass;
                larger_share = largerBody.invMass / total_inv_mass;
            }
            smallerBody.position -= correction * smaller_share;
            largerBody.position += correction * larger_share;
            largerBody.previousPosition = largerBody.position;
            smallerBody.previousPosition = smallerBody.position;
        }
    }
}

void PhysicsSystem::handleAccretion(GravitationalBodyPair& grav_body_pair)
{
    GravitationalBody& heavier = *grav_body_pair.heavierBody;
    GravitationalBody& lighter = *grav_body_pair.lighterBody;

    if (heavier.isParticle)
    {
        return;
    }

    // Opposite signs subtract. The merge rule needs a 10x size ratio, so at least 90% of the heavier mass is always
    // left: new_mass never gets near 0 (dividing the momentum by it below stays safe) and never flips sign.
    double new_mass = heavier.mass + lighter.mass;
    // Force-static bodies count as infinitely heavy wherever momentum is exchanged (GravitationalBody::isForceStatic)
    if (heavier.isForceStatic)
    {
        // The absorber is "infinitely heavy": its velocity can't change, the absorbed momentum vanishes into it
    }
    else if (lighter.isForceStatic)
    {
        // The absorber now contains an "infinite" mass: it takes the static body's velocity and becomes static itself
        heavier.velocity = lighter.velocity;
        heavier.isForceStatic = true;
    }
    else
    {
        heavier.velocity = (heavier.velocity * heavier.mass + lighter.velocity * lighter.mass) / new_mass;
    }
    heavier.radius *= DynamoEngine::deterministicCbrt(new_mass / heavier.mass); // radius grows as the cube root of mass
    heavier.mass = new_mass;
    heavier.invMass = 1.0 / heavier.mass;

    // Disabled because disabling promotion //TODO: Review?
    // if (!heavier.isMacro && heavier.radius >= PARTICLE_PROMOTION_RADIUS_THRESHOLD)
    // {
    //     heavier.isMarkedForPromotion = true;
    // }

    lighter.isMarkedForDeletion = true;
}

// TODO: Prune? currently uncalled, see commented invocation in updateSystemFrame
void PhysicsSystem::promoteOversizedParticles(UniverseState& universe_state)
{
    auto& particles = universe_state.particles;
    auto& macro_bodies = universe_state.macro_bodies;

    for (auto& particle : particles)
    {
        if (!particle.isMarkedForPromotion)
        {
            continue;
        }

        universe_state.newest_spawn_id++;

        GravitationalBody promoted;
        promoted.position = particle.position;
        promoted.previousPosition = particle.previousPosition;
        promoted.velocity = particle.velocity;
        promoted.netForce = particle.netForce;
        promoted.prevForce = particle.prevForce;
        promoted.mass = particle.mass;
        promoted.invMass = particle.invMass;
        promoted.radius = particle.radius;
        promoted.macroIdentifier = universe_state.newest_spawn_id;

        promoted.isMacro = true;
        promoted.isAccretable = true;
        promoted.isCollidable = true;
        promoted.isShatterable = true;
        promoted.isPlanet = true;

        macro_bodies.push_back(promoted);

        // Defer actual removal to the existing end-of-tick particle cleanup rather than
        // erasing here, mid-iteration over the same particles vector.
        particle.isMarkedForDeletion = true;
    }
}

void PhysicsSystem::substituteWithParticles(GravitationalBody& original_body,
                                            std::vector<GravitationalBody>& fragments_out,
                                            uint32_t target_fragment_count, uint64_t random_key)
{
    size_t first_new_fragment = fragments_out.size();
    layOutAndPopulateFragments(original_body, target_fragment_count, fragments_out);

    // Each fragment keeps roughly the original's velocity, randomly scaled so the pieces spread out. "Random" means
    // hashed from the body's key and the fragment's number: the same every run.
    uint64_t speed_key = DynamoEngine::combineHash(random_key, RANDOM_PURPOSE_FRAGMENT_SPEED);
    for (size_t i = first_new_fragment; i < fragments_out.size(); ++i)
    {
        uint64_t fragment_number = i - first_new_fragment;
        double speed_scale = DynamoEngine::hashToRange(DynamoEngine::combineHash(speed_key, fragment_number), 0.8, 1.1);
        fragments_out[i].velocity = original_body.velocity * speed_scale;
    }

    original_body.isMarkedForDeletion = true;
}

void PhysicsSystem::substituteWithParticlesFromImpact(GravitationalBody& original_body,
                                                      std::vector<GravitationalBody>& fragments_out,
                                                      uint32_t target_fragment_count,
                                                      const DynamoEngine::Vector2D& impactPoint, uint64_t random_key)
{
    const double R = original_body.radius;
    const double original_mass = original_body.mass;

    size_t start_index = fragments_out.size(); // this body's fragments are the ones added after this point

    substituteWithParticles(original_body, fragments_out, target_fragment_count, random_key);

    // Grow radius with distance from the impact point: near-impact fragments stay small and
    // pulverized, far-side fragments stay large and coherent -- keyed off the actual contact
    // point (not the body's own center), so a graze reads differently from a direct hit.
    double total_weighted_volume = 0.0;
    uint64_t size_key = DynamoEngine::combineHash(random_key, RANDOM_PURPOSE_FRAGMENT_SIZE);
    for (size_t i = start_index; i < fragments_out.size(); ++i)
    {
        GravitationalBody& p = fragments_out[i];
        double distance = (p.position - impactPoint).magnitude();
        double t = std::clamp(distance / (2.0 * R), 0.0, 1.0);
        uint64_t fragment_number = i - start_index;
        double random_growth = DynamoEngine::hashToRange(DynamoEngine::combineHash(size_key, fragment_number), 0.5, 0.9);
        double size_multiplier = 1.0 + t * (IMPACT_SKEW_GROWTH_FACTOR - 1.0) * random_growth;
        p.radius *= size_multiplier;
        total_weighted_volume += p.radius * p.radius * p.radius;
    }

    // Re-normalize mass so the biased fragments still sum to originalMass, weighted by volume
    // (radius^3) to stay consistent with handleAccretion's mass<->radius law.
    for (size_t i = start_index; i < fragments_out.size(); ++i)
    {
        GravitationalBody& p = fragments_out[i];
        p.mass = original_mass * (p.radius * p.radius * p.radius) / total_weighted_volume;
        p.invMass = 1.0 / p.mass;
    }
}

// --------- GRAVITY --------- //

void PhysicsSystem::recomputeAllForces(UniverseState& universe_state)
{
    addGravityBetweenGravBodies(universe_state);
    addGravityBetweenShipAndBodies(universe_state);
    // Update other forces
}

void PhysicsSystem::addGravityBetweenGravBodies(UniverseState& universe_state)
{
    std::vector<GravitationalBody>& macro_bodies = universe_state.macro_bodies;
    std::vector<GravitationalBody>& particles = universe_state.particles;
    size_t num_macro_bodies = macro_bodies.size();
    size_t num_particles = particles.size();

    // Macro-Macro gravity
    {
        DynamoEngine::ProfileScope timer(m_profiler, m_profile.planet_gravity, m_profiling_real_tick);
        for (size_t i = 0; i < num_macro_bodies; i++)
        {
            for (size_t j = i + 1; j < num_macro_bodies; ++j)
            {
                addGravityBetweenGravBodyPair(macro_bodies[i], macro_bodies[j]);
            }
        }
    }

    // Particle-Macro gravity. Particles still don't gravitate with each other -- that stays
    // an explicit simplification (would be O(n^2) or need Barnes-Hut, both out of scope here).
    {
        DynamoEngine::ProfileScope timer(m_profiler, m_profile.debris_planet_gravity, m_profiling_real_tick);
        for (size_t i = 0; i < num_particles; ++i)
        {
            for (size_t j = 0; j < num_macro_bodies; ++j)
            {
                addGravityBetweenGravBodyPair(particles[i], macro_bodies[j]);
            }
        }
    }

    // Particle-Particle gravity.
    // TODO: Review? -- Works now with elastic collision fix, need to majorly optimize
    // {
    //     DynamoEngine::ProfileScope timer(m_profiler, m_profile.debris_debris_gravity, m_profiling_real_tick);
    //     for (size_t i = 0; i < num_particles; ++i)
    //     {
    //         for (size_t j = i + 1; j < num_particles; ++j)
    //         {
    //             addGravityBetweenGravBodyPair(particles[i], particles[j]);
    //         }
    //     }
    // }
}
void PhysicsSystem::addGravityBetweenShipAndBodies(UniverseState& universe_state)
{
    DynamoEngine::ProfileScope timer(m_profiler, m_profile.ship_gravity, m_profiling_real_tick);
    // Two-way gravity between the ship and every planet (Newton's third law): the same pull acts on both, in
    // opposite directions. The ship is so much lighter than any planet that planets barely notice it.
    if (!universe_state.has_ship)
    {
        return;
    }
    Starship& ship = universe_state.ship;
    DynamoEngine::Vector2D ship_center = ship.center();
    DynamoEngine::Vector2D total_force_on_ship = {0.0, 0.0};

    for (GravitationalBody& body : universe_state.macro_bodies)
    {
        if (firstWithinEpsilonOfSecond(body.mass, 0.0))
        {
            continue;
        }

        // Same softening rule as addGravityBetweenGravBodyPair (half the two sizes), so the pull stays finite up close
        double softening = (ship.halfSize() + body.radius) / 2.0;
        DynamoEngine::Vector2D direction = body.position - ship_center;
        double softened_distance = std::sqrt(direction.squareMagnitude() + softening * softening);
        double softened_distance_cubed = softened_distance * softened_distance * softened_distance;

        // The pull on the ship, toward the planet (direction isn't normalized: its length r is folded into r^3)
        DynamoEngine::Vector2D force =
            direction * (GRAVITATIONAL_CONSTANT * ship.mass() * body.mass / softened_distance_cubed);
        total_force_on_ship += force;
        if (!body.isForceStatic)
        {
            body.netForce -= force; // the equal and opposite pull on the planet
        }
    }

    ship.setGravityForce(total_force_on_ship);
}

void PhysicsSystem::addGravityBetweenGravBodyPair(GravitationalBody& firstBody, GravitationalBody& secondBody)
{
    if (firstWithinEpsilonOfSecond(firstBody.mass, 0.0) || firstWithinEpsilonOfSecond(secondBody.mass, 0.0))
    {
        return;
    }
    // TODO: Review?
    // if (firstWithinEpsilonOfSecond(firstBody.radius, 0.0) || firstWithinEpsilonOfSecond(secondBody.radius, 0.0))
    // {
    //     return;
    // }
    const double G = GRAVITATIONAL_CONSTANT;

    double softening_constant = (firstBody.radius + secondBody.radius) / 2.0;
    double epsilon_squared = softening_constant * softening_constant;

    // 1. Calculate direction Vector
    DynamoEngine::Vector2D direction_vector = secondBody.position - firstBody.position;

    // 2. Calculate Distance Squared
    double r_squared = direction_vector.squareMagnitude();

    // 3. Calculate the Denominator Term (r^2 + epsilon^2)^(3/2)
    double denominator_1 = sqrt(r_squared + epsilon_squared);
    double denominator_3 = denominator_1 * denominator_1 * denominator_1;

    // 4. Calculate Coefficient (F = direction vector * G * m1 * m2 /
    // Denominator)
    double coefficient = (G * firstBody.mass * secondBody.mass) / denominator_3;

    // Force vector is C * direction_vector (r)
    DynamoEngine::Vector2D force = direction_vector * coefficient;

    // 6. Apply Forces (Newton's Third Law)
    if (!firstBody.isForceStatic)
        firstBody.netForce += force;

    if (!secondBody.isForceStatic)
        secondBody.netForce -= force;
}

// --------- GRAVITATIONAL BODY CREATION --------- //

// Fills in a new planet (or a cluster's template body) from where it was aimed and the spawn panel's settings
static inline void populateSpawnedBodyProperties(GravitationalBody& grav_body, const SpawnProperties& spawn_aim,
                                                 const SpawnSettings& spawn_settings)
{
    grav_body.mass = spawn_settings.mass;
    grav_body.invMass = 1 / grav_body.mass;
    grav_body.radius = spawn_settings.radius;
    grav_body.position = spawn_aim.world_position;
    grav_body.previousPosition = grav_body.position;
    grav_body.velocity = spawn_aim.world_velocity; // 0 unless it was aimed with Shift

    // Type flag
    grav_body.isMacro = true;

    // Property Flags
    grav_body.isAccretable = spawn_settings.is_accretable;
    grav_body.isBounce = spawn_settings.is_bounce;
    grav_body.isCollidable = !spawn_settings.is_ghost;
    grav_body.isForceStatic = spawn_settings.is_force_static;
    grav_body.isFragment = false;   // not a child of a collision
    grav_body.isMacroGhost = false; // no spawn setting for it (physics still honors the flag)
    grav_body.isPreview = false;
    grav_body.isShatterable = spawn_settings.is_shatterable;
    grav_body.isTransient = false;

    // Visual Identifier Flags
    grav_body.isDust = false;
    grav_body.isGas = false;
    grav_body.isGravStar = false;
    grav_body.isMoon = false;
    grav_body.isPlanet = true; // Default for now until I implement procedural texture generation.
}

bool PhysicsSystem::makeSpawnPreviewBody(const SpawnProperties& spawn_aim, const SpawnSettings& spawn_settings,
                                         GravitationalBody& body_out) const
{
    // Same refusals as createMacroBody (the particle budget aside: a refused spawn simply won't happen)
    if (spawn_settings.radius <= 1.0 || firstWithinEpsilonOfSecond(spawn_settings.mass, 0.0))
    {
        return false;
    }
    body_out = GravitationalBody();
    populateSpawnedBodyProperties(body_out, spawn_aim, spawn_settings);
    return true;
}

void PhysicsSystem::createMacroBody(UniverseState& universe_state, const SpawnProperties& spawn_aim,
                                    const SpawnSettings& spawn_settings)
{
    std::vector<GravitationalBody>& macro_bodies = universe_state.macro_bodies;
    if (spawn_settings.radius <= 1.0)
    {
        return;
    }
    if (firstWithinEpsilonOfSecond(spawn_settings.mass, 0.0))
    {
        return;
    }
    // Budget: a new macro body might shatter later, so it needs room for all of its fragments up front
    if (!hasRoomForParticles(universe_state, DEFAULT_FRAGMENT_COUNT))
    {
        return;
    }
    universe_state.newest_spawn_id++;
    int new_macro_body_id = universe_state.newest_spawn_id;

    GravitationalBody macro_body;
    macro_body.macroIdentifier = new_macro_body_id;

    populateSpawnedBodyProperties(macro_body, spawn_aim, spawn_settings);

    // Now with populated flags, nudge particles out?
    macro_bodies.push_back(macro_body);
}

void PhysicsSystem::createParticleCluster(UniverseState& universe_state, const SpawnProperties& spawn_aim,
                                          const SpawnSettings& spawn_settings)
{
    std::vector<GravitationalBody>& particles = universe_state.particles;
    if (spawn_settings.radius <= 1.0)
    {
        return;
    }
    if (firstWithinEpsilonOfSecond(spawn_settings.mass, 0.0))
    {
        return;
    }

    GravitationalBody macro_body;

    populateSpawnedBodyProperties(macro_body, spawn_aim, spawn_settings);

    // Load limit: refuse the whole cluster rather than add part of it (needs the radius, so it comes after populate)
    if (!hasRoomForParticles(universe_state, survivableFragmentCount(macro_body, DEFAULT_FRAGMENT_COUNT)))
    {
        return;
    }
    // A spawn ID like a planet's: its particles inherit it, so Delete can undo the whole cluster at once
    // (and each cluster gets its own shader seed)
    universe_state.newest_spawn_id++;
    macro_body.macroIdentifier = universe_state.newest_spawn_id;
    size_t first_new_particle = particles.size();
    // No loop is running: add them directly
    substituteWithParticles(macro_body, particles, DEFAULT_FRAGMENT_COUNT, bodyRandomKey(universe_state, macro_body));
    for (size_t i = first_new_particle; i < particles.size(); ++i)
    {
        particles[i].isAccretable = spawn_settings.is_accretable;
        particles[i].isBounce = spawn_settings.is_bounce;
        particles[i].isCollidable = !spawn_settings.is_ghost;
    }
}
// --------- CLEANUP --------- //

void PhysicsSystem::cleanupParticles(UniverseState& universe_state)
{
    auto& particles = universe_state.particles;

    // 1. Use std::remove_if to move all elements marked for deletion
    //    to the end of the vector. It returns an iterator to the new
    //    end.
    auto new_end = std::remove_if(
        particles.begin(), particles.end(), [](const GravitationalBody& p)
        { return p.isMarkedForDeletion || firstWithinEpsilonOfSecond(p.mass, 0.0) || p.radius < MIN_PARTICLE_RADIUS; });

    // 2. Use vector::erase to destroy the elements in the range
    // [new_end, particles.end())
    //    This efficiently shrinks the vector to the correct size.
    particles.erase(new_end, particles.end());
}

void PhysicsSystem::cleanupMacroBodies(UniverseState& universe_state)
{
    auto& particles = universe_state.macro_bodies;

    // 1. Use std::remove_if to move all elements marked for deletion
    //    to the end of the vector. It returns an iterator to the new
    //    end.
    auto new_end = std::remove_if(particles.begin(), particles.end(), [](const GravitationalBody& b)
                                  { return b.isMarkedForDeletion || firstWithinEpsilonOfSecond(b.mass, 0.0); });

    // 2. Use vector::erase to destroy the elements in the range
    // [new_end, particles.end())
    //    This efficiently shrinks the vector to the correct size.
    particles.erase(new_end, particles.end());
}

void PhysicsSystem::adjustThrottleWhilePlanning(UniverseState& universe_state, const ShipControls& ship_controls,
                                                float frame_seconds)
{
    if (universe_state.has_ship)
    {
        universe_state.ship.adjustThrottle(ship_controls.thrust, frame_seconds);
    }
}

void PhysicsSystem::applyShipInputPhysics(UniverseState& universe_state, const ShipControls& ship_controls)
{
    if (!universe_state.has_ship)
    {
        return;
    }
    universe_state.ship.applyRotation(ship_controls);
    universe_state.ship.applyThrust(ship_controls);
}

bool PhysicsSystem::hasRoomForParticles(const UniverseState& universe_state, size_t new_particle_count) const
{
    size_t potential_particles = potentialParticleCount(universe_state);
    return (potential_particles + new_particle_count) <= MAX_LIVE_PARTICLES;
}

size_t PhysicsSystem::potentialParticleCount(const UniverseState& universe_state) const
{
    // A flat DEFAULT_FRAGMENT_COUNT per macro, not survivableFragmentCount: a macro that grows by accretion would
    // raise its exact count after it was budgeted, but it can never raise the flat one
    return liveParticleCount(universe_state) + universe_state.macro_bodies.size() * DEFAULT_FRAGMENT_COUNT;
}

// True if the body's position and velocity are real numbers (not infinity, not NaN = "not a number")
static bool hasFiniteMotion(const GravitationalBody& body)
{
    return std::isfinite(body.position.x_val) && std::isfinite(body.position.y_val) &&
           std::isfinite(body.velocity.x_val) && std::isfinite(body.velocity.y_val);
}

// Safety net: marks every body whose motion has stopped being a real number for deletion, and returns how many.
// Such a body is a bug somewhere upstream (a division by ~0), and it spreads: anything that touches or pulls on it
// turns NaN too. It also wrecks the frame rate: the grid files every NaN position under one cell, and NaN fails the
// "too far apart to touch" test (any comparison with NaN is false), so every NaN pair counts as a contact.
// Deleting them costs a few bodies instead.
static size_t markNonFiniteBodiesForDeletion(std::vector<GravitationalBody>& bodies)
{
    size_t marked = 0;
    for (GravitationalBody& body : bodies)
    {
        if (!body.isMarkedForDeletion && !hasFiniteMotion(body))
        {
            body.isMarkedForDeletion = true;
            marked += 1;
        }
    }
    return marked;
}

void PhysicsSystem::cleanupGravBodies(UniverseState& universe_state)
{
    DynamoEngine::ProfileScope timer(m_profiler, m_profile.cleanup, m_profiling_real_tick);
    size_t non_finite_count = markNonFiniteBodiesForDeletion(universe_state.macro_bodies) +
                              markNonFiniteBodiesForDeletion(universe_state.particles);
    // Said out loud so the bug behind it gets found (real ticks only: predictions run on copies of the same bodies)
    if (non_finite_count > 0 && m_profiling_real_tick)
    {
        std::cerr << "PhysicsSystem: deleted " << non_finite_count
                  << " bodies whose position or velocity was infinite or NaN (a physics bug upstream)\n";
    }
    cleanupMacroBodies(universe_state);
    cleanupParticles(universe_state);
}