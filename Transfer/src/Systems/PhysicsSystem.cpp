// File: Transfer/src/Systems/PhysicsSystem.cpp

#include "PhysicsSystem.hpp"
#include "Utilities/Constants/PhysicsConstants.hpp"

PhysicsSystem::PhysicsSystem()
{
    // Initialize physics system variables if needed
}

PhysicsSystem::~PhysicsSystem() {}

void PhysicsSystem::UpdateSystemFrame(GameState& game_state, UIState& ui_state)
{
    // Mental Model of System

    // Game updates gravitational body instantiations outside of Physics loop because the instantiations should be
    // rendered immediately (even with stopped time)

    // First, handle collisions, then integrate, then update forces, then integrate again, then finally cleanup

    handleCollisions(game_state);

    updatePlayerPhysics(game_state, ui_state);

    // promoteOversizedParticles(game_state); //TODO: Review?

    // Update forces (grav, later will add electromagnetic)
    integrateForwardsVelocityVerletPhase1(game_state);
    updateAllForces(game_state);
    integrateForwardsVelocityVerletPhase2(game_state);

    cleanupMacroBodies(game_state);
    cleanupParticles(game_state);

    // Verify calculation doesn't cook us
    // calculateTotalEnergy(game_state); //TODO: Review?
}

void PhysicsSystem::CleanUp()
{
    // Any necessary cleanup code for the physics system
}

void PhysicsSystem::UpdateGravBodyInstantiations(GameState& game_state, UIState& ui_state)
{
    DEPRECATED_InputState& input_state = ui_state.getMutableDEPRECATED_InputState();
    if (!input_state.UIInputConsumed)
    {
        if (input_state.isCreatingMacro)
        {
            createMacroBody(
                game_state, input_state,
                ui_state.spawnSettings()); // can inline replace with other create* methods for test. isCreatingMacro
                                           // should eventually be unique to just creating planets, etc.
            input_state.resetTransientFlags();
        }
        else if (input_state.isCreatingParticleCluster)
        {
            createParticleCluster(game_state, input_state, ui_state.spawnSettings());
            input_state.resetTransientFlags();
        }
    }
    else
    {
        // run limiters?
    }
    // Undo the newest spawn (plain Delete)
    if (input_state.removeNewestSpawn)
    {
        removeNewestSpawn(game_state);
        input_state.removeNewestSpawn = false;
    }
    // Check if all Gravitational Bodies are supposed to be wiped
    if (input_state.clearAll)
    {
        game_state.getMacroBodiesMutable().clear();
        game_state.getParticlesMutable().clear();
        input_state.clearAll = false;
    }
}
void PhysicsSystem::removeNewestSpawn(GameState& game_state)
{
    std::vector<GravitationalBody>& macro_bodies = game_state.getMacroBodiesMutable();
    std::vector<GravitationalBody>& particles = game_state.getParticlesMutable();

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
    double distance = r_vector.magnitude();
    DynamoEngine::Vector2D unit_normal_vector =
        (distance > 1e-8) ? (r_vector / distance) : DynamoEngine::Vector2D(1.0, 0.0);

    DynamoEngine::Vector2D relative_velocity_vector = b.velocity - a.velocity;
    double normal_speed = relative_velocity_vector.dot(unit_normal_vector);
    double abs_normal_speed = std::abs(normal_speed);
    bool should_collide = (distance < b.radius + a.radius);
    double approaching_speed = -normal_speed; // positive if bodies are moving towards each other
    bool should_blow_up = (approaching_speed >= MIN_SHATTER_SPEED);
    return {distance,       unit_normal_vector, relative_velocity_vector, normal_speed, abs_normal_speed,
            should_collide, should_blow_up};
}

static inline GravitationalBodyPair pickMassPair(GravitationalBody& a, GravitationalBody& b)
{
    if (abs(a.mass) == abs(b.mass))
        return {&a, &b, 1.0, true};
    if (abs(a.mass) > abs(b.mass))
        return {&a, &b, abs(a.mass / b.mass), false};
    return {&b, &a, abs(b.mass / a.mass), false};
}

void PhysicsSystem::handleCollisions(GameState& game_state)
{
    // Shatters during these passes put their fragments in m_pending_fragments, not in particles
    handleMacroMacroCollisions(game_state);
    handleMacroParticleCollisions(game_state);
    handleParticleParticleCollisions(game_state);
    handleShipCollisions(game_state);
    handleParticleParticleCollisions(game_state);

    // All loops over particles are finished: now it's safe to grow the vector
    std::vector<GravitationalBody>& particles = game_state.getParticlesMutable();
    particles.insert(particles.end(), m_pending_fragments.begin(), m_pending_fragments.end());
    m_pending_fragments.clear(); // clear() keeps the memory, so next tick's shatters don't allocate again
}

size_t PhysicsSystem::liveParticleCount(const GameState& game_state) const
{
    return game_state.getParticles().size() + m_pending_fragments.size();
}

void PhysicsSystem::handleMacroMacroCollisions(GameState& game_state)
{
    std::vector<GravitationalBody>& macro_body_list = game_state.getMacroBodiesMutable();
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
            handleDynamicCollision(grav_body_pair, collision_info, game_state);
            // first_body may have just shattered or been absorbed: it's gone, so it can't hit anything else this tick
            if (first_body.isMarkedForDeletion)
            {
                break;
            }
        }
    }
}

void PhysicsSystem::handleMacroParticleCollisions(GameState& game_state)
{
    std::vector<GravitationalBody>& particles = game_state.getParticlesMutable();
    std::vector<GravitationalBody>& macro_bodies = game_state.getMacroBodiesMutable();

    for (auto& particle : particles)
    {
        if (!particle.isCollidable || particle.isMarkedForDeletion)
        {
            continue;
        }

        for (auto& macro_body : macro_bodies)
        {
            if (!macro_body.isCollidable || macro_body.isMacroGhost || macro_body.isMarkedForDeletion)
            {
                continue;
            }

            CollisionInfo collision_info = getCollisionInfo(particle, macro_body);
            if (!collision_info.shouldCollide)
            {
                continue;
            }

            GravitationalBodyPair grav_body_pair = pickMassPair(particle, macro_body);
            handleDynamicCollision(grav_body_pair, collision_info, game_state);
            // The particle may have just been absorbed: it's gone, so it can't hit anything else this tick
            if (particle.isMarkedForDeletion)
            {
                break;
            }
        }
    }
}

void PhysicsSystem::handleParticleParticleCollisions(GameState& game_state)
{
    std::vector<GravitationalBody>& particles = game_state.getParticlesMutable();

    particleGrid.build(particles);

    std::vector<size_t> candidates;
    size_t num_particles = particles.size();

    for (size_t i = 0; i < num_particles; ++i)
    {
        if (!particles[i].isCollidable || particles[i].isMarkedForDeletion)
        {
            continue;
        }

        particleGrid.queryCandidates(i, particles, candidates);
        for (size_t j : candidates)
        {
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
            handleDynamicCollision(grav_body_pair, collision_info, game_state);
            // particles[i] may have just been absorbed: it's gone, so it can't hit anything else this tick
            if (particles[i].isMarkedForDeletion)
            {
                break;
            }
        }
    }
}
void PhysicsSystem::handleShipCollisions(GameState& game_state)
{
    Starship& ship = game_state.getPlayerMutable().starship;
    double impact_energy_this_tick = 0.0;

    // Planets
    for (GravitationalBody& body : game_state.getMacroBodiesMutable())
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
    for (GravitationalBody& particle : game_state.getParticlesMutable())
    {
        if (!particle.isCollidable || particle.isMarkedForDeletion)
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
                                           GameState& game_state)
{
    GravitationalBody& heavier = *grav_body_pair.heavierBody;
    GravitationalBody& lighter = *grav_body_pair.lighterBody;

    if (heavier.isBounce || lighter.isBounce)
    {
        handleElasticCollisions(lighter, heavier);
        return;
    }

    if (collisionInfo.shouldBlowUp && lighter.isShatterable)
    {
        if (!lighter.isShatterable)
        {
            handleElasticCollisions(lighter, heavier);
            return;
        }

        DynamoEngine::Vector2D toward_lighter = (lighter.position - heavier.position).normalize();
        DynamoEngine::Vector2D impact_point = heavier.position + toward_lighter * heavier.radius;

        if (heavier.isShatterable && grav_body_pair.ratio <= MUTUAL_SHATTER_MASS_RATIO_THRESHOLD)
        {
            substituteWithParticlesFromImpact(heavier, m_pending_fragments, DEFAULT_FRAGMENT_COUNT, impact_point);
        }
        substituteWithParticlesFromImpact(lighter, m_pending_fragments, DEFAULT_FRAGMENT_COUNT, impact_point);
        return;
    }

    if (collisionInfo.absNormalSpeed >= MAX_ACCRETION_COLLISION_SPEED)
    {
        // Too fast to cleanly merge, not fast enough to shatter: always bounce.
        handleElasticCollisions(lighter, heavier);
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
            substituteWithParticlesFromImpact(lighter, m_pending_fragments, DEFAULT_FRAGMENT_COUNT / 3, contact_point);
        }
        else if (heavier.isParticle)
        {
            // Particles can't accrete, so we just handle an elastic collision.
            handleElasticCollisions(lighter, heavier);
        }
        else
        {
            handleAccretion(grav_body_pair);
        }
    }
    else
    {
        handleElasticCollisions(lighter, heavier);
    }
}

void PhysicsSystem::handleElasticCollisions(GravitationalBody& smallerBody, GravitationalBody& largerBody)
{
    if (smallerBody.isForceStatic && largerBody.isForceStatic)
    {
        // Two "infinitely heavy" bodies: neither outweighs the other, so they share the bounce like two bodies of EQUAL
        // mass. Each takes half the push-apart and half the velocity change. Both stay static.
        DynamoEngine::Vector2D offset = largerBody.position - smallerBody.position;
        double distance = offset.magnitude();
        if (firstWithinEpsilonOfSecond(distance, 0.0))
        {
            return; // exactly on top of each other: no direction to push in
        }
        DynamoEngine::Vector2D n = offset / distance; // points from the smaller body to the larger one

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

        DynamoEngine::Vector2D offset = dyn.position - stat.position;
        double distance = offset.magnitude();
        if (firstWithinEpsilonOfSecond(distance, 0.0))
        {
            return; // exactly on top of each other: no direction to push in
        }
        DynamoEngine::Vector2D n = offset / distance; // points from the static body to the other one

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

    DynamoEngine::Vector2D r_vector = largerBody.position - smallerBody.position;
    double distance = r_vector.magnitude();

    if (firstWithinEpsilonOfSecond(distance, 0.0))
    {
        return;
    }

    DynamoEngine::Vector2D normal_vector = r_vector / distance; // points from the smaller body to the larger one
    double v_smaller_n = smallerBody.velocity.dot(normal_vector);
    double v_larger_n = largerBody.velocity.dot(normal_vector);
    double m_smaller = smallerBody.mass;
    double m_larger = largerBody.mass;

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
            double totalInvMass = smallerBody.invMass + largerBody.invMass;
            if (totalInvMass > DynamoEngine::EPSILON)
            {
                double smaller_share = smallerBody.invMass / totalInvMass;
                double larger_share = largerBody.invMass / totalInvMass;
                smallerBody.position -= correction * smaller_share;
                largerBody.position += correction * larger_share;
            }
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
    heavier.radius *= pow(new_mass / heavier.mass, 1.0 / 3.0);
    heavier.mass = new_mass;
    heavier.invMass = 1.0 / heavier.mass;

    // Disabled because disabling promotion //TODO: Review?
    // if (!heavier.isMacro && heavier.radius >= PARTICLE_PROMOTION_RADIUS_THRESHOLD)
    // {
    //     heavier.isMarkedForPromotion = true;
    // }

    lighter.isMarkedForDeletion = true;
}

// TODO: Prune? currently uncalled, see commented invocation in UpdateSystemFrame
void PhysicsSystem::promoteOversizedParticles(GameState& game_state)
{
    auto& particles = game_state.getParticlesMutable();
    auto& macro_bodies = game_state.getMacroBodiesMutable();

    for (auto& particle : particles)
    {
        if (!particle.isMarkedForPromotion)
        {
            continue;
        }

        game_state.incrementMaxIDInstantiated();

        GravitationalBody promoted;
        promoted.position = particle.position;
        promoted.previousPosition = particle.previousPosition;
        promoted.velocity = particle.velocity;
        promoted.netForce = particle.netForce;
        promoted.prevForce = particle.prevForce;
        promoted.mass = particle.mass;
        promoted.invMass = particle.invMass;
        promoted.radius = particle.radius;
        promoted.macroIdentifier = game_state.getMaxIDInstantiated();

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
                                            std::vector<GravitationalBody>& fragments_out, uint32_t targetFragmentCount)
{
    // uint32_t num_particles = std::max<uint32_t>(1, targetFragmentCount);
    uint32_t num_particles = survivableFragmentCount(original_body, targetFragmentCount);

    const double R = original_body.radius;
    const DynamoEngine::Vector2D center = original_body.position;
    const double original_mass = original_body.mass;
    const DynamoEngine::Vector2D original_velocity = original_body.velocity;

    double density_factor = (DynamoEngine::PI * R * R) / num_particles;
    const double fragment_radius = OVERLAP_MARGIN * sqrt(density_factor / DynamoEngine::PI);
    const double particle_mass = original_mass / num_particles;

    for (uint32_t k = 0; k < num_particles; ++k)
    {
        double r_k = R * sqrt((k + 0.5) / num_particles);
        double theta_k = k * GOLDEN_ANGLE;
        DynamoEngine::Vector2D pos_k = center + DynamoEngine::Vector2D{r_k * cos(theta_k), r_k * sin(theta_k)};

        GravitationalBody p;
        p.mass = particle_mass;
        p.invMass = 1.0 / particle_mass;
        p.radius = fragment_radius;
        p.position = pos_k;
        p.previousPosition = pos_k;
        p.isFragment = true;
        p.isAccretable = true;
        p.isCollidable = true;
        p.macroIdentifier = original_body.macroIdentifier;
        p.isParticle = true;
        p.velocity = original_velocity * randomDouble(0.8, 1.1);

        fragments_out.push_back(p);
    }

    original_body.isMarkedForDeletion = true;
}

void PhysicsSystem::substituteWithParticlesFromImpact(GravitationalBody& original_body,
                                                      std::vector<GravitationalBody>& fragments_out,
                                                      uint32_t targetFragmentCount,
                                                      const DynamoEngine::Vector2D& impactPoint)
{
    const double R = original_body.radius;
    const double original_mass = original_body.mass;

    size_t start_index = fragments_out.size(); // this body's fragments are the ones added after this point

    substituteWithParticles(original_body, fragments_out, targetFragmentCount);

    // Grow radius with distance from the impact point: near-impact fragments stay small and
    // pulverized, far-side fragments stay large and coherent -- keyed off the actual contact
    // point (not the body's own center), so a graze reads differently from a direct hit.
    double total_weighted_volume = 0.0;
    for (size_t i = start_index; i < fragments_out.size(); ++i)
    {
        GravitationalBody& p = fragments_out[i];
        double distance = (p.position - impactPoint).magnitude();
        double t = std::clamp(distance / (2.0 * R), 0.0, 1.0);
        double size_multiplier = 1.0 + t * (IMPACT_SKEW_GROWTH_FACTOR - 1.0) * randomDouble(0.5, 0.9);
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

void PhysicsSystem::updateAllForces(GameState& game_state)
{
    updateGravityForSystem(game_state);
    updateShipGravity(game_state);
    // Update other forces
}

void PhysicsSystem::updateGravityForSystem(GameState& game_state)
{
    std::vector<GravitationalBody>& macro_bodies = game_state.getMacroBodiesMutable();
    std::vector<GravitationalBody>& particles = game_state.getParticlesMutable();
    size_t num_macro_bodies = macro_bodies.size();
    size_t num_particles = particles.size();

    // Macro-Macro gravity
    for (size_t i = 0; i < num_macro_bodies; i++)
    {
        for (size_t j = i + 1; j < num_macro_bodies; ++j)
        {
            calculateGravity(macro_bodies[i], macro_bodies[j]);
        }
    }

    // Particle-Macro gravity. Particles still don't gravitate with each other -- that stays
    // an explicit simplification (would be O(n^2) or need Barnes-Hut, both out of scope here).
    for (size_t i = 0; i < num_particles; ++i)
    {
        for (size_t j = 0; j < num_macro_bodies; ++j)
        {
            calculateGravity(particles[i], macro_bodies[j]);
        }
    }

    // TODO: Review? -- Works now with elastic collision fix, need to majorly optimize
    // for (size_t i = 0; i < num_particles; ++i)
    // {
    //     for (size_t j = i + 1; j < num_particles; ++j)
    //     {
    //         calculateGravity(particles[i], particles[j]);
    //     }
    // }
}
void PhysicsSystem::updateShipGravity(GameState& game_state)
{
    // Two-way gravity between the ship and every planet (Newton's third law): the same pull acts on both, in
    // opposite directions. The ship is so much lighter than any planet that planets barely notice it.
    Starship& ship = game_state.getPlayerMutable().starship;
    DynamoEngine::Vector2D ship_center = ship.center();
    DynamoEngine::Vector2D total_force_on_ship = {0.0, 0.0};

    for (GravitationalBody& body : game_state.getMacroBodiesMutable())
    {
        if (firstWithinEpsilonOfSecond(body.mass, 0.0))
        {
            continue;
        }

        // Same softening rule as calculateGravity (half the two sizes), so the pull stays finite up close
        double softening = (ship.halfSize() + body.radius) / 2.0;
        DynamoEngine::Vector2D direction = body.position - ship_center;
        double softened_distance = std::sqrt(direction.squareMagnitude() + softening * softening);
        double softened_distance_cubed = softened_distance * softened_distance * softened_distance;

        // The pull on the ship, toward the planet (direction isn't normalised: its length r is folded into r^3)
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

void PhysicsSystem::calculateGravity(GravitationalBody& firstBody, GravitationalBody& secondBody)
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

// --------- INTEGRATION (VELOCITY VERLET) --------- //
// Position integration from previous frame force.

// Honestly check if we need the velocity integration calculation at all. Are we not updated that at end of
// integration phase 2 of prev frame?

void PhysicsSystem::integrateForwardsVelocityVerletPhase1(GameState& game_state)
{
    std::vector<GravitationalBody>& particles = game_state.getParticlesMutable();
    std::vector<GravitationalBody>& macro_bodies = game_state.getMacroBodiesMutable();

    for (auto& particle : particles)
    {
        applyVelocityVerletPhase1(particle);
    }
    for (auto& macro_body : macro_bodies)
    {
        applyVelocityVerletPhase1(macro_body);
    }
    game_state.getPlayerMutable().starship.applyVelocityVerletPhase1();
}

void PhysicsSystem::applyVelocityVerletPhase1(GravitationalBody& grav_body)
{
    grav_body.previousPosition = grav_body.position;
    grav_body.prevForce = grav_body.netForce;
    bool has_mass = !(firstWithinEpsilonOfSecond(grav_body.mass, 0.0));
    if (has_mass && !grav_body.isForceStatic)
    {
        // Calculate the acceleration from the previous frame's final force
        DynamoEngine::Vector2D acceleration = grav_body.netForce * grav_body.invMass;
        grav_body.velocity += acceleration * (PHYSICS_TIME_STEP / 2); // Half of a full integrated step
    }
    // Step the position
    grav_body.position += grav_body.velocity * PHYSICS_TIME_STEP;
    // Reset to force 0 for next frame
    grav_body.netForce = DynamoEngine::Vector2D(0.0, 0.0);
}

void PhysicsSystem::integrateForwardsVelocityVerletPhase2(GameState& game_state)
{
    std::vector<GravitationalBody>& particles = game_state.getParticlesMutable();
    std::vector<GravitationalBody>& macro_bodies = game_state.getMacroBodiesMutable();

    for (auto& particle : particles)
    {
        applyVelocityVerletPhase2(particle);
    }
    for (auto& macro_body : macro_bodies)
    {
        applyVelocityVerletPhase2(macro_body);
    }
    game_state.getPlayerMutable().starship.applyVelocityVerletPhase2();
}

void PhysicsSystem::applyVelocityVerletPhase2(GravitationalBody& grav_body)
{
    bool has_mass = !(firstWithinEpsilonOfSecond(grav_body.mass, 0.0));

    if (!has_mass || grav_body.isForceStatic)
    {
        return;
    }
    else
    {
        DynamoEngine::Vector2D acceleration = grav_body.netForce * grav_body.invMass;
        grav_body.velocity += acceleration * (PHYSICS_TIME_STEP / 2.0); // Other half of full integrated step
    }
}

// --------- GRAVITATIONAL BODY CREATION --------- //

static inline void populateGravBodyPropertiesFromDEPRECATED_InputState(GravitationalBody& grav_body,
                                                                       GameState& game_state,
                                                                       DEPRECATED_InputState& input_state,
                                                                       const SpawnSettings& spawn_settings)
{
    // Default sets, position may be overwrriten if isCreatingWithInitialVelocity set to true
    grav_body.mass = input_state.selectedMass;
    grav_body.invMass = 1 / grav_body.mass;
    grav_body.radius = input_state.selectedRadius;
    grav_body.position = ScreenToWorldCoordinates(input_state.mouseCurrPosition, game_state.getCameraState());
    grav_body.previousPosition = grav_body.position;

    // Flags (will be moved to control from within the input_state)
    // Type flag
    grav_body.isMacro = true;

    // Property Flags
    grav_body.isAccretable = spawn_settings.is_accretable;
    grav_body.isBounce = spawn_settings.is_bounce;
    grav_body.isCollidable = true;
    grav_body.isForceStatic = spawn_settings.is_force_static;
    grav_body.isFragment = false; // not a child of a collision
    grav_body.isMacroGhost = false;
    grav_body.isPreview = false;
    grav_body.isShatterable = spawn_settings.is_shatterable;
    grav_body.isTransient = false;

    // Visual Identifier Flags
    grav_body.isDust = false;
    grav_body.isGas = false;
    grav_body.isGravStar = false;
    grav_body.isMoon = false;
    grav_body.isPlanet = true; // Default for now until I implement procedural texture generation.

    if (input_state.isCreatingWithInitialVelocity)
    {
        grav_body.position = ScreenToWorldCoordinates(input_state.mouseDragStartPosition, game_state.getCameraState());
        grav_body.previousPosition = grav_body.position;
        grav_body.velocity =
            (input_state.mouseCurrPosition - input_state.mouseDragStartPosition) / game_state.getCameraState().zoom;
    }
}

void PhysicsSystem::createMacroBody(GameState& game_state, DEPRECATED_InputState& input_state,
                                    const SpawnSettings& spawn_settings)
{
    std::vector<GravitationalBody>& macro_bodies = game_state.getMacroBodiesMutable();
    if (input_state.selectedRadius <= 1.0)
    {
        return;
    }
    if (firstWithinEpsilonOfSecond((input_state.selectedMass), 0.0))
    {
        return;
    }
    // Budget: a new macro body might shatter later, so it needs room for all of its fragments up front
    if (!hasRoomForParticles(game_state, DEFAULT_FRAGMENT_COUNT))
    {
        return;
    }
    game_state.incrementMaxIDInstantiated();
    int new_macro_body_id = game_state.getMaxIDInstantiated();

    GravitationalBody macro_body;
    macro_body.macroIdentifier = new_macro_body_id;

    // Pass flags from input_state as possible.
    populateGravBodyPropertiesFromDEPRECATED_InputState(macro_body, game_state, input_state, spawn_settings);

    // Now with populated flags, nudge particles out?
    macro_bodies.push_back(macro_body);
}

void PhysicsSystem::createParticleCluster(GameState& game_state, DEPRECATED_InputState& input_state,
                                          const SpawnSettings& spawn_settings)
{
    std::vector<GravitationalBody>& particles = game_state.getParticlesMutable();
    if (input_state.selectedRadius <= 1.0)
    {
        return;
    }
    if (firstWithinEpsilonOfSecond((input_state.selectedMass), 0.0))
    {
        return;
    }

    GravitationalBody macro_body;

    // Pass flags from input_state as possible.
    populateGravBodyPropertiesFromDEPRECATED_InputState(macro_body, game_state, input_state, spawn_settings);

    // Load limit: refuse the whole cluster rather than add part of it (needs the radius, so it comes after populate)
    if (!hasRoomForParticles(game_state, survivableFragmentCount(macro_body, DEFAULT_FRAGMENT_COUNT)))
    {
        return;
    }
    // A spawn ID like a planet's: its particles inherit it, so Delete can undo the whole cluster at once
    // (and each cluster gets its own shader seed)
    game_state.incrementMaxIDInstantiated();
    macro_body.macroIdentifier = game_state.getMaxIDInstantiated();
    size_t first_new_particle = particles.size();
    substituteWithParticles(macro_body, particles, DEFAULT_FRAGMENT_COUNT); // no loop is running: add them directly
    for (size_t i = first_new_particle; i < particles.size(); ++i)
    {
        particles[i].isAccretable = spawn_settings.is_accretable;
        particles[i].isBounce = spawn_settings.is_bounce;
    }
}

// --------- UTILITY --------- //

// TODO: Prune?
void PhysicsSystem::calculateTotalEnergy(GameState& game_state)
{
    auto& macro_bodies = game_state.getMacroBodies();
    size_t num_macro_bodies = macro_bodies.size();

    double total_energy = 0.0;
    for (int i = 0; i < num_macro_bodies; ++i)
    {
        total_energy += macro_bodies[i].mass * macro_bodies[i].velocity.squareMagnitude() / 2.0;
    }

    for (size_t i = 0; i < num_macro_bodies; i++)
    {
        for (size_t j = i + 1; j < num_macro_bodies; j++)
        {
            double epsilon = (macro_bodies[i].radius + macro_bodies[j].radius) / 2.0; // Simple average radius
            double epsilon_sq = epsilon * epsilon;

            // 1. Calculate Distance Vector
            DynamoEngine::Vector2D distance = macro_bodies[i].position - macro_bodies[j].position;

            // 2. Calculate Distance Squared (r^2)
            double r_sq = distance.squareMagnitude();

            // 3. Calculate the Denominator Term (r^2 + epsilon^2)^(3/2)
            // The term inside the parenthesis: rSq + epsilonSq
            // The final term in the denominator: pow(rSq +
            // epsilonSq, 1.5)
            double denominator = sqrt(r_sq);

            total_energy -= GRAVITATIONAL_CONSTANT * macro_bodies[i].mass * macro_bodies[j].mass / denominator;
        }
    }
    std::cout << "Total E: " << total_energy << std::endl;
}

// --------- CLEANUP --------- //

void PhysicsSystem::cleanupParticles(GameState& game_state)
{
    auto& particles = game_state.getParticlesMutable();

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

void PhysicsSystem::cleanupMacroBodies(GameState& game_state)
{
    auto& particles = game_state.getMacroBodiesMutable();

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

void PhysicsSystem::updatePlayerPhysics(GameState& game_state, UIState& ui_state)
{
    game_state.getPlayerMutable().starship.applyRotation(ui_state);
    game_state.getPlayerMutable().starship.applyThrust(ui_state);
}

uint32_t PhysicsSystem::survivableFragmentCount(const GravitationalBody& body, uint32_t maxCount)
{
    double max_survivable = (OVERLAP_MARGIN * body.radius) / MIN_PARTICLE_RADIUS;
    return std::clamp(static_cast<uint32_t>(max_survivable * max_survivable), 1u, maxCount);
}

bool PhysicsSystem::hasRoomForParticles(const GameState& game_state, size_t new_particle_count) const
{
    size_t potential_particles = potentialParticleCount(game_state);
    return (potential_particles + new_particle_count) <= MAX_LIVE_PARTICLES;
}

size_t PhysicsSystem::potentialParticleCount(const GameState& game_state) const
{
    // A flat DEFAULT_FRAGMENT_COUNT per macro, not survivableFragmentCount: a macro that grows by accretion would
    // raise its exact count after it was budgeted, but it can never raise the flat one
    return liveParticleCount(game_state) + game_state.getMacroBodies().size() * DEFAULT_FRAGMENT_COUNT;
}