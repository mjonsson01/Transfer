// File: Transfer/src/Systems/PhysicsSystem.cpp

#include "PhysicsSystem.hpp"

PhysicsSystem::PhysicsSystem()
{
    // Initialize physics system variables if needed
}

PhysicsSystem::~PhysicsSystem() {}

// --------- ANONYMOUS NAMESPACE HELPER STRUCTS --------- //
struct CollisionInfo
{
    double distance;
    Vector2D unitNormalVector;       // unit normal (from bodyA to bodyB)
    Vector2D relativeVelocityVector; // vB - vA
    double normalSpeed;              // signed speed along normal vector
    double absNormalSpeed;           // abs value of signed speed along normal vector
    bool shouldCollide;
    bool shouldBlowUp;
};

static inline CollisionInfo getCollisionInfo(const GravitationalBody& a, const GravitationalBody& b)
{
    Vector2D r_vector = b.position - a.position;
    double distance = r_vector.magnitude();
    Vector2D unit_normal_vector = (distance > 1e-8) ? (r_vector / distance) : Vector2D(1.0, 0.0);

    Vector2D relative_velocity_vector = b.velocity - a.velocity;
    double normal_speed = relative_velocity_vector.dot(unit_normal_vector);
    double abs_normal_speed = std::abs(normal_speed);
    bool shouldCollide = (distance < b.radius + a.radius);
    bool shouldBlowUp = (abs_normal_speed >= MIN_SHATTER_SPEED);
    return {distance,      unit_normal_vector, relative_velocity_vector, normal_speed, abs_normal_speed,
            shouldCollide, shouldBlowUp};
}

struct GravitationalBodyPair
{
    GravitationalBody* heavy;
    GravitationalBody* light;
    double ratio; // heavy/light (>= 1)
    bool equal;
};

static inline GravitationalBodyPair pickMassPair(GravitationalBody& a, GravitationalBody& b)
{
    if (a.mass == b.mass)
        return {&a, &b, 1.0, true};
    if (a.mass > b.mass)
        return {&a, &b, a.mass / b.mass, false};
    return {&b, &a, b.mass / a.mass, false};
}

// --------- SYSTEM-LEVEL METHOD --------- //

void PhysicsSystem::UpdateSystemFrame(GameState& gameState, UIState& UIState)
{
    // Check if the Physics System is paused for early exit
    // ------------------------ MAIN PHYSICS LOOP ------------------------ //

    // Handle all the collisions
    handleCollisions(gameState);

    // Advance Collision artifact lifetime
    auto& macro_bodies = gameState.getMacroBodiesMutable();
    for (auto& body : macro_bodies)
    {
        if (body.isTransient)
        {
            if (body.lifetime > 0)
            {
                --body.lifetime; // decrement lifetime
                body.mass *= 0.95;
                body.radius *= 0.999;
                body.velocity *= 0.999;
            }

            if (body.lifetime == 0)
            {
                body.isMarkedForDeletion = true; // flag for removal
                // body.mass = 0;
                // body.isCollidable = false;
            }
        }
    }
    // Integrate first half step
    integrateForwardsPhase1(gameState);

    // Update all the gravity
    updateGravityForSystem(gameState);

    // Integrate second half step
    integrateForwardsPhase2(gameState);

    // Check total energy
    // calculateTotalEnergy(gameState);

    // Clean up the marked for deletion particles and macro bodies
    cleanupParticles(gameState);
    cleanupMacroBodies(gameState);
}

// --------- CLEANUP METHOD --------- //

void PhysicsSystem::CleanUp()
{
    // Any necessary cleanup code for the physics system
}

// --------- BODY INSTANTIATION METHOD --------- //

void PhysicsSystem::UpdateGravBodyInstantiations(GameState& gameState, UIState& UIState)
{
    InputState& input_state = UIState.getMutableInputState();
    if (!input_state.UIInputConsumed)
    {
        if (input_state.isCreatingMacro)
        {
            createMacroBody(gameState, input_state);
            input_state.resetTransientFlags();
        }
    }
    else
    {
        // macroLimiter.reset();
        // particleLimiter.reset();
        // clusterLimiter.reset();
    }
    // Check if all Gravitational Bodies are supposed to be wiped
    if (input_state.clearAll)
    {
        gameState.getMacroBodiesMutable().clear();
        gameState.getParticlesMutable().clear();
        input_state.clearAll = false;
    }
}

// --------- GRAVITY METHODS --------- //

void PhysicsSystem::updateGravityForSystem(GameState& gameState)
{
    // Get Particles and Macro Bodies
    auto& particles = gameState.getParticlesMutable();
    int num_particles = particles.size();
    auto& macro_bodies = gameState.getMacroBodiesMutable();
    int num_macro_bodies = macro_bodies.size();

    // Update gravity for macro-macro interactions
    for (size_t i = 0; i < num_macro_bodies; i++)
    {
        for (size_t j = i + 1; j < num_macro_bodies; j++)
        {
            calculateGravity(macro_bodies[i], macro_bodies[j]);
        }
    }
    // Update gravity for particle-macro interactions
    for (size_t i = 0; i < num_particles; i++)
    {
        for (size_t j = 0; j < num_macro_bodies; j++)
        {
            // if (particles[i].macroIdentifier ==
            // macro_bodies[j].macroIdentifier)
            // {
            //     continue;
            // }
            calculateGravity(particles[i], macro_bodies[j]);
        }
    }
}

void PhysicsSystem::calculateGravity(GravitationalBody& body1, GravitationalBody& body2)
{
    if (firstWithinEpsilonOfSecond(body1.mass, 0.0) || firstWithinEpsilonOfSecond(body2.mass, 0.0))
    {
        return;
    }
    static const double G = GRAVITATIONAL_CONSTANT;

    if (body1.radius == 0 || body2.radius == 0) // this one might be okay because its cast to int? not sure
        return;
    double proxyScale = 1.0;
    // if (body1.macroIdentifier == body2.macroIdentifier)
    // {
    //     proxyScale = 0.10;
    // }
    // Define softening epsilon
    double epsilon = (body1.radius + body2.radius) / 2.0; // Simple average radius
    double epsilon_squared = epsilon * epsilon;

    // 1. Calculate direction Vector
    Vector2D direction_vector = body2.position - body1.position;

    // 2. Calculate Distance Squared
    double r_squared = direction_vector.square_magnitude();

    // 3. Calculate the Denominator Term (r^2 + epsilon^2)^(3/2)
    double denominator_1 = sqrt(r_squared + epsilon_squared);
    double denominator_3 = denominator_1 * denominator_1 * denominator_1;

    // 4. Calculate Coefficient (F = direction vector * G * m1 * m2 /
    // Denominator)
    double coefficient = proxyScale * (G * body1.mass * body2.mass) / denominator_3;

    // Force vector is C * direction_vector (r)
    Vector2D force = direction_vector * coefficient;

    // 6. Apply Forces (Newton's Third Law)
    if (!body1.isStatic)
        body1.netForce += force;

    if (!body2.isStatic)
        body2.netForce -= force;
}

// --------- INTEGRATION METHODS --------- //

void PhysicsSystem::integrateForwardsPhase1(GameState& gameState)
{
    // Get particles and macro bodies
    auto& particles = gameState.getParticlesMutable();
    auto& macro_bodies = gameState.getMacroBodiesMutable();

    // Phase 1 integrate particles
    for (auto& particle : particles)
    {
        bool hasMass = abs(particle.mass) > EPSILON;

        Vector2D acceleration(0.0, 0.0);

        if (!particle.isStatic && hasMass)
        {
            acceleration = particle.netForce / particle.mass;

            // Kick 1 (only if mass exists)
            particle.velocity += acceleration * (PHYSICS_TIME_STEP / 2.0);
        }

        // Drift (ALWAYS happens)
        particle.previousPosition = particle.position;
        particle.position += particle.velocity * PHYSICS_TIME_STEP;

        // Reset force
        particle.netForce = Vector2D(0.0, 0.0);
    }

    // Phase 1 integrate macro bodies
    for (auto& body : macro_bodies)
    {
        Vector2D current_acceleration(0.0, 0.0);
        if (body.isStatic)
        {
            // do nothing
        }
        else
        {
            if (abs(body.mass) >= EPSILON)
            {
                current_acceleration = body.netForce / body.mass;
            }
            body.velocity += current_acceleration * (PHYSICS_TIME_STEP / 2.0);
        }

        body.previousPosition = body.position;
        body.position += body.velocity * PHYSICS_TIME_STEP;

        // Reset netForce for the next step's force accumulation
        body.netForce = Vector2D(0.0, 0.0);
    }
}

void PhysicsSystem::integrateForwardsPhase2(GameState& gameState)
{
    // Get particles and macro bodies
    auto& particles = gameState.getParticlesMutable();
    auto& macro_bodies = gameState.getMacroBodiesMutable();

    // Phase 2 integrate particles
    for (auto& particle : particles)
    {
        if (particle.isStatic)
            continue;

        if (abs(particle.mass) > EPSILON)
        {
            Vector2D acceleration = particle.netForce / particle.mass;

            // Kick 2
            particle.velocity += acceleration * (PHYSICS_TIME_STEP / 2.0);
        }

        // else: zero-mass → no acceleration, velocity unchanged
    }
    // Phase 2 integrate macro bodies
    for (auto& body : macro_bodies)
    {
        Vector2D new_acceleration(0.0, 0.0);
        if (body.isStatic)
        {
            // no force integration, but velocity still carried through
            body.velocity += Vector2D(0.0, 0.0);
        }
        else
        {
            if (abs(body.mass) >= EPSILON)
            {
                new_acceleration = body.netForce / body.mass;
            }

            body.velocity += new_acceleration * (PHYSICS_TIME_STEP / 2.0);
        }
    }
}

void PhysicsSystem::handleCollisions(GameState& gameState)
{
    auto& macro_bodies = gameState.getMacroBodiesMutable();
    size_t num_macro_bodies = macro_bodies.size();

    // Macro-Macro collision pairs
    for (size_t i = 0; i < num_macro_bodies; ++i)
    {
        auto& body_a = macro_bodies[i];
        // skip if deleting
        if (body_a.isMarkedForDeletion)
        {
            continue;
        }
        // skip if not collidable
        if (!body_a.isCollidable)
        {
            continue;
        }

        // skip if they are ghost, since this is a macro macro loop
        if (body_a.isMacroGhost)
        {
            continue;
        }

        // compare with all other macro bodies
        for (size_t j = i + 1; j < num_macro_bodies; ++j)
        {
            auto& body_b = macro_bodies[j];
            // skip if deleting
            if (body_b.isMarkedForDeletion)
            {
                // std::cout << "Collision Branch Hit: 1" << std::endl;
                continue;
            }
            // skip if not collidable
            if (!body_b.isCollidable)
            {
                // std::cout << "Collision Branch Hit: 2" << std::endl;
                continue;
            }
            // skip if they are ghost, since this is a macro macro loop
            if (body_b.isMacroGhost)
            {
                // std::cout << "Collision Branch Hit: 3" << std::endl;
                continue;
            }

            CollisionInfo collision_info = getCollisionInfo(body_a, body_b);

            // if collision distance condition not satisfied skip
            if (!collision_info.shouldCollide)
            {
                // std::cout << "Collision Branch Hit: 4" << std::endl;
                continue;
            }
            // both are collidable and not macro ghosting
            GravitationalBodyPair pair = pickMassPair(body_a, body_b);

            // if both are bounce, elastic collide
            if (pair.light->isBounce && pair.heavy->isBounce)
            {
                // std::cout << "Collision Branch Hit: 5" << std::endl;
                handleElasticCollisions(*pair.light, *pair.heavy);
                continue;
            }

            if (collision_info.shouldBlowUp)
            {
                if (pair.heavy->isBounce) // even if going fast, should bounce
                // elastically, so shatter out the
                // particles of the light, which we
                // know is not a Bounce type
                {
                    if (pair.light->isShatterable)
                    {
                        //  std::cout << "Collision Branch Hit: 6" << std::endl;
                        substituteWithParticles(*pair.light, gameState);
                        continue;
                    }
                    else
                    {
                        //  std::cout << "Collision Branch Hit: 7" << std::endl;
                        handleElasticCollisions(*pair.light, *pair.heavy);
                        continue;
                    }
                }
                else if (pair.light->isBounce)
                {
                    if (pair.heavy->isShatterable)
                    {
                        //  std::cout << "Collision Branch Hit: 8" << std::endl;
                        substituteWithParticles(*pair.heavy, gameState);
                        continue;
                    }
                    else
                    {
                        // std::cout << "Collision Branch Hit: 9" << std::endl;
                        handleElasticCollisions(*pair.light, *pair.heavy);
                        continue;
                    }
                }
                else // neither are of bounce type
                {
                    if (pair.ratio <= 1.4)
                    {
                        // close in size, blow both of em up
                        // std::cout << "Collision Branch Hit: 10" << std::endl;
                        handleDynamicExplosionCollision(*pair.light, *pair.heavy, gameState);
                        continue;
                    }
                    else // just blow the smaller one into bits
                    {
                        //  std::cout << "Collision Branch Hit: 11" << std::endl;
                        substituteWithParticles(*pair.light, gameState);
                        continue;
                    }
                }
            }
            // Now we know not going fast enough to blow up, we can attempt to
            // accrete
            else
            {
                if (pair.light->isAccretable && pair.ratio >= MIN_BODY_BODY_ACCRETION_THRESHOLD_RATIO)
                {
                    // std::cout << "Collision Branch Hit: 12" << std::endl;
                    // handleAccretion(*pair.light, *pair.heavy);
                    substituteWithParticles(*pair.light, gameState);
                }
                else
                {
                    // std::cout << "Collision Branch Hit: 13" << std::endl;
                    handleElasticCollisions(*pair.light, *pair.heavy);
                }
            }
        }
    }

    auto& particles = gameState.getParticlesMutable();
    for (auto& particle : particles)
    {
        if (particle.isMarkedForDeletion)
        {
            // std::cout << "Collision Branch Hit: 14" << std::endl;
            continue;
        }
        if (!particle.isCollidable)
        {
            // std::cout << "Collision Branch Hit: 15" << std::endl;
            continue;
        }
        for (auto& macro_body : macro_bodies)
        {
            if (particle.macroIdentifier == macro_body.macroIdentifier)
            {
                continue;
            };
            if (macro_body.isMarkedForDeletion)
            {
                // std::cout << "Collision Branch Hit: 16" << std::endl;
                continue;
            }
            if (!macro_body.isCollidable)
            {
                // std::cout << "Collision Branch Hit: 17" << std::endl;
                continue;
            }
            CollisionInfo collision_info = getCollisionInfo(macro_body, particle);
            if (collision_info.shouldCollide)
            {
                if (macro_body.isBounce)
                {
                    // std::cout << "Collision Branch Hit: 18" << std::endl;
                    handleElasticCollisions(particle, macro_body);
                    continue;
                }
                else
                {
                    if (particle.isAccretable)
                    {
                        // std::cout << "Collision Branch Hit: 19" << std::endl;
                        if (particle.mass > macro_body.mass)
                        {
                            handleElasticCollisions(particle, macro_body);
                        }
                        else
                        {
                            handleAccretion(particle, macro_body);
                        }
                        break; // don't accrete the same particle twice
                    }
                }
            }
            else
            {
                // std::cout << "Collision Branch Hit: 20" << std::endl;
                continue;
            }
        }
    }
}

// newest elastic collision

void PhysicsSystem::handleElasticCollisions(GravitationalBody& smallerBody, GravitationalBody& largerBody)
{
    if (smallerBody.macroIdentifier == largerBody.macroIdentifier)
    {
        return;
    }
    if (smallerBody.isStatic && largerBody.isStatic)
    {
        return;
    }
    // if (smallerBody.isStatic xor largerBody.isStatic)
    // {
    //     GravitationalBody& dyn = smallerBody.isStatic ? largerBody : smallerBody;
    //     GravitationalBody& stat = smallerBody.isStatic ? smallerBody : largerBody;

    //     Vector2D n = (dyn.position - stat.position).normalize();
    //     double v_n = dyn.velocity.dot(n);

    //     if (v_n < 0.0)
    //     {
    //         dyn.velocity -= n * (1.0 + ELASTIC_LOSS_FACTOR) * v_n;
    //     }
    //     return;
    // }

    if (firstWithinEpsilonOfSecond(smallerBody.mass, 0.0) || firstWithinEpsilonOfSecond(largerBody.mass, 0.0))
    {
        return;
    }

    Vector2D r_vector = largerBody.position - smallerBody.position;
    double distance = r_vector.magnitude();

    if (firstWithinEpsilonOfSecond(distance, 0.0))
    {
        return;
    }

    Vector2D normal_vector = r_vector / distance; // collision normal
    Vector2D v_smaller = smallerBody.velocity;
    Vector2D v_larger = largerBody.velocity;
    double m_smaller = smallerBody.mass;
    double m_larger = largerBody.mass;
    double v_smaller_n = v_smaller.dot(normal_vector);
    double v_larger_n = v_larger.dot(normal_vector);

    // 1D elastic collision formula (normal direction only)
    double v_smaller_n_new =
        (v_smaller_n * (m_smaller - m_larger) + 2 * m_larger * v_larger_n) / (m_smaller + m_larger);
    double v_larger_n_new =
        (v_larger_n * (m_larger - m_smaller) + 2 * m_smaller * v_smaller_n) / (m_smaller + m_larger);

    Vector2D v_smaller_change = normal_vector * (v_smaller_n_new - v_smaller_n) * ELASTIC_LOSS_FACTOR;
    Vector2D v_larger_change = normal_vector * (v_larger_n_new - v_larger_n) * ELASTIC_LOSS_FACTOR;

    // Apply
    smallerBody.velocity = v_smaller + v_smaller_change;
    largerBody.velocity = v_larger + v_larger_change;

    double penetration = (smallerBody.radius + largerBody.radius) - distance;

    if (penetration > 0.0)
    {
        constexpr double percent = 0.8; // correction strength
        constexpr double slop = 0.01;   // small tolerance

        double correction_magnitude = std::max(penetration - slop, 0.0) * percent;

        Vector2D correction = normal_vector * correction_magnitude;

        // static handling
        if (smallerBody.isStatic && !largerBody.isStatic)
        {
            largerBody.position += correction;
            largerBody.previousPosition = largerBody.position;
        }
        else if (!smallerBody.isStatic && largerBody.isStatic)
        {
            smallerBody.position -= correction;
            smallerBody.previousPosition = smallerBody.position;
        }
        else if (!smallerBody.isStatic && !largerBody.isStatic)
        {
            smallerBody.position -= correction * 0.5;
            largerBody.position += correction * 0.5;
            largerBody.previousPosition = largerBody.position;
            smallerBody.previousPosition = smallerBody.position;
        }
    }
}

void PhysicsSystem::handleDynamicExplosionCollision(GravitationalBody& macroBody1, GravitationalBody& macroBody2,
                                                    GameState& gameState)
{
    // Buffers to prevent iterator invalidation
    std::vector<GravitationalBody> pending_macros;

    // instantiate collision proxy for both bodies. Consider a shrinking
    // proxy?Add later?
    double distance = (macroBody1.position - macroBody2.position).magnitude();
    // int proxyBody1LifeTime = TARGET_FPS * distance / macroBody1.velocity.magnitude();

    // int proxyBody2LifeTime = TARGET_FPS * distance / macroBody2.velocity.magnitude();
    int proxyBody1LifeTime = 200;
    int proxyBody2LifeTime = 200;

    // int proxyBody2LifeTime = proxyBody1LifeTime;

    // shattered bodies will contain child particles with
    // the same id. Then a skipout of the logic if id
    // matches macro (in elastic collsion)
    GravitationalBody proxyForMB1;
    gameState.incrementMaxIDInstantiated();
    populateCollisionProxyFromMacroBody(macroBody1, proxyForMB1);
    proxyForMB1.lifetime = proxyBody1LifeTime;
    // proxyForMB1.lifetime = 1000;
    GravitationalBody proxyForMB2;
    gameState.incrementMaxIDInstantiated();
    populateCollisionProxyFromMacroBody(macroBody2, proxyForMB2);
    proxyForMB2.lifetime = proxyBody2LifeTime;
    // proxyForMB2.lifetime = 1000;
    pending_macros.push_back(proxyForMB1);
    pending_macros.push_back(proxyForMB2);

    substituteWithParticles(macroBody1, gameState);
    substituteWithParticles(macroBody2, gameState);

    auto& macros = gameState.getMacroBodiesMutable();
    macros.insert(macros.end(), pending_macros.begin(), pending_macros.end());
}

void PhysicsSystem::handleAccretion(GravitationalBody& particle, GravitationalBody& body)
{
    // Accrete mass
    body.mass += particle.mass;
    // Add to radius
    body.radius = body.radius * pow((body.mass + particle.mass) / body.mass, 1.0 / 3.0);
    particle.isMarkedForDeletion = true;
}

void PhysicsSystem::createMacroBody(GameState& gameState, InputState& inputState)
{
    if (inputState.selectedRadius <= 1.0)
    {
        // need to throw error toast or something somehow
        std::cout << "early exit" << std::endl;
        return;
    }
    gameState.incrementMaxIDInstantiated();
    int newMacroBodyID = gameState.getMaxIDInstantiated();

    GravitationalBody body;
    body.mass = inputState.selectedMass;
    body.radius = inputState.selectedRadius;
    body.position = ScreenToWorldCoordinates(inputState.mouseCurrPosition, gameState.getCameraState());
    body.previousPosition = body.position;
    // TODO PASS FLAGS HERE
    body.isPlanet = true;
    body.isMacro = true;
    body.isShatterable = true;
    body.isCollidable = inputState.isCreatingCollidable;
    body.isStatic = inputState.isCreatingStatic;
    body.isMacroGhost = inputState.isCreatingMacroGhost;
    body.macroIdentifier = newMacroBodyID;
    body.isAccretable = true;

    if (inputState.isCreatingWithInitialVelocity)
    {
        body.position = ScreenToWorldCoordinates(inputState.mouseDragStartPosition, gameState.getCameraState());
        body.previousPosition = body.position;
        body.velocity =
            (inputState.mouseCurrPosition - inputState.mouseDragStartPosition) / gameState.getCameraState().zoom;
    }
    // Nudge fragments out of the new body's radius
    auto& particles = gameState.getParticlesMutable();
    double nudge_factor = 1.01; // Nudge fragments out by 1% more than the radius

    for (auto& particle : particles)
    {
        Vector2D displacement = particle.position - body.position;
        double dist = displacement.magnitude();
        double min_dist = body.radius + particle.radius;

        if (dist < min_dist)
        {
            // Calculate the required outward displacement
            // We want the new distance to be min_dist * nudge_factor
            double required_separation = (min_dist * nudge_factor) - dist;

            // Normalize the displacement vector (safety check for
            // dist=0, though unlikely here)
            Vector2D direction = (dist == 0) ? Vector2D(1.0, 0.0) : displacement / dist;

            // Apply the positional nudge
            particle.position += direction * required_separation;

            // Crucial for stability in Verlet integration:
            // Set prevPosition to the new position to prevent the
            // next
            // frame's velocity calculation from being massive and
            // inaccurate.
            particle.previousPosition = particle.position;
        }
    }

    gameState.getMacroBodiesMutable().push_back(body);
}

void PhysicsSystem::createParticle(GameState& gameState, InputState& inputState)
{
    GravitationalBody body;
    body.mass = inputState.selectedMass;
    body.radius = inputState.selectedRadius;
    body.position = ScreenToWorldCoordinates(inputState.mouseCurrPosition, gameState.getCameraState());
    body.previousPosition = body.position;
    body.isDust = inputState.isCreatingDust;
    body.isStatic = inputState.isCreatingStatic;
    body.isCollidable = true;
    inputState.isCreatingDust = false;
    gameState.getParticlesMutable().push_back(body);
}

void PhysicsSystem::createParticleCluster(GameState& gameState, InputState& inputState)
{
    GravitationalBody body;
    body.mass = inputState.selectedMass;
    body.radius = inputState.selectedRadius;
    body.position = ScreenToWorldCoordinates(inputState.mouseCurrPosition, gameState.getCameraState());
    body.previousPosition = body.position;
    body.isPlanet = true;
    body.isStatic = inputState.isCreatingStatic;

    substituteWithParticles(body, gameState);
}
// --------- TOTAL ENERGY CALCULATION METHOD --------- //

void PhysicsSystem::calculateTotalEnergy(GameState& gameState)
{
    // Calculate total energy logic here
    auto& bodies = gameState.getMacroBodies();
    int num_bodies = bodies.size();
    double totalEnergy = 0.0;
    for (int i = 0; i < num_bodies; ++i)
    {
        totalEnergy += bodies[i].mass * bodies[i].velocity.square_magnitude() / 2.0;
    }

    for (size_t i = 0; i < num_bodies; i++)
    {
        for (size_t j = i + 1; j < num_bodies; j++)
        {
            double epsilon = (bodies[i].radius + bodies[j].radius) / 2.0; // Simple average radius
            double epsilonSq = epsilon * epsilon;

            // 1. Calculate Distance Vector
            Vector2D distance = bodies[i].position - bodies[j].position;

            // 2. Calculate Distance Squared (r^2)
            double rSq = distance.square_magnitude();

            // 3. Calculate the Denominator Term (r^2 + epsilon^2)^(3/2)
            // The term inside the parenthesis: rSq + epsilonSq
            // The final term in the denominator: pow(rSq +
            // epsilonSq, 1.5)
            double denominator = sqrt(rSq + epsilonSq);

            totalEnergy -= GRAVITATIONAL_CONSTANT * bodies[i].mass * bodies[j].mass / denominator;
        }
    }
    std::cout << "Total E: " << totalEnergy << std::endl;
}

// --------- PARTICLE SUBSTITUTION METHOD --------- //

void PhysicsSystem::substituteWithParticles(GravitationalBody& originalBody, GameState& gameState)
{
    // A simplified value for the collision radius (R)
    const double R = originalBody.radius;
    const Vector2D center = originalBody.position;
    const double originalMass = originalBody.mass;
    const Vector2D originalVelocity = originalBody.velocity;

    std::vector<Vector2D> pixelPositions;

    // 1. RASTERIZATION (Find all pixel positions)
    // Assuming a 1-to-1 mapping where 1 unit = 1 pixel for simplicity.
    int R_int = static_cast<int>(R);
    for (int i = -R_int; i <= R_int; ++i)
    {
        for (int j = -R_int; j <= R_int; ++j)
        {
            // Check if (i, j) is inside the circle
            if ((i * i) + (j * j) <= (R * R))
            {
                // Store the particle's center position
                pixelPositions.push_back(center + Vector2D{(double)i, (double)j});
            }
        }
    }

    // Check if we found any pixels (safety)
    if (pixelPositions.empty())
    {
        return;
    }

    // 2. MASS CALCULATION
    const double particleMass = originalMass / pixelPositions.size();

    // 3. PARTICLE INSTANTIATION
    for (const auto& pos : pixelPositions)
    {
        GravitationalBody p;
        p.mass = particleMass;
        p.radius = 1.0;
        p.position = pos;
        p.previousPosition = pos; // No offset here
        p.isFragment = true;
        p.isAccretable = true;
        p.isCollidable = true;
        p.macroIdentifier = originalBody.macroIdentifier;

        // inherit velocity exactly
        p.velocity = originalVelocity;
        // p.velocity = originalVelocity * randomDouble(0.95, 1.05);
        // p.velocity = originalVelocity * randomDouble(0.5, 0.85);

        gameState.getParticlesMutable().push_back(p);
    }

    // 4. CLEANUP
    originalBody.isMarkedForDeletion = true; // Destroy the original macro body
}

void PhysicsSystem::populateCollisionProxyFromMacroBody(GravitationalBody& originalMacroBody,
                                                        GravitationalBody& proxyBody)
{
    proxyBody.position = originalMacroBody.position;
    proxyBody.previousPosition = proxyBody.position;
    proxyBody.isMacroGhost = true;
    proxyBody.isMacro = true;
    // proxyBody.lifetime = TARGET_FPS; // ~1 second for now, still need
    // to fix.

    // proxy body lifetime should be the amount of time it is predicted
    // to take to reach the collision center. needs to be instantiated
    // at the collision level.
    proxyBody.mass = 1.0 * originalMacroBody.mass;
    proxyBody.radius = originalMacroBody.radius * 1.0;
    proxyBody.velocity = originalMacroBody.velocity * 1.0;
    proxyBody.isStatic = true; // should follow a single trajectory
    proxyBody.isBounce = true;
    proxyBody.isShatterable = false;
    proxyBody.isCollidable = true;
    proxyBody.isTransient = true;
    proxyBody.netForce = {0.0, 0.0};
    proxyBody.prevForce = {0.0, 0.0};
    proxyBody.visible = false; 
    proxyBody.macroIdentifier = originalMacroBody.macroIdentifier;
}
// --------- CLEANUP GRAVITATIONAL BODIES METHODS --------- //

void PhysicsSystem::cleanupParticles(GameState& gameState)
{
    auto& particles = gameState.getParticlesMutable();

    // 1. Use std::remove_if to move all elements marked for deletion
    //    to the end of the vector. It returns an iterator to the new
    //    end.
    auto new_end = std::remove_if(particles.begin(), particles.end(),
                                  [](const GravitationalBody& p)
                                  {
                                      // The predicate returns true for
                                      // elements to be 'removed' (moved
                                      // to end)
                                      return p.isMarkedForDeletion;
                                  });

    // 2. Use vector::erase to destroy the elements in the range
    // [new_end, particles.end())
    //    This efficiently shrinks the vector to the correct size.
    particles.erase(new_end, particles.end());
}

void PhysicsSystem::cleanupMacroBodies(GameState& gameState)
{
    auto& particles = gameState.getMacroBodiesMutable();

    // 1. Use std::remove_if to move all elements marked for deletion
    //    to the end of the vector. It returns an iterator to the new
    //    end.
    auto new_end = std::remove_if(particles.begin(), particles.end(),
                                  [](const GravitationalBody& b)
                                  {
                                      // The predicate returns true for
                                      // elements to be 'removed' (moved
                                      // to end)
                                      return b.isMarkedForDeletion;
                                  });

    // 2. Use vector::erase to destroy the elements in the range
    // [new_end, particles.end())
    //    This efficiently shrinks the vector to the correct size.
    particles.erase(new_end, particles.end());
}
