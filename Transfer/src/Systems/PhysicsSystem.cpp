// File: Transfer/src/Systems/PhysicsSystem.cpp

#include "PhysicsSystem.hpp"

PhysicsSystem::PhysicsSystem()
{
    // Initialize physics system variables if needed
}

PhysicsSystem::~PhysicsSystem() {}

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

static inline GravitationalBodyPair pickMassPair(GravitationalBody& a, GravitationalBody& b)
{
    if (abs(a.mass) == abs(b.mass))
        return {&a, &b, 1.0, true};
    if (abs(a.mass) > abs(b.mass))
        return {&a, &b, abs(a.mass / b.mass), false};
    return {&b, &a, abs(b.mass / a.mass), false};
}

static inline void populateGravBodyPropertiesFromInputState(GravitationalBody& gravBody, GameState& gameState,
                                                            InputState& inputState)
{
    // Default sets, position may be overwrriten if isCreatingWithInitialVelocity set to true
    gravBody.mass = inputState.selectedMass;
    gravBody.invMass = 1 / gravBody.mass;
    gravBody.radius = inputState.selectedRadius;
    gravBody.position = ScreenToWorldCoordinates(inputState.mouseCurrPosition, gameState.getCameraState());
    gravBody.previousPosition = gravBody.position;

    // Flags (will be moved to control from within the inputState)
    // Type flag
    gravBody.isMacro = true;

    // Property Flags
    gravBody.isAccretable = true;
    gravBody.isBounce = false;
    gravBody.isCollidable = true;
    gravBody.isForceStatic = false;
    gravBody.isFragment = false; // not a child of a collision
    gravBody.isMacroGhost = false;
    gravBody.isPreview = false;
    gravBody.isShatterable = true;
    gravBody.isTransient = false;

    // Visual Identifier Flags
    gravBody.isDust = false;
    gravBody.isGas = false;
    gravBody.isGravStar = false;
    gravBody.isMoon = false;
    gravBody.isPlanet = true; // Default for now until I implement procedural texture generation.

    if (inputState.isCreatingWithInitialVelocity)
    {
        gravBody.position = ScreenToWorldCoordinates(inputState.mouseDragStartPosition, gameState.getCameraState());
        gravBody.previousPosition = gravBody.position;
        gravBody.velocity =
            (inputState.mouseCurrPosition - inputState.mouseDragStartPosition) / gameState.getCameraState().zoom;
    }
}

// --------- CLEANUP METHOD --------- //

void PhysicsSystem::CleanUp()
{
    // Any necessary cleanup code for the physics system
}

void PhysicsSystem::UpdateGravBodyInstantiations(GameState& gameState, UIState& UIState)
{
    InputState& input_state = UIState.getMutableInputState();
    if (!input_state.UIInputConsumed)
    {
        if (input_state.isCreatingMacro)
        {
            createMacroBody(gameState,
                            input_state); // can inline replace with other create* methods for test. isCreatingMacro
                                          // should eventually be unique to just creating planets, etc.
            input_state.resetTransientFlags();
        }
    }
    else
    {
        // run limiters?
    }
    // Check if all Gravitational Bodies are supposed to be wiped
    if (input_state.clearAll)
    {
        gameState.getMacroBodiesMutable().clear();
        gameState.getParticlesMutable().clear();
        input_state.clearAll = false;
    }
}

void PhysicsSystem::updateAllForces(GameState& gameState)
{
    updateGravityForSystem(gameState);
    // Update other forces
}

void PhysicsSystem::updateGravityForSystem(GameState& gameState)
{
    std::vector<GravitationalBody>& macro_bodies = gameState.getMacroBodiesMutable();
    std::vector<GravitationalBody>& particles = gameState.getParticlesMutable();
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

    // for (size_t i = 0; i < num_particles; ++i)
    // {
    //     for (size_t j = i + 1; j < num_particles; ++j)
    //     {
    //         calculateGravity(particles[i], particles[j]);
    //     }
    // }
}

void PhysicsSystem::calculateGravity(GravitationalBody& first_body, GravitationalBody& second_body)
{
    if (firstWithinEpsilonOfSecond(first_body.mass, 0.0) || firstWithinEpsilonOfSecond(second_body.mass, 0.0))
    {
        return;
    }
    // if (firstWithinEpsilonOfSecond(first_body.radius, 0.0) || firstWithinEpsilonOfSecond(second_body.radius, 0.0))
    // {
    //     return;
    // }
    const double G = GRAVITATIONAL_CONSTANT;

    double softening_constant = (first_body.radius + second_body.radius) / 2.0;
    double epsilon_squared = softening_constant * softening_constant;

    // 1. Calculate direction Vector
    Vector2D direction_vector = second_body.position - first_body.position;

    // 2. Calculate Distance Squared
    double r_squared = direction_vector.square_magnitude();

    // 3. Calculate the Denominator Term (r^2 + epsilon^2)^(3/2)
    double denominator_1 = sqrt(r_squared + epsilon_squared);
    double denominator_3 = denominator_1 * denominator_1 * denominator_1;

    // 4. Calculate Coefficient (F = direction vector * G * m1 * m2 /
    // Denominator)
    double coefficient = (G * first_body.mass * second_body.mass) / denominator_3;

    // Force vector is C * direction_vector (r)
    Vector2D force = direction_vector * coefficient;

    // 6. Apply Forces (Newton's Third Law)
    if (!first_body.isForceStatic)
        first_body.netForce += force;

    if (!second_body.isForceStatic)
        second_body.netForce -= force;
}

void PhysicsSystem::UpdateSystemFrame(GameState& gameState, UIState& UIState)
{
    // Mental Model of System

    // Game updates gravitational body instantiations outside of Physics loop because the instantiations should be
    // rendered immediately (even with stopped time)

    // First, handle collisions, then integrate, then update forces, then integrate again, then finally cleanup

    handleCollisions(gameState);

    // promoteOversizedParticles(gameState);

    // Update forces (grav, later will add electromagnetic)
    integrateForwardsVelocityVerletPhase1(gameState);
    updateAllForces(gameState);
    integrateForwardsVelocityVerletPhase2(gameState);

    cleanupMacroBodies(gameState);
    cleanupParticles(gameState);

    // Verify calculation doesn't cook us
    // calculateTotalEnergy(gameState);
}
// --------- INTEGRATION METHODS --------- //
// Position integration from previous frame force.

// Honestly check if we need the velocity integration calculation at all. Are we not updated that at end of
// integration phase 2 of prev frame?

void PhysicsSystem::integrateForwardsVelocityVerletPhase1(GameState& gameState)
{
    std::vector<GravitationalBody>& particles = gameState.getParticlesMutable();
    std::vector<GravitationalBody>& macro_bodies = gameState.getMacroBodiesMutable();

    for (auto& particle : particles)
    {
        applyVelocityVerletPhase1(particle);
    }
    for (auto& macro_body : macro_bodies)
    {
        applyVelocityVerletPhase1(macro_body);
    }
}

void PhysicsSystem::applyVelocityVerletPhase1(GravitationalBody& gravBody)
{
    gravBody.previousPosition = gravBody.position;
    gravBody.prevForce = gravBody.netForce;
    bool hasMass = !(firstWithinEpsilonOfSecond(gravBody.mass, 0.0));
    if (hasMass && !gravBody.isForceStatic)
    {
        // Calculate the acceleration from the previous frame's final force
        Vector2D acceleration = gravBody.netForce * gravBody.invMass;
        gravBody.velocity += acceleration * (PHYSICS_TIME_STEP / 2); // Half of a full integrated step
    }
    // Step the position
    gravBody.position += gravBody.velocity * PHYSICS_TIME_STEP;
    // Reset to force 0 for next frame
    gravBody.netForce = Vector2D(0.0, 0.0);
}

void PhysicsSystem::integrateForwardsVelocityVerletPhase2(GameState& gameState)
{
    std::vector<GravitationalBody>& particles = gameState.getParticlesMutable();
    std::vector<GravitationalBody>& macro_bodies = gameState.getMacroBodiesMutable();

    for (auto& particle : particles)
    {
        applyVelocityVerletPhase2(particle);
    }
    for (auto& macro_body : macro_bodies)
    {
        applyVelocityVerletPhase2(macro_body);
    }
}

void PhysicsSystem::applyVelocityVerletPhase2(GravitationalBody& gravBody)
{
    bool hasMass = !(firstWithinEpsilonOfSecond(gravBody.mass, 0.0));

    if (!hasMass || gravBody.isForceStatic)
    {
        return;
    }
    else
    {
        Vector2D acceleration = gravBody.netForce * gravBody.invMass;
        gravBody.velocity += acceleration * (PHYSICS_TIME_STEP / 2.0); // Other half of full integrated step
    }
}

void PhysicsSystem::handleCollisions(GameState& gameState)
{
    auto& particles = gameState.getParticlesMutable();
    particles.reserve(particles.size() +
                      static_cast<size_t>(MAX_SIMULTANEOUS_SHATTERS_PER_TICK) * DEFAULT_FRAGMENT_COUNT * 2);

    handleMacroMacroCollisions(gameState);
    handleMacroParticleCollisions(gameState);
    handleParticleParticleCollisions(gameState);
}

void PhysicsSystem::handleMacroMacroCollisions(GameState& gameState)
{
    std::vector<GravitationalBody>& macroBodyList = gameState.getMacroBodiesMutable();
    size_t num_macro_bodies = macroBodyList.size();

    for (size_t i = 0; i < num_macro_bodies; ++i)
    {
        GravitationalBody& first_body = macroBodyList[i];
        if (!first_body.isCollidable || first_body.isMacroGhost || first_body.isMarkedForDeletion)
        {
            continue;
        }

        // j = i + 1: each unordered pair is visited exactly once (previously visited twice)
        for (size_t j = i + 1; j < num_macro_bodies; ++j)
        {
            GravitationalBody& second_body = macroBodyList[j];
            if (!second_body.isCollidable || second_body.isMacroGhost || second_body.isMarkedForDeletion)
            {
                continue;
            }

            CollisionInfo collision_info = getCollisionInfo(first_body, second_body);
            if (!collision_info.shouldCollide)
            {
                continue;
            }

            GravitationalBodyPair gravBodyPair = pickMassPair(first_body, second_body);
            handleDynamicCollision(gravBodyPair, collision_info, gameState);
        }
    }
}

void PhysicsSystem::handleMacroParticleCollisions(GameState& gameState)
{
    std::vector<GravitationalBody>& particles = gameState.getParticlesMutable();
    std::vector<GravitationalBody>& macroBodies = gameState.getMacroBodiesMutable();

    for (auto& particle : particles)
    {
        if (!particle.isCollidable || particle.isMarkedForDeletion)
        {
            continue;
        }

        for (auto& macroBody : macroBodies)
        {
            if (!macroBody.isCollidable || macroBody.isMacroGhost || macroBody.isMarkedForDeletion)
            {
                continue;
            }

            CollisionInfo collision_info = getCollisionInfo(particle, macroBody);
            if (!collision_info.shouldCollide)
            {
                continue;
            }

            GravitationalBodyPair gravBodyPair = pickMassPair(particle, macroBody);
            handleDynamicCollision(gravBodyPair, collision_info, gameState);
        }
    }
}

void PhysicsSystem::handleParticleParticleCollisions(GameState& gameState)
{
    std::vector<GravitationalBody>& particles = gameState.getParticlesMutable();

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

            GravitationalBodyPair gravBodyPair = pickMassPair(particles[i], particles[j]);
            handleDynamicCollision(gravBodyPair, collision_info, gameState);
        }
    }
}

void PhysicsSystem::handleDynamicCollision(GravitationalBodyPair& gravBodyPair, const CollisionInfo& collisionInfo,
                                           GameState& gameState)
{
    GravitationalBody& heavier = *gravBodyPair.heavierBody;
    GravitationalBody& lighter = *gravBodyPair.lighterBody;

    if (heavier.isBounce && lighter.isBounce)
    {
        handleElasticCollisions(lighter, heavier);
        return;
    }

    if (collisionInfo.shouldBlowUp)
    {
        if (!lighter.isShatterable)
        {
            handleElasticCollisions(lighter, heavier);
            return;
        }

        // Hard safety net: never let a shatter push live particle count past a cap, regardless
        // of any tuning elsewhere that might otherwise cascade (see point 2).
        if (gameState.getParticlesMutable().size() >= MAX_LIVE_PARTICLES)
        {
            handleElasticCollisions(lighter, heavier);
            return;
        }

        Vector2D towardLighter = (lighter.position - heavier.position).normalize();
        Vector2D impactPoint = heavier.position + towardLighter * heavier.radius;

        if (heavier.isShatterable && gravBodyPair.ratio <= MUTUAL_SHATTER_MASS_RATIO_THRESHOLD)
        {
            substituteWithParticlesFromImpact(heavier, gameState, DEFAULT_FRAGMENT_COUNT, impactPoint);
        }
        substituteWithParticlesFromImpact(lighter, gameState, DEFAULT_FRAGMENT_COUNT, impactPoint);
        return;
    }

    if (collisionInfo.absNormalSpeed >= MAX_ACCRETION_COLLISION_SPEED)
    {
        // Too fast to cleanly merge, not fast enough to shatter: always bounce.
        handleElasticCollisions(lighter, heavier);
        return;
    }

    // Gentle contact: merge if there's enough size disparity to look right, else bounce.
    double accretionRatioThreshold;
    if (heavier.isMacro && lighter.isMacro)
    {
        accretionRatioThreshold = MIN_BODY_BODY_ACCRETION_THRESHOLD_RATIO;
    }
    else if (heavier.isMacro || lighter.isMacro)
    {
        accretionRatioThreshold = MIN_BODY_PARTICLE_ACCRETION_THRESHOLD_RATIO;
    }
    else
    {
        accretionRatioThreshold = MIN_PARTICLE_PARTICLE_ACCRETION_THRESHOLD_RATIO;
    }

    if (lighter.isAccretable && gravBodyPair.ratio >= accretionRatioThreshold)
    {
        handleAccretion(gravBodyPair);
    }
    else
    {
        handleElasticCollisions(lighter, heavier);
    }
}

void PhysicsSystem::promoteOversizedParticles(GameState& gameState)
{
    auto& particles = gameState.getParticlesMutable();
    auto& macroBodies = gameState.getMacroBodiesMutable();

    for (auto& particle : particles)
    {
        if (!particle.isMarkedForPromotion)
        {
            continue;
        }

        gameState.incrementMaxIDInstantiated();

        GravitationalBody promoted;
        promoted.position = particle.position;
        promoted.previousPosition = particle.previousPosition;
        promoted.velocity = particle.velocity;
        promoted.netForce = particle.netForce;
        promoted.prevForce = particle.prevForce;
        promoted.mass = particle.mass;
        promoted.invMass = particle.invMass;
        promoted.radius = particle.radius;
        promoted.macroIdentifier = gameState.getMaxIDInstantiated();

        promoted.isMacro = true;
        promoted.isAccretable = true;
        promoted.isCollidable = true;
        promoted.isShatterable = true;
        promoted.isPlanet = true;

        macroBodies.push_back(promoted);

        // Defer actual removal to the existing end-of-tick particle cleanup rather than
        // erasing here, mid-iteration over the same particles vector.
        particle.isMarkedForDeletion = true;
    }
}
void PhysicsSystem::handleElasticCollisions(GravitationalBody& smallerBody, GravitationalBody& largerBody)
{
    if (smallerBody.isForceStatic && largerBody.isForceStatic)
    {
        return;
    }
    if (smallerBody.isForceStatic != largerBody.isForceStatic)
    {
        GravitationalBody& dyn = smallerBody.isForceStatic ? largerBody : smallerBody;
        GravitationalBody& stat = smallerBody.isForceStatic ? smallerBody : largerBody;

        Vector2D n = (dyn.position - stat.position).normalize();
        double v_n = dyn.velocity.dot(n);

        if (v_n < 0.0)
        {
            dyn.velocity -= n * (1.0 + ELASTIC_LOSS_FACTOR) * v_n;
        }
        return;
    }

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

    Vector2D normal_vector = r_vector / distance;
    double v_smaller_n = smallerBody.velocity.dot(normal_vector);
    double v_larger_n = largerBody.velocity.dot(normal_vector);
    double m_smaller = smallerBody.mass;
    double m_larger = largerBody.mass;

    double v_smaller_n_new =
        (v_smaller_n * (m_smaller - m_larger) + 2 * m_larger * v_larger_n) / (m_smaller + m_larger);
    double v_larger_n_new =
        (v_larger_n * (m_larger - m_smaller) + 2 * m_smaller * v_smaller_n) / (m_smaller + m_larger);

    smallerBody.velocity += normal_vector * (v_smaller_n_new - v_smaller_n) * ELASTIC_LOSS_FACTOR;
    largerBody.velocity += normal_vector * (v_larger_n_new - v_larger_n) * ELASTIC_LOSS_FACTOR;

    double penetration = (smallerBody.radius + largerBody.radius) - distance;
    if (penetration > 0.0)
    {
        constexpr double percent = 0.8;
        constexpr double slop = 0.01;
        double correction_magnitude = std::max(penetration - slop, 0.0) * percent;
        Vector2D correction = normal_vector * correction_magnitude;

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
            if (totalInvMass > EPSILON)
            {
                double smallerShare = smallerBody.invMass / totalInvMass;
                double largerShare = largerBody.invMass / totalInvMass;
                smallerBody.position -= correction * smallerShare;
                largerBody.position += correction * largerShare;
            }
            largerBody.previousPosition = largerBody.position;
            smallerBody.previousPosition = smallerBody.position;
        }
    }
}

void PhysicsSystem::handleAccretion(GravitationalBodyPair& gravBodyPair)
{
    GravitationalBody& heavier = *gravBodyPair.heavierBody;
    GravitationalBody& lighter = *gravBodyPair.lighterBody;

    if (heavier.isParticle)
    {
        return;
    }

    double newMass = heavier.mass + lighter.mass;
    heavier.radius *= pow(newMass / heavier.mass, 1.0 / 3.0);
    heavier.mass = newMass;
    heavier.invMass = 1.0 / heavier.mass;

    // Disabled because disabling promotion
    // if (!heavier.isMacro && heavier.radius >= PARTICLE_PROMOTION_RADIUS_THRESHOLD)
    // {
    //     heavier.isMarkedForPromotion = true;
    // }

    lighter.isMarkedForDeletion = true;
}

void PhysicsSystem::substituteWithParticles(GravitationalBody& originalBody, GameState& gameState,
                                            uint32_t targetFragmentCount)
{
    uint32_t num_particles = std::max<uint32_t>(1, targetFragmentCount);

    const double R = originalBody.radius;
    const Vector2D center = originalBody.position;
    const double originalMass = originalBody.mass;
    const Vector2D originalVelocity = originalBody.velocity;

    double densityFactor = (PI * R * R) / num_particles;
    const double fragmentRadius = OVERLAP_MARGIN * sqrt(densityFactor / PI);
    const double particleMass = originalMass / num_particles;

    auto& particles = gameState.getParticlesMutable();
    for (uint32_t k = 0; k < num_particles; ++k)
    {
        double r_k = R * sqrt((k + 0.5) / num_particles);
        double theta_k = k * GOLDEN_ANGLE;
        Vector2D pos_k = center + Vector2D{r_k * cos(theta_k), r_k * sin(theta_k)};

        GravitationalBody p;
        p.mass = particleMass;
        p.invMass = 1.0 / particleMass;
        p.radius = fragmentRadius;
        p.position = pos_k;
        p.previousPosition = pos_k;
        p.isFragment = true;
        p.isAccretable = true;
        p.isCollidable = true;
        p.macroIdentifier = originalBody.macroIdentifier;
        p.isParticle = true;
        p.velocity = originalVelocity * randomDouble(0.8, 1.1);

        particles.push_back(p);
    }

    originalBody.isMarkedForDeletion = true;
}

void PhysicsSystem::substituteWithParticlesFromImpact(GravitationalBody& originalBody, GameState& gameState,
                                                      uint32_t targetFragmentCount, const Vector2D& impactPoint)
{
    const double R = originalBody.radius;
    const double originalMass = originalBody.mass;

    auto& particles = gameState.getParticlesMutable();
    size_t start_index = particles.size();

    substituteWithParticles(originalBody, gameState, targetFragmentCount);

    // Grow radius with distance from the impact point: near-impact fragments stay small and
    // pulverized, far-side fragments stay large and coherent -- keyed off the actual contact
    // point (not the body's own center), so a graze reads differently from a direct hit.
    double totalWeightedVolume = 0.0;
    for (size_t i = start_index; i < particles.size(); ++i)
    {
        GravitationalBody& p = particles[i];
        double distance = (p.position - impactPoint).magnitude();
        double t = std::clamp(distance / (2.0 * R), 0.0, 1.0);
        double sizeMultiplier = 1.0 + t * (IMPACT_SKEW_GROWTH_FACTOR - 1.0) * randomDouble(0.5, 0.9);
        p.radius *= sizeMultiplier;
        totalWeightedVolume += p.radius * p.radius * p.radius;
    }

    // Re-normalize mass so the biased fragments still sum to originalMass, weighted by volume
    // (radius^3) to stay consistent with handleAccretion's mass<->radius law.
    for (size_t i = start_index; i < particles.size(); ++i)
    {
        GravitationalBody& p = particles[i];
        p.mass = originalMass * (p.radius * p.radius * p.radius) / totalWeightedVolume;
        p.invMass = 1.0 / p.mass;
    }
}
void PhysicsSystem::createMacroBody(GameState& gameState, InputState& inputState)
{
    std::vector<GravitationalBody>& macro_bodies = gameState.getMacroBodiesMutable();
    if (inputState.selectedRadius <= 1.0)
    {
        return;
    }
    if (firstWithinEpsilonOfSecond((inputState.selectedMass), 0.0))
    {
        return;
    }
    gameState.incrementMaxIDInstantiated();
    int newMacroBodyID = gameState.getMaxIDInstantiated();

    GravitationalBody macro_body;
    macro_body.macroIdentifier = newMacroBodyID;

    // Pass flags from inputState as possible.
    populateGravBodyPropertiesFromInputState(macro_body, gameState, inputState);

    // Now with populated flags, nudge particles out?
    macro_bodies.push_back(macro_body);
}

// --------- TOTAL ENERGY CALCULATION METHOD --------- //
void PhysicsSystem::calculateTotalEnergy(GameState& gameState)
{
    auto& macro_bodies = gameState.getMacroBodies();
    size_t num_macro_bodies = macro_bodies.size();

    double totalEnergy = 0.0;
    for (int i = 0; i < num_macro_bodies; ++i)
    {
        totalEnergy += macro_bodies[i].mass * macro_bodies[i].velocity.square_magnitude() / 2.0;
    }

    for (size_t i = 0; i < num_macro_bodies; i++)
    {
        for (size_t j = i + 1; j < num_macro_bodies; j++)
        {
            double epsilon = (macro_bodies[i].radius + macro_bodies[j].radius) / 2.0; // Simple average radius
            double epsilonSq = epsilon * epsilon;

            // 1. Calculate Distance Vector
            Vector2D distance = macro_bodies[i].position - macro_bodies[j].position;

            // 2. Calculate Distance Squared (r^2)
            double rSq = distance.square_magnitude();

            // 3. Calculate the Denominator Term (r^2 + epsilon^2)^(3/2)
            // The term inside the parenthesis: rSq + epsilonSq
            // The final term in the denominator: pow(rSq +
            // epsilonSq, 1.5)
            double denominator = sqrt(rSq);

            totalEnergy -= GRAVITATIONAL_CONSTANT * macro_bodies[i].mass * macro_bodies[j].mass / denominator;
        }
    }
    std::cout << "Total E: " << totalEnergy << std::endl;
}

// --------- CLEANUP GRAVITATIONAL BODIES METHODS --------- //

void PhysicsSystem::cleanupParticles(GameState& gameState)
{
    auto& particles = gameState.getParticlesMutable();

    // 1. Use std::remove_if to move all elements marked for deletion
    //    to the end of the vector. It returns an iterator to the new
    //    end.
    auto new_end =
        std::remove_if(particles.begin(), particles.end(), [](const GravitationalBody& p)
                       { return p.isMarkedForDeletion || firstWithinEpsilonOfSecond(p.mass, 0.0) || p.radius < 1; });

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
    auto new_end = std::remove_if(particles.begin(), particles.end(), [](const GravitationalBody& b)
                                  { return b.isMarkedForDeletion || firstWithinEpsilonOfSecond(b.mass, 0.0); });

    // 2. Use vector::erase to destroy the elements in the range
    // [new_end, particles.end())
    //    This efficiently shrinks the vector to the correct size.
    particles.erase(new_end, particles.end());
}
