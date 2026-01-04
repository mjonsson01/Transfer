// File: Transfer/src/Systems/PhysicsSystem.cpp

#include "PhysicsSystem.h"

PhysicsSystem::PhysicsSystem()
{
    // Initialize physics system variables if needed
}

PhysicsSystem::~PhysicsSystem()
{
    CleanUp();
}

// --------- ANONYMOUS NAMESPACE HELPER STRUCTS --------- //
struct CollisionInfo {
    double distance;
    Vector2D unitNormalVector;          // unit normal (from bodyA to bodyB)
    Vector2D relativeVelocityVector;    // vB - vA
    double normalSpeed;                 // signed speed along normal vector
    double absNormalSpeed;              // abs value of signed speed along normal vector
};

static inline CollisionInfo getCollisionInfo(const GravitationalBody& a,
                                             const GravitationalBody& b)
{
    Vector2D r_vector = b.position - a.position;
    double distance = r_vector.magnitude();
    Vector2D unit_normal_vector = (distance > 1e-8) ? (r_vector / distance) : Vector2D(1.0, 0.0);

    Vector2D relative_velocity_vector = b.velocity - a.velocity;
    double normal_speed = relative_velocity_vector.dot(unit_normal_vector);
    double abs_normal_speed = std::abs(normal_speed);
    return { distance, unit_normal_vector, relative_velocity_vector, normal_speed, abs_normal_speed };
}


struct GravitationalBodyPair {
    GravitationalBody* heavy;
    GravitationalBody* light;
    double ratio;   // heavy/light (>= 1)
    bool equal;
};

static inline GravitationalBodyPair pickMassPair(GravitationalBody& a, GravitationalBody& b)
{
    if (a.mass == b.mass) return { &a, &b, 1.0, true };
    if (a.mass > b.mass)  return { &a, &b, a.mass / b.mass, false };
    return { &b, &a, b.mass / a.mass, false };
}

// --------- SYSTEM-LEVEL METHOD --------- //

void PhysicsSystem::UpdateSystemFrame(GameState& gameState, UIState& UIState)
{

    // Check for dirty input that requires immediate updates to the Physics System
    updateGravBodyInstantiations(gameState, UIState);

    // Check if the Physics System is paused for early exit
    if (UIState.getInputState().isPhysicsPaused) return; 

    // ------------------------ MAIN PHYSICS LOOP ------------------------ //
    // Integrate first half step
    integrateForwardsPhase1(gameState);

    // Update all the gravity
    updateGravityForSystem(gameState);

    // Integrate second half step
    integrateForwardsPhase2(gameState);

    // Check total energy
    // calculateTotalEnergy(gameState);
    
    // Handle all the collisions
    handleCollisions(gameState);

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

void PhysicsSystem::updateGravBodyInstantiations(GameState& gameState, UIState& UIState)
{
    InputState& input_state = UIState.getMutableInputState();
    if (input_state.dirty)
    {
        if (input_state.isCreatingMacro)
            {
                createMacroBody(gameState, input_state);
                input_state.resetTransientFlags();
            }
        else if (input_state.isCreatingParticle)
            {
                createParticle(gameState, input_state);
                input_state.resetTransientFlags();
            }
        // else if (input_state.isCreatingParticleCluster)
        //     {
        //         createMacroBody(gameState, input_state);
        //         substituteWithParticles(gameState.getMacroBodiesMutable().back(), gameState);
        //         input_state.resetTransientFlags();
        //     }
    }
    // Check if all Gravitational Bodies are supposed to be wiped
    if (input_state.clearAll){
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
            calculateGravity(particles[i], macro_bodies[j]);
        }
    }
}

void PhysicsSystem::calculateGravity(GravitationalBody& body1, GravitationalBody& body2)
{   
    if (body1.mass == 0 || body2.mass == 0) return;
    static const double G = GRAVITATIONAL_CONSTANT;  

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

    // 4. Calculate Coefficient (F = direction vector * G * m1 * m2 / Denominator)
    double coefficient = (G * body1.mass * body2.mass) / denominator_3;

    // Force vector is C * direction_vector (r)
    Vector2D force = direction_vector * coefficient;

    // 6. Apply Forces (Newton's Third Law)
    body1.netForce += force;
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
        if (particle.isStatic) continue; // Skip any Static particles

        // Calculate current acceleration (a_t)
        Vector2D current_acceleration = particle.netForce / particle.mass;

        // Step 1. Kick 1: Update velocity by half the acceleration over one PHYSICS TIME STEP (dt)
        particle.velocity += current_acceleration * (PHYSICS_TIME_STEP / 2.0);

        // Step 2. Drift: Update position using the half-step velocity
        particle.previousPosition = particle.position;
        particle.position += particle.velocity * PHYSICS_TIME_STEP;
        
        // Reset netForce for the next step's force accumulation
        particle.netForce = Vector2D(0.0, 0.0); 
    }

    // Phase 1 integrate macro bodies
    for (auto& body : macro_bodies)
    {
        if (body.isStatic) continue;

        // Calculate current acceleration (a_t)
        Vector2D current_acceleration = body.netForce / body.mass;

        // Step 1. Kick 1: Update velocity by half the acceleration over one PHYSICS TIME STEP (dt)
        body.velocity += current_acceleration * (PHYSICS_TIME_STEP / 2.0);

        // Step 2. Drift: Update position using the half-step velocity
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
        if (particle.isStatic) continue;

        // Calculate new acceleration (a_t + dt) using the newly calculated netForce
        Vector2D new_acceleration = particle.netForce / particle.mass;

        // Step 3. Kick 2: Update velocity by half the acceleration over one PHYSICS TIME STEP (dt)
        particle.velocity += new_acceleration * (PHYSICS_TIME_STEP / 2.0);
    }
    
    // Phase 2 integrate macro bodies
    for (auto& body : macro_bodies)
    {
        if (body.isStatic) continue;

        // Calculate new acceleration (a_t + dt) using the newly calculated netForce
        Vector2D new_acceleration = body.netForce / body.mass;

        // Step 3. Kick 2: Update velocity by half the acceleration over one PHYSICS TIME STEP (dt)
        body.velocity += new_acceleration * (PHYSICS_TIME_STEP / 2.0);
    }
}


// --------- COLLISION METHODS --------- //

void PhysicsSystem::handleCollisions(GameState& gameState)
{
    // Get bodies and deal with body-body collisions
    auto& macro_bodies = gameState.getMacroBodiesMutable();
    size_t num_macro_bodies = macro_bodies.size();

    // Loop through all macro body collision pairs
    for (size_t i = 0; i < num_macro_bodies; ++i)
    {
        auto& body_a = macro_bodies[i];
        if (body_a.isMarkedForDeletion) continue; // ignore bodies that will be deleted this frame

        for (size_t j = i + 1; j < num_macro_bodies; ++j)
        {
            auto& body_b = macro_bodies[j];
            if (body_b.isMarkedForDeletion) continue; // ignore bodies that will be deleted this frame
            
            CollisionInfo info = getCollisionInfo(body_a, body_b);
           
            if (info.distance >= body_a.radius + body_b.radius) continue;

            GravitationalBodyPair pair = pickMassPair(body_a,body_b);

            if (info.absNormalSpeed < MIN_SHATTER_SPEED)
            {
                if (pair.ratio >= MIN_BODY_BODY_ACCRETION_THRESHOLD_RATIO)
                {
                    handleAccretion(*pair.light, *pair.heavy);
                }
                else
                {
                    handleElasticCollisions(*pair.light, *pair.heavy);
                }
            }
            else
            {
                if (pair.equal)
                {
                    handleDynamicExplosionCollision(body_a, body_b, gameState);
                }
                else
                {
                    substituteWithParticles(*pair.light, gameState);
                }
            }
        }
    }
    
    auto& particles = gameState.getParticlesMutable();
    for (auto& particle : particles)
    {
        if (particle.isMarkedForDeletion) continue;

        for (auto& body : macro_bodies)
        {
            if (body.isMarkedForDeletion) continue;

            auto info = getCollisionInfo(particle, body);
            if (info.distance >= particle.radius + body.radius) continue;

            if (info.absNormalSpeed > MAX_ACCRETION_COLLISION_SPEED)
            {
                // handleElasticCollisions(particle, body);
                handleAccretion(particle,body);
            }
            else
            {
                handleAccretion(particle, body);
                break; // Don't accrete same particle into multiple bodies this frame
            }
        }
    }
}


void PhysicsSystem::handleElasticCollisions(GravitationalBody& smallerBody, GravitationalBody& largerBody)
{
    Vector2D r_vector = largerBody.position - smallerBody.position;
    double distance = r_vector.magnitude();

    if (distance == 0) return;

    if (smallerBody.mass == 0 || largerBody.mass == 0) return;

    Vector2D normal_vector = r_vector / distance; // collision normal
    Vector2D v_smaller = smallerBody.velocity;
    Vector2D v_larger = largerBody.velocity;
    double m_smaller = smallerBody.mass;
    double m_larger = largerBody.mass;
    double v_smaller_n = v_smaller.dot(normal_vector);
    double v_larger_n = v_larger.dot(normal_vector);
    
    // 1D elastic collision formula (normal direction only)
    double v_smaller_n_new = (v_smaller_n*(m_smaller - m_larger) + 2*m_larger*v_larger_n) / (m_smaller + m_larger);
    double v_larger_n_new = (v_larger_n*(m_larger - m_smaller) + 2*m_smaller*v_smaller_n) / (m_smaller + m_larger);

    // Change in normal components
    Vector2D v_smaller_change = normal_vector * (v_smaller_n_new - v_smaller_n)*ELASTIC_LOSS_FACTOR;
    Vector2D v_larger_change = normal_vector * (v_larger_n_new - v_larger_n)*ELASTIC_LOSS_FACTOR;

    // Apply
    smallerBody.velocity = v_smaller + v_smaller_change;
    largerBody.velocity = v_larger + v_larger_change;
    
    double overlap = smallerBody.radius + largerBody.radius - distance;

    // Positional correction to disperse overlapping bodies
    if (abs(overlap) > 0){
        // Small slop to avoid micro jitter
        const double slop = 0.01;
        const double percent = 0.8; // correction strength

        double correction_mag = std::max(overlap - slop, 0.0) * percent;
        Vector2D correction = normal_vector * correction_mag;

        double totalMass = smallerBody.mass + largerBody.mass;

        smallerBody.position -= correction * (largerBody.mass / totalMass);
        largerBody.position += correction * (smallerBody.mass / totalMass);

        // Prevent fake velocity from Verlet-style position jumps
        smallerBody.previousPosition = smallerBody.position;
        largerBody.previousPosition = largerBody.position;
    }
}

void PhysicsSystem::handleDynamicExplosionCollision(GravitationalBody& body1, GravitationalBody& body2, GameState& gameState)
{
    const double R_1 = body1.radius;
    const double R_2 = body2.radius;
    const Vector2D center_1 = body1.position;
    const Vector2D center_2 = body2.position;
    const Vector2D center_of_explosion = (center_1 + center_2)/ 2.0;
    const double original_mass_1 = body1.mass;
    const double original_mass_2 = body2.mass;
    const Vector2D original_velocity_1 = body1.velocity;
    const Vector2D original_velocity_2 = body2.velocity;
    
    std::vector<Vector2D> pixel_positions_1; 

    // 1. RASTERIZATION (Find all pixel positions)
    // Assuming a 1-to-1 mapping where 1 unit = 1 pixel for simplicity.
    int R_int1 = static_cast<int>(R_1);
    for (int i = -R_int1; i <= R_int1; ++i) {
        for (int j = -R_int1; j <= R_int1; ++j) {
            // Check if (i, j) is inside the circle
            if ((i * i) + (j * j) <= (R_1 * R_1)) {
                // Store the particle's center position
                pixel_positions_1.push_back(center_1 + Vector2D{(double)i, (double)j});
            }
        }
    }
    std::vector<Vector2D> pixel_positions_2; 
    int R_int2 = static_cast<int>(R_2);
    for (int i = -R_int2; i <= R_int2; ++i) {
        for (int j = -R_int2; j <= R_int2; ++j) {
            // Check if (i, j) is inside the circle
            if ((i * i) + (j * j) <= (R_2 * R_2)) {
                // Store the particle's center position
                pixel_positions_2.push_back(center_2 + Vector2D{(double)i, (double)j});
            }
        }
    }
    
    // Check if we found any pixels for body2 (safety)
    if (pixel_positions_1.empty()){
        return;
    } 
        // Check if we found any pixels for body2 (safety)
    if (pixel_positions_2.empty()){
        return;
    } 

    // 2. MASS CALCULATION
    const double particleMass1 = original_mass_1 / pixel_positions_1.size();
    const double particleMass2 = original_mass_2 / pixel_positions_2.size();

    // 3. PARTICLE INSTANTIATION
    for (const auto& pos : pixel_positions_1) {
        GravitationalBody p;
        p.mass = particleMass1;
        p.radius = 1.0; 
        p.position = pos;
        p.isFragment = true;
        
        Vector2D r = pos - center_of_explosion;
        double dist = r.magnitude();

        Vector2D n = (dist > 1e-8) ? (r / dist) : Vector2D(1, 0);
        Vector2D t = Vector2D(-n.yVal, n.xVal);

        // Pick a characteristic radius for falloff.
        // This makes the falloff scale with the explosion size.
        double falloffRadius = 0.5 * (R_1 + R_2);          // or tweak: 0.25*(R_1+R_2), etc.
        double d = dist / (falloffRadius + 1e-8);

        // falloff in (0,1], where 1 at center, ~0 far away
        double falloff = 1.0 / (1.0 + d * d);            // nice smooth curve

        double factor = original_velocity_1.magnitude();
        double blast = factor * falloff;
        double shear = factor * falloff;

        double jitter = randomDouble(-1.0, 1.0) * falloff;

        // p.velocity = lerp(original_velocity_2, original_velocity_1, 0.5) + n * (blast + jitter) + t * (shear * randomDouble(-1.0, 1.0));
        Vector2D originalVel1 = original_velocity_1; // or original_velocity_2 depending on loop

        // Near center: use blended/mixed base. Far: keep original body velocity.
        Vector2D mixedBase1 = lerp(original_velocity_2, original_velocity_1, 0.5);
        Vector2D baseVel1   = lerp(originalVel1, mixedBase1, falloff);
        
        p.velocity = baseVel1 + n * (blast + jitter) + t * (shear * randomDouble(-1.0, 1.0));


        p.previousPosition = p.position - p.velocity * PHYSICS_TIME_STEP;

        gameState.getParticlesMutable().push_back(p);
    }
    for (const auto& pos : pixel_positions_2) {
        GravitationalBody p;
        p.mass = particleMass2;
        p.radius = 1.0; 
        p.position = pos;
        p.isFragment = true;
        
        Vector2D r = pos - center_of_explosion;
        double dist = r.magnitude();

        Vector2D n = (dist > 1e-8) ? (r / dist) : Vector2D(1, 0);
        Vector2D t = Vector2D(-n.yVal, n.xVal);

        // Pick a characteristic radius for falloff.
        // This makes the falloff scale with the explosion size.
        double falloffRadius = 0.5 * (R_1 + R_2);          // or tweak: 0.25*(R_1+R_2), etc.
        double d = dist / (falloffRadius + 1e-8);

        // falloff in (0,1], where 1 at center, ~0 far away
        double falloff = 1.0 / (1.0 + d * d);            // nice smooth curve
        double factor = original_velocity_2.magnitude();
        double blast = factor * falloff;
        double shear = factor * falloff;

        double jitter = randomDouble(-1.0, 1.0) * falloff;

        // p.velocity = lerp(original_velocity_1, original_velocity_2, 0.5) + n * (blast + jitter) + t * (shear * randomDouble(-1.0, 1.0));
        Vector2D originalVel2 = original_velocity_2; // or original_velocity_2 depending on loop

        // Near center: use blended/mixed base. Far: keep original body velocity.
        Vector2D mixedBase2 = lerp(original_velocity_2, original_velocity_1, 0.5);
        Vector2D baseVel2   = lerp(originalVel2, mixedBase2, falloff);
        p.velocity = baseVel2 + n * (blast + jitter) + t * (shear * randomDouble(-1.0, 1.0));

        p.previousPosition = p.position - p.velocity * PHYSICS_TIME_STEP;
        
        
        gameState.getParticlesMutable().push_back(p);
    }

    // Cleanup
    body1.isMarkedForDeletion = true; // Destroy the original macro body
    body2.isMarkedForDeletion = true; // Destroy the original macro body
}

void PhysicsSystem::handleAccretion(GravitationalBody& particle, GravitationalBody& body)
{
    // Accrete mass
    body.mass += particle.mass;
    // Add to radius
    body.radius = body.radius * pow((body.mass + particle.mass) / body.mass, 1.0/3.0);
    particle.isMarkedForDeletion = true;
}


void PhysicsSystem::createMacroBody(GameState& gameState, InputState& inputState)
{
    GravitationalBody body;
    body.mass = inputState.selectedMass;
    body.radius = inputState.selectedRadius;
    body.position = inputState.mouseCurrPosition;
    body.previousPosition = body.position;
    body.isPlanet = true;
    body.isStatic = inputState.isCreatingStatic;

    // --- NEW LOGIC: Nudge fragments out of the new body's radius ---
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

            // Normalize the displacement vector (safety check for dist=0, though unlikely here)
            Vector2D direction = (dist == 0) ? Vector2D(1.0, 0.0) : displacement / dist;

            // Apply the positional nudge
            particle.position += direction * required_separation;
            
            // Crucial for stability in Verlet integration:
            // Set prevPosition to the new position to prevent the next frame's 
            // velocity calculation from being massive and inaccurate.
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
    body.position = inputState.mouseCurrPosition;
    body.previousPosition = body.position;
    body.isDust = inputState.isCreatingDust;
    body.isStatic = inputState.isCreatingStatic;
    inputState.isCreatingDust = false;
    gameState.getParticlesMutable().push_back(body);
    inputState.dirty = false;
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
        totalEnergy += bodies[i].mass*bodies[i].velocity.square_magnitude()/2.0; 
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
            // The final term in the denominator: pow(rSq + epsilonSq, 1.5)
            double denominator = sqrt(rSq + epsilonSq); 

            totalEnergy -= GRAVITATIONAL_CONSTANT * bodies[i].mass*bodies[j].mass/denominator;
        }
    }
    std::cout<<"Total E: "<< totalEnergy << std::endl;
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
    for (int i = -R_int; i <= R_int; ++i) {
        for (int j = -R_int; j <= R_int; ++j) {
            // Check if (i, j) is inside the circle
            if ((i * i) + (j * j) <= (R * R)) {
                // Store the particle's center position
                pixelPositions.push_back(center + Vector2D{(double)i, (double)j});
            }
        }
    }
    
    // Check if we found any pixels (safety)
    if (pixelPositions.empty()){
        // std::cout<<"no pixels"<<std::endl;
        return;
    } 

    // 2. MASS CALCULATION
    const double particleMass = originalMass / pixelPositions.size();

    // 3. PARTICLE INSTANTIATION
    for (const auto& pos : pixelPositions) {
        GravitationalBody p;
        p.mass = particleMass;
        p.radius = 1.0; 
        p.position = pos;
        p.previousPosition = pos; // Initialize for Verlet integration
        p.isFragment = true;
        
        // Start with the original body's momentum
        p.velocity = originalVelocity;
        
        // Add a random outward "explosion" vector
        // This is necessary to make the cloud expand.
        Vector2D explosionVector = (pos - center); // Vector from center to particle
        // Normalize and scale by a collision factor (e.g., 0.1)
        explosionVector = explosionVector.normalize() * (0.0 * originalVelocity.magnitude()); 
        
        p.velocity += explosionVector;
        p.previousPosition = p.position - p.velocity * PHYSICS_TIME_STEP;


        gameState.getParticlesMutable().push_back(p);
    }

    // 4. CLEANUP
    originalBody.isMarkedForDeletion = true; // Destroy the original macroscopic body
}


// --------- CLEANUP GRAVITATIONAL BODIES METHODS --------- //

void PhysicsSystem::cleanupParticles(GameState& gameState)
{
    auto& particles = gameState.getParticlesMutable();

    // 1. Use std::remove_if to move all elements marked for deletion 
    //    to the end of the vector. It returns an iterator to the new end.
    auto new_end = std::remove_if(particles.begin(), particles.end(), 
        [](const GravitationalBody& p) {
            // The predicate returns true for elements to be 'removed' (moved to end)
            return p.isMarkedForDeletion; 
        }
    );

    // 2. Use vector::erase to destroy the elements in the range [new_end, particles.end())
    //    This efficiently shrinks the vector to the correct size.
    particles.erase(new_end, particles.end());
}

void PhysicsSystem::cleanupMacroBodies(GameState& gameState)
{
    auto& particles = gameState.getMacroBodiesMutable();

    // 1. Use std::remove_if to move all elements marked for deletion 
    //    to the end of the vector. It returns an iterator to the new end.
    auto new_end = std::remove_if(particles.begin(), particles.end(), 
        [](const GravitationalBody& b) {
            // The predicate returns true for elements to be 'removed' (moved to end)
            return b.isMarkedForDeletion; 
        }
    );

    // 2. Use vector::erase to destroy the elements in the range [new_end, particles.end())
    //    This efficiently shrinks the vector to the correct size.
    particles.erase(new_end, particles.end());
}
