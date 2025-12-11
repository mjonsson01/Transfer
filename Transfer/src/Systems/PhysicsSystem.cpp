// File: Transfer/src/Systems/PhysicsSystem.cpp

#include "PhysicsSystem.h"

PhysicsSystem::PhysicsSystem()
{
    // Initialize physics system variables if needed
}

PhysicsSystem::~PhysicsSystem()
{
    // Cleanup if necessary
}
void PhysicsSystem::CleanUp()
{
    // Any necessary cleanup code for the physics system
}

void PhysicsSystem::UpdateSystemFrame(GameState& state, UIState& UIState)
{
    InputState& inputState = UIState.getMutableInputState();
    if (inputState.dirty)
    {
        if (inputState.isCreatingCluster)
            {
                inputState.dirty = false;
                inputState.isCreatingCluster = false;
            }
        if (inputState.isCreatingPlanet)
            {
                createPlanet(state, inputState);
                inputState.dirty = false;
                inputState.isCreatingPlanet = false;
            }
        if (inputState.isCreatingDust)
            {
                // Create dust logic here
                inputState.dirty = false;
                inputState.isCreatingDust = false;
            }
    }
    integrateForwards_Phase1(state);
    
    auto& particles = state.getParticlesMutable();
    int num_particles = particles.size();

    auto& bodies = state.getMacroBodiesMutable();
    int num_bodies = bodies.size();

    // std::cout<<"num bodies: "<<num_bodies<<", num particles: " << num_particles<<std::endl;

    for (size_t i = 0; i < num_bodies; i++)
    {
        for (size_t j = i + 1; j < num_bodies; j++)
        {
            calculateGravity(bodies[i], bodies[j]);
            // std::cout<<"i: "<<i<<"j: "<<j<<std::endl;
        }
    }

    for (size_t i = 0; i < num_particles; i++)
        for (size_t j = 0; j < num_bodies; j++)
        {
            calculateGravityForSmallFragments(particles[i], bodies[j]);
        }


    // integrateForwards(state);
    integrateForwards_Phase2(state);

    // check total energy
    // calculateTotalEnergy(state);

    handleCollisions(state);
    cleanupParticles(state);
    cleanupMacroBodies(state);

}


void PhysicsSystem::handleElasticCollisions(GravitationalBody& particle, GravitationalBody& body)
{
    Vector2D d = body.position - particle.position;
    double dist = d.magnitude();

    if (dist == 0) return;
    if (particle.mass == 0 || body.mass == 0) return;
    Vector2D n = d / dist; // collision normal
    Vector2D vP = particle.velocity;
    Vector2D vB = body.velocity;
    double mP = particle.mass;
    double mB = body.mass;
    double vP_n = vP.dot(n);
    double vB_n = vB.dot(n);
    // 1D elastic collision formula (normal direction only)
    double vP_n_new = (vP_n*(mP - mB) + 2*mB*vB_n) / (mP + mB);
    double vB_n_new = (vB_n*(mB - mP) + 2*mP*vP_n) / (mP + mB);

//     // Change in normal components
    Vector2D vP_change = n * (vP_n_new - vP_n)*0.9;
    Vector2D vB_change = n * (vB_n_new - vB_n)*0.9;

//     // Apply
    particle.velocity = vP + vP_change;
    body.velocity = vB + vB_change;

    // // *** NEW VELOCITY-BASED POSITIONAL CORRECTION ***
    // double overlap = particle.radius + body.radius - dist;
    // if (overlap > 0) {
    //     // Calculate the relative velocity needed to cancel the overlap *over the next time step*
    //     Vector2D separation_velocity = n * (overlap / PHYSICS_TIME_STEP);

    //     // Apply this additional velocity to push them apart
    //     particle.velocity -= separation_velocity * 0.5; // push particle away
    //     body.velocity += separation_velocity * 0.5; // push body away
        
    //     // You might need to make PHYSICS_TIME_STEP accessible or pass it in.
    //     // Also, you may need a small epsilon check on overlap to avoid division by zero if overlap is tiny.
    // }
    // *** VELOCITY-BASED POSITIONAL CORRECTION ***
    double overlap = particle.radius + body.radius - dist;
    
    if (overlap > 0.0) {
        // Calculate the required velocity (V_sep) to correct the overlap (P) over 1 timestep (dt)
        // V_sep = P / dt. The factor of 1.05 gives a slight margin of separation.
        double sep_velocity_mag = (overlap / PHYSICS_TIME_STEP) * 1.05; 
        
        Vector2D separation_velocity = n * sep_velocity_mag; 

        // Apply the separation velocity, weighted by mass
        double mSum = particle.mass + body.mass;
        
        // Push particle away from the body
        particle.velocity -= separation_velocity * (body.mass / mSum); 
        // Push body away from the particle
        body.velocity += separation_velocity * (particle.mass / mSum); 
    }
}

void PhysicsSystem::calculateTotalEnergy(GameState& state)
{
    // Calculate total energy logic here
    auto& bodies = state.getMacroBodies();
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
void PhysicsSystem::cleanupParticles(GameState& state)
{
    auto& particles = state.getParticlesMutable();

    // 1. Use std::remove_if to move all elements marked for deletion 
    //    to the end of the vector. It returns an iterator to the new end.
    auto new_end = std::remove_if(particles.begin(), particles.end(), 
        [](const GravitationalBody& p) {
            // The predicate returns true for elements to be 'removed' (moved to end)
            return p.markedForDeletion; 
        }
    );

    // 2. Use vector::erase to destroy the elements in the range [new_end, particles.end())
    //    This efficiently shrinks the vector to the correct size.
    particles.erase(new_end, particles.end());
}

void PhysicsSystem::cleanupMacroBodies(GameState& state)
{
    auto& particles = state.getMacroBodiesMutable();

    // 1. Use std::remove_if to move all elements marked for deletion 
    //    to the end of the vector. It returns an iterator to the new end.
    auto new_end = std::remove_if(particles.begin(), particles.end(), 
        [](const GravitationalBody& b) {
            // The predicate returns true for elements to be 'removed' (moved to end)
            return b.markedForDeletion; 
        }
    );

    // 2. Use vector::erase to destroy the elements in the range [new_end, particles.end())
    //    This efficiently shrinks the vector to the correct size.
    particles.erase(new_end, particles.end());
}
void PhysicsSystem::handleCollisions(GameState& state)
{
    auto& bodies = state.getMacroBodiesMutable();
    for (size_t i = 0; i < bodies.size(); ++i)
    {
        for (size_t j = i + 1; j < bodies.size(); ++j)
        {
                double distance = (bodies[i].position - bodies[j].position).magnitude();
                if (distance < (bodies[i].radius + bodies[j].radius))
                {
                    // Rework logic later?
                    substituteWithParticles(bodies[j], state);
                    substituteWithParticles(bodies[i], state);
                }
        }
    }
    auto& particles = state.getParticlesMutable();
    int num_particles = particles.size();

    for (auto& particle : particles)
    {
        for (auto& body : bodies){
            double distance = (particle.position - body.position).magnitude();
            
            if (distance < (particle.radius + body.radius))
            {
                // Call the function containing the fixed logic
                handleElasticCollisions(particle, body); 
            }
            // std::cout<<"Elastic Collision"<<std::endl;
        }
    }
}



// 1. Kick (Half-Step Velocity Update) & Drift (Position Update)
void PhysicsSystem::integrateForwards_Phase1(GameState& state)
{
    auto& particles = state.getParticlesMutable();
    auto& bodies = state.getMacroBodiesMutable();
    for (auto& particle : particles)
    {
        if (particle.isStatic) continue;

        // Calculate current acceleration (a_t)
        Vector2D current_acceleration = particle.netForce / particle.mass;

        // 1. Kick 1: Update velocity by half the acceleration
        particle.velocity += current_acceleration * (PHYSICS_TIME_STEP / 2.0);

        // 2. Drift: Update position using the half-step velocity
        particle.prevPosition = particle.position;
        particle.position += particle.velocity * PHYSICS_TIME_STEP;

        // Reset netForce for the next step's force accumulation
        particle.netForce = Vector2D(0.0, 0.0); 
    }

    for (auto& body : bodies)
    {
        if (body.isStatic) continue;

        // Calculate current acceleration (a_t)
        Vector2D current_acceleration = body.netForce / body.mass;

        // 1. Kick 1: Update velocity by half the acceleration
        body.velocity += current_acceleration * (PHYSICS_TIME_STEP / 2.0);

        // 2. Drift: Update position using the half-step velocity
        body.prevPosition = body.position;
        body.position += body.velocity * PHYSICS_TIME_STEP;

        // Reset netForce for the next step's force accumulation
        body.netForce = Vector2D(0.0, 0.0); 
    }
}

// 2. FORCE/ACCELERATION RECALCULATION
// This function needs to iterate through ALL particles to calculate N^2 gravity
// This calculates the new netForce for the new positions (netForce at t + dt)
// You would call your gravity calculation routine here:
// PhysicsSystem::calculateNetForces(state);

// 3. Kick (Final Half-Step Velocity Update)
void PhysicsSystem::integrateForwards_Phase2(GameState& state)
{   
    auto& particles = state.getParticlesMutable();
    auto& bodies = state.getMacroBodiesMutable();
    for (auto& particle : particles)
    {
        if (particle.isStatic) continue;

        // Calculate new acceleration (a_t + dt) using the newly calculated netForce
        Vector2D new_acceleration = particle.netForce / particle.mass;

        // 3. Kick 2: Update velocity by the remaining half of the acceleration
        particle.velocity += new_acceleration * (PHYSICS_TIME_STEP / 2.0);
    }
    for (auto& body : bodies)
    {
        if (body.isStatic) continue;

        // Calculate new acceleration (a_t + dt) using the newly calculated netForce
        Vector2D new_acceleration = body.netForce / body.mass;

        // 3. Kick 2: Update velocity by the remaining half of the acceleration
        body.velocity += new_acceleration * (PHYSICS_TIME_STEP / 2.0);
    }
}

void PhysicsSystem::calculateGravityForSmallFragments(GravitationalBody& particle, GravitationalBody& body)
{
     if (particle.mass == 0 || body.mass == 0)
    {
        return;
    }

    // Define SOFTENING_EPSILON_SQUARED outside the function as a static const double
    // Let's assume you define: static const double SOFTENING_EPSILON_SQUARED = 10.0;
    static const double G = GRAVITATIONAL_CONSTANT; 
    // static const double epsilonSq = SOFTENING_EPSILON_SQUARED;
    // static const double epsilonSq = 10.0; 
    double epsilon = (particle.radius + body.radius) / 2.0; // Simple average radius
    double epsilonSq = epsilon * epsilon;

    // 1. Calculate Distance Vector
    Vector2D distance = body.position - particle.position;

    // 2. Calculate Distance Squared (r^2)
    double rSq = distance.square_magnitude();

    // 3. Calculate the Denominator Term (r^2 + epsilon^2)^(3/2)
    // The term inside the parenthesis: rSq + epsilonSq
    // The final term in the denominator: pow(rSq + epsilonSq, 1.5)
    double denominator = pow(rSq + epsilonSq, 1.5); 

    // 4. Calculate Force Magnitude (F = G * m1 * m2 / Denominator)
    // double forceMagnitude = (G * particle.mass * body.mass) / denominator;
    double coefficient = (G * particle.mass * body.mass) / denominator;

    // Force vector is C * distance vector (r)
    Vector2D force = distance * coefficient;

    // 6. Apply Forces (Newton's Third Law)
    particle.netForce += force;
    body.netForce -= force;
//    std::cout<<"body1 net force: "<<body1.netForce<<std::endl;
}

void PhysicsSystem::calculateGravity(GravitationalBody& body1, GravitationalBody& body2)
{   
    if (body1.mass == 0 || body2.mass == 0)
    {
        return;
    }

    // Define SOFTENING_EPSILON_SQUARED outside the function as a static const double
    // Let's assume you define: static const double SOFTENING_EPSILON_SQUARED = 10.0;
    static const double G = GRAVITATIONAL_CONSTANT; 
    // static const double epsilonSq = SOFTENING_EPSILON_SQUARED;
    // static const double epsilonSq = 10.0; 
    double epsilon = (body1.radius + body2.radius) / 2.0; // Simple average radius
    double epsilonSq = epsilon * epsilon;

    // 1. Calculate Distance Vector
    Vector2D distance = body2.position - body1.position;

    // 2. Calculate Distance Squared (r^2)
    double rSq = distance.square_magnitude();

    // 3. Calculate the Denominator Term (r^2 + epsilon^2)^(3/2)
    // The term inside the parenthesis: rSq + epsilonSq
    // The final term in the denominator: pow(rSq + epsilonSq, 1.5)
    double denominator = pow(rSq + epsilonSq, 1.5); 

    // 4. Calculate Force Magnitude (F = G * m1 * m2 / Denominator)
    // double forceMagnitude = (G * body1.mass * body2.mass) / denominator;
    double coefficient = (G * body1.mass * body2.mass) / denominator;

    // Force vector is C * distance vector (r)
    Vector2D force = distance * coefficient;

    // 6. Apply Forces (Newton's Third Law)
    body1.netForce += force;
    body2.netForce -= force;
//    std::cout<<"body1 net force: "<<body1.netForce<<std::endl;
}
// void PhysicsSystem::createPlanet(GameState& state, InputState& inputState)
// {
//     GravitationalBody body;
//     body.mass = inputState.selectedMass;
//     body.radius = inputState.selectedRadius;
//     body.position = inputState.mouseCurrPosition;
//     body.prevPosition = body.position;
//     body.isPlanet = true;
//     body.isStatic = inputState.isCreatingStatic;

//     state.getMacroBodiesMutable().push_back(body);
// }
void PhysicsSystem::createPlanet(GameState& state, InputState& inputState)
{
    GravitationalBody body;
    body.mass = inputState.selectedMass;
    body.radius = inputState.selectedRadius;
    body.position = inputState.mouseCurrPosition;
    body.prevPosition = body.position;
    body.isPlanet = true;
    body.isStatic = inputState.isCreatingStatic;

    // --- NEW LOGIC: Nudge fragments out of the new body's radius ---
    auto& particles = state.getParticlesMutable();
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
            particle.prevPosition = particle.position;
        }
    }
    // --- END NEW LOGIC ---

    state.getMacroBodiesMutable().push_back(body);
}

void PhysicsSystem::substituteWithParticles(GravitationalBody& originalBody, GameState& state)
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
                // Store the particle's centerc position
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
        // GravitationalBody radius is tiny, often 0.5 (half a pixel) for collision checks
        p.radius = 1.0; 
        p.position = pos;
        p.prevPosition = pos; // Initialize for Verlet integration
        p.isFragment = true;
        
        // Start with the original body's momentum
        p.velocity = originalVelocity;
        
        // Add a random outward "explosion" vector
        // This is necessary to make the cloud expand.
        Vector2D explosionVector = (pos - center); // Vector from center to particle
        // Normalize and scale by a collision factor (e.g., 0.1)
        explosionVector = explosionVector.normalize() * (0.01 * originalVelocity.magnitude()); 
        
        p.velocity += explosionVector;

        state.getParticlesMutable().push_back(p);
    }

    // 4. CLEANUP
    originalBody.markedForDeletion = true; // Destroy the original macroscopic body
}