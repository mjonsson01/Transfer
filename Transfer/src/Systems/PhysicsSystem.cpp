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
                inputState.resetTransientFlags();
            }
        if (inputState.isCreatingDust)
            {
                createDust(state, inputState);
                inputState.resetTransientFlags();
            }
    }

    if (inputState.clearAll){
        state.getMacroBodiesMutable().clear();
        state.getParticlesMutable().clear();
        inputState.clearAll = false;
    }
    if (!inputState.isPaused){
            integrateForwards_Phase1(state);
        
        auto& particles = state.getParticlesMutable();
        int num_particles = particles.size();
        // std::cout<<"num particles"<<num_particles<<std::endl;

        auto& bodies = state.getMacroBodiesMutable();
        int num_bodies = bodies.size();
        // std::cout<<"num bodies: "<<num_bodies<<", num particles: " << num_particles<<std::endl;
        // if (num_bodies > 0)
        //     std::cout <<"Body Speed"<<2*bodies[0].velocity.magnitude()<<std::endl; 

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
    Vector2D vP_change = n * (vP_n_new - vP_n)*ELASTIC_LOSS_FACTOR;
    Vector2D vB_change = n * (vB_n_new - vB_n)*ELASTIC_LOSS_FACTOR;

//     // Apply
    particle.velocity = vP + vP_change;
    body.velocity = vB + vB_change;
    
    double overlap = particle.radius + body.radius - dist;

    // std::cout<<"Grav Body:\n"<<body<<std::endl;
    // std::cout<<"Particle:\n"<<particle<<std::endl;

    if (abs(overlap) > 0){
        // Small slop to avoid micro jitter
        const double slop = 0.01;
        const double percent = 0.8; // correction strength

        double correctionMag = std::max(overlap - slop, 0.0) * percent;
        Vector2D correction = n * correctionMag;

        double totalMass = particle.mass + body.mass;

        particle.position -= correction * (body.mass / totalMass);
        body.position += correction * (particle.mass / totalMass);

        // Critical: prevent fake velocity from Verlet-style position jumps
        particle.prevPosition = particle.position;
        body.prevPosition = body.position;
    }
    // std::cout<<"overlap: "<<overlap<<std::endl;
    // if (overlap > 1.0) // more than a pixel overlap
    // {
    //     // Vector2D separation_velocity = n * (overlap / PHYSICS_TIME_STEP);
    //     // particle.velocity -= separation_velocity * 0.5; // push particle away
    //     // body.velocity += separation_velocity * 0.5; // push body away
    // }

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
    // double overlap = particle.radius + body.radius - dist;
    
    // if (overlap > 0.0) {
    //     // Calculate the required velocity (V_sep) to correct the overlap (P) over 1 timestep (dt)
    //     // V_sep = P / dt. The factor of 1.05 gives a slight margin of separation.
    //     double sep_velocity_mag = (overlap / PHYSICS_TIME_STEP) * 0.07; 
        
    //     Vector2D separation_velocity = n * sep_velocity_mag; 

    //     // Apply the separation velocity, weighted by mass
    //     double mSum = particle.mass + body.mass;
        
    //     // Push particle away from the body
    //     particle.velocity -= separation_velocity * (body.mass / mSum); 
    //     // Push body away from the particle
    //     body.velocity += separation_velocity * (particle.mass / mSum); 
    // }
}

void PhysicsSystem::handleDynamicExplosionCollision(GravitationalBody& body1, GravitationalBody& body2, GameState& state)
{

    // end with substituting both with particles
    // substituteWithParticles(body1, state);
    // substituteWithParticles(body2, state);
      // A simplified value for the collision radius (R)
    const double R1 = body1.radius;
    const double R2 = body2.radius;
    const Vector2D center1 = body1.position;
    const Vector2D center2 = body2.position;
    const Vector2D center_of_explosion = (center1 + center2)/ 2.0;
    const double originalMass1 = body1.mass;
    const double originalMass2 = body2.mass;
    const Vector2D originalVelocity1 = body1.velocity;
    const Vector2D originalVelocity2 = body2.velocity;
    
    std::vector<Vector2D> pixelPositions1; 

    // 1. RASTERIZATION (Find all pixel positions)
    // Assuming a 1-to-1 mapping where 1 unit = 1 pixel for simplicity.
    int R_int1 = static_cast<int>(R1);
    for (int i = -R_int1; i <= R_int1; ++i) {
        for (int j = -R_int1; j <= R_int1; ++j) {
            // Check if (i, j) is inside the circle
            if ((i * i) + (j * j) <= (R1 * R1)) {
                // Store the particle's center position
                pixelPositions1.push_back(center1 + Vector2D{(double)i, (double)j});
            }
        }
    }
    std::vector<Vector2D> pixelPositions2; 
    int R_int2 = static_cast<int>(R2);
    for (int i = -R_int2; i <= R_int2; ++i) {
        for (int j = -R_int2; j <= R_int2; ++j) {
            // Check if (i, j) is inside the circle
            if ((i * i) + (j * j) <= (R2 * R2)) {
                // Store the particle's center position
                pixelPositions2.push_back(center2 + Vector2D{(double)i, (double)j});
            }
        }
    }
    
    // Check if we found any pixels (safety)
    if (pixelPositions1.empty()){
        // std::cout<<"no pixels"<<std::endl;
        return;
    } 
        // Check if we found any pixels (safety)
    if (pixelPositions2.empty()){
        // std::cout<<"no pixels"<<std::endl;
        return;
    } 

    // 2. MASS CALCULATION
    const double particleMass1 = originalMass1 / pixelPositions1.size();
    const double particleMass2 = originalMass2 / pixelPositions2.size();

    // 3. PARTICLE INSTANTIATION
    for (const auto& pos : pixelPositions1) {
        GravitationalBody p;
        p.mass = particleMass1;
        p.radius = 1.0; 
        p.position = pos;
        p.isFragment = true;
        
        Vector2D r = pos - center_of_explosion;
        double dist = r.magnitude();

        Vector2D n = (dist > 1e-8) ? (r / dist) : Vector2D(1, 0);
        Vector2D t = Vector2D(-n.y_val, n.x_val);

        // Pick a characteristic radius for falloff.
        // This makes the falloff scale with the explosion size.
        double falloffRadius = 0.5 * (R1 + R2);          // or tweak: 0.25*(R1+R2), etc.
        double d = dist / (falloffRadius + 1e-8);

        // falloff in (0,1], where 1 at center, ~0 far away
        double falloff = 1.0 / (1.0 + d * d);            // nice smooth curve

        double factor = originalVelocity1.magnitude();
        double blast = factor * falloff;
        double shear = factor * falloff;

        double jitter = randomDouble(-1.0, 1.0) * falloff;

        // p.velocity = lerp(originalVelocity2, originalVelocity1, 0.5) + n * (blast + jitter) + t * (shear * randomDouble(-1.0, 1.0));
        Vector2D originalVel1 = originalVelocity1; // or originalVelocity2 depending on loop

        // Near center: use blended/mixed base. Far: keep original body velocity.
        Vector2D mixedBase1 = lerp(originalVelocity2, originalVelocity1, 0.5);
        Vector2D baseVel1   = lerp(originalVel1, mixedBase1, falloff);
        
        p.velocity = baseVel1 + n * (blast + jitter) + t * (shear * randomDouble(-1.0, 1.0));


        p.prevPosition = p.position - p.velocity * PHYSICS_TIME_STEP;

        state.getParticlesMutable().push_back(p);
    }
    for (const auto& pos : pixelPositions2) {
        GravitationalBody p;
        p.mass = particleMass2;
        p.radius = 1.0; 
        p.position = pos;
        p.isFragment = true;
        
        Vector2D r = pos - center_of_explosion;
        double dist = r.magnitude();

        Vector2D n = (dist > 1e-8) ? (r / dist) : Vector2D(1, 0);
        Vector2D t = Vector2D(-n.y_val, n.x_val);

        // Pick a characteristic radius for falloff.
        // This makes the falloff scale with the explosion size.
        double falloffRadius = 0.5 * (R1 + R2);          // or tweak: 0.25*(R1+R2), etc.
        double d = dist / (falloffRadius + 1e-8);

        // falloff in (0,1], where 1 at center, ~0 far away
        double falloff = 1.0 / (1.0 + d * d);            // nice smooth curve
        double factor = originalVelocity2.magnitude();
        double blast = factor * falloff;
        double shear = factor * falloff;

        double jitter = randomDouble(-1.0, 1.0) * falloff;

        // p.velocity = lerp(originalVelocity1, originalVelocity2, 0.5) + n * (blast + jitter) + t * (shear * randomDouble(-1.0, 1.0));
        Vector2D originalVel2 = originalVelocity2; // or originalVelocity2 depending on loop

        // Near center: use blended/mixed base. Far: keep original body velocity.
        Vector2D mixedBase2 = lerp(originalVelocity2, originalVelocity1, 0.5);
        Vector2D baseVel2   = lerp(originalVel2, mixedBase2, falloff);
        p.velocity = baseVel2 + n * (blast + jitter) + t * (shear * randomDouble(-1.0, 1.0));

        p.prevPosition = p.position - p.velocity * PHYSICS_TIME_STEP;
        
        
        state.getParticlesMutable().push_back(p);
    }

    // 4. CLEANUP
    body1.markedForDeletion = true; // Destroy the original macroscopic body
    body2.markedForDeletion = true;
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
        if (bodies[i].markedForDeletion) continue;
        for (size_t j = i + 1; j < bodies.size(); ++j)
        {
            if (bodies[j].markedForDeletion) continue;
            double distance = (bodies[i].position - bodies[j].position).magnitude();
            // double net_speed = (bodies[i].velocity - bodies[j].velocity).magnitude();
            Vector2D r_vector = (bodies[i].position - bodies[j].position);
            Vector2D r_unit_vector = r_vector.normalize();
            Vector2D net_velocity= (bodies[i].velocity - bodies[j].velocity);
            double normal_speed = net_velocity.dot(r_unit_vector);
            double abs_normal_speed = std::abs(normal_speed);
            // add helper function to return more massive body and get rid of this case check
            int case_check = 0;
            if (bodies[i].mass > bodies[j].mass){
                case_check = 1; // mass 1 GREATER
            }
            else if (bodies[i].mass < bodies[j].mass){
                case_check = 2; // mass 2 GREATER
            }
            else
                case_check = 3; // equal
            double mass_ratio = 0; // prob rewrite as helper here too

            if (case_check == 0)
                break;
            else if (case_check == 1)
                mass_ratio = bodies[i].mass / bodies[j].mass;
            else if (case_check == 2)
                mass_ratio = bodies[j].mass / bodies[i].mass;
            else
                mass_ratio = 1.0;



            if (distance < (bodies[i].radius + bodies[j].radius)){
                if (abs_normal_speed < MIN_SHATTER_SPEED){ //not in shatter speed territory
                    if (mass_ratio >= MIN_BODY_BODY_ACCRETION_THRESHOLD_RATIO) // if ratio is bigger than min body accretion threshold ratio (one body dwarfs the other by enough)
                        {
                            if (case_check == 1)
                                handleAccretion(bodies[j],bodies[i]);// j is lighter
                            else
                                handleAccretion(bodies[i],bodies[j]);// i is lighter
                        }
                    else if (mass_ratio < MIN_BODY_BODY_ACCRETION_THRESHOLD_RATIO)
                        // refactor again later
                        if (case_check == 1)
                            handleElasticCollisions(bodies[j], bodies[i]); // no shatter logic yet, just bounce the bodies off each other
                        else
                            handleElasticCollisions(bodies[i], bodies[j]); // no shatter logic yet, just bounce the bodies off each other
                }
                else // in shatter speed territory
                    if (case_check == 1)
                        substituteWithParticles(bodies[j], state);
                    else if (case_check == 2)
                        substituteWithParticles(bodies[i], state);
                    else {
                        // handleElasticCollisions(bodies[i], bodies[j]); // equal mass, just bounce off each other but immediately substitute with particles
                        // handleDynamicExplosionCollision(bodies[i], bodies[j], state);
                        substituteWithParticles(bodies[i], state);
                        substituteWithParticles(bodies[j], state);
                    }
                }
                //     if (net_speed < MIN_SHATTER_SPEED)
                //     {
                //         handleElasticCollisions(bodies[i], bodies[j]); 
                //     }
                //     else
                //     {
                //         // shatter more massive body
                //         if (bodies[i].mass > bodies[j].mass)
                //             substituteWithParticles(bodies[j], state);
                //         else if (bodies[j].mass > bodies[i].mass)
                //             substituteWithParticles(bodies[i], state);
                //         else {
                //             substituteWithParticles(bodies[j], state);
                //             substituteWithParticles(bodies[i], state);
                //         }    
                //     }
        }
    }
    auto& particles = state.getParticlesMutable();
    int num_particles = particles.size();

    for (auto& particle : particles)
    {
        if (particle.markedForDeletion) continue;
        for (auto& body : bodies){
            if (body.markedForDeletion) continue;
            double distance = (particle.position - body.position).magnitude();
            
            if (distance < (particle.radius + body.radius))
            {
                // Call the function containing the fixed logic
                // double net_speed = (particle.velocity - body.velocity).magnitude();
                Vector2D r_vector = (particle.position - body.position);
                Vector2D r_unit_vector = r_vector.normalize();
                Vector2D net_velocity= (particle.velocity - body.velocity);
                // double normal_velocity = net_velocity.dot(r_unit_vector);
                // double net_normal_speed = (normal_velocity).magnitude();
                double normal_speed = net_velocity.dot(r_unit_vector);
                double abs_normal_speed = std::abs(normal_speed);

                if (abs_normal_speed > MAX_ACCRETION_COLLISION_SPEED)// && body.mass / particle.mass > MIN_BODY_PARTICLE_ACCRETION_THRESHOLD_RATIO) //
                    {
                        // std::cout<<"elastic"<<std::endl;
                        // handleElasticCollisions(particle, body);
                        handleAccretion(particle, body);
                    }
                else 
                    handleAccretion(particle, body);
            }
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

void PhysicsSystem::handleAccretion(GravitationalBody& particle, GravitationalBody& body)
{
    // Accrete mass
    body.mass += particle.mass;
    // Add to radius
    body.radius = body.radius * pow((body.mass + particle.mass) / body.mass, 1.0/3.0); // investigate caching this calculation?
    particle.markedForDeletion = true;
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

    state.getMacroBodiesMutable().push_back(body);
}

void PhysicsSystem::createDust(GameState& state, InputState& inputState)
{
    GravitationalBody body;
    body.mass = inputState.selectedMass;
    body.radius = inputState.selectedRadius;
    body.position = inputState.mouseCurrPosition;
    body.prevPosition = body.position;
    body.isDust = inputState.isCreatingDust;
    body.isStatic = inputState.isCreatingStatic;
    inputState.isCreatingDust = false;
    state.getParticlesMutable().push_back(body);
    inputState.dirty = false;
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
        p.prevPosition = pos; // Initialize for Verlet integration
        p.isFragment = true;
        
        // Start with the original body's momentum
        p.velocity = originalVelocity;
        
        // Add a random outward "explosion" vector
        // This is necessary to make the cloud expand.
        Vector2D explosionVector = (pos - center); // Vector from center to particle
        // Normalize and scale by a collision factor (e.g., 0.1)
        explosionVector = explosionVector.normalize() * (0 * originalVelocity.magnitude()); 
        
        p.velocity += explosionVector;
        p.prevPosition = p.position - p.velocity * PHYSICS_TIME_STEP;


        state.getParticlesMutable().push_back(p);
    }

    // 4. CLEANUP
    originalBody.markedForDeletion = true; // Destroy the original macroscopic body
}