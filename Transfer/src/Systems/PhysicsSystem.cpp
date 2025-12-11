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
    // ... (Cluster creation remains the same)
    if (UIState.getInputState().dirty) {
        createNewGravitationalCluster(UIState, state);
    }

    std::vector<Particle>& particles = state.getParticlesMutable();
// 1) LEAPFROG VELOCITY HALF-STEP COMPLETION (v_n)
    // This completes the velocity update v_n = v_{n-1/2} + 0.5 * a_{n-1} * dt
    // using the acceleration a_{n-1} stored in p.prevForce from the last frame.
    for (Particle& p : particles) {
        if (p.mass != 0) {
            Vector2D acc_prev = p.prevForce / p.mass;
            p.velocity += acc_prev * 0.5 * PHYSICS_TIME_STEP; // v_n is now accurate
        }
    } // CRITICAL: After this, p.velocity contains v_n

    // 2) Reset net forces for everyone at frame start (FOR F_n CALCULATION)
    // ... (body reset omitted)
    for (Particle& p : particles) {
        p.netForce = Vector2D(0.0, 0.0); // Reset particle net force
    }

    // 3) Calculate Gravity pairwise (Computes F_n using r_n)
    // ... (Force calculations, including inter-particle and inter-cluster)
    // CRITICAL: After these calls, p.netForce contains F_n
    int size = particles.size();
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = i + 1; j < size; ++j) {
            calculateGravForceBetweenParticles(particles[i], particles[j]);
        }
    }
    // applyInterClusterGravity(state); // NEED TO REWRITE
    // 4) Integrate the system forward (Updates r_{n+1} and v_{n+1/2})
    double totalEnergy = 0.0;
    for (Particle& p : particles) {
        // Use the current F_n to move position to r_{n+1} and velocity to v_{n+1/2}
        integrateParticleLeapfrog(p);
        totalEnergy += 0.5*p.mass*p.velocity.square_magnitude();
    }
    
    // NOTE: In the *next* frame, the final force F_{n+1} will be calculated
    // using the new position r_{n+1}, and then Step 1 will complete the velocity
    // update to v_{n+1} = v_{n+1/2} + 0.5 * a_{n+1} * dt
    // ... (Debugging/collision/load management remain the same)
    handleCollisions(state);
    manageLoad(state);
    std::cout<<"Total system energy"<<totalEnergy<<std::endl;
}
void PhysicsSystem::manageLoad(GameState& state)
{
    std::vector<GravitationalBody>& bodies = state.getBodiesMutable();
    // If we're within limit, nothing to do
    if (bodies.size() <= MAX_BODIES_ALLOWED)
        return;

    // Number of bodies to keep (the largest ones)
    const size_t keep = static_cast<size_t>(MAX_BODIES_ALLOWED);

    // Partition the vector so that the first `keep` elements are the largest masses (unordered among themselves)
    // Using nth_element is O(n) on average and avoids costly sorts or repeated erases.
    std::nth_element(bodies.begin(), bodies.begin() + keep, bodies.end(),
        [](const GravitationalBody& a, const GravitationalBody& b) {
            return a.getMass() > b.getMass(); // descending: larger masses first
        }
    );

    // Erase all elements beyond the `keep` threshold
    bodies.erase(bodies.begin() + keep, bodies.end());
}


void PhysicsSystem::integrateParticleForwards(Particle& p )
{
    
    // F = m a = m dv/dt. F*dt/m = dv. v_new = v + dv. v_new * dt = dx
    Vector2D currentVel = p.velocity;
    Vector2D currentPos = p.position;
    Vector2D currentForce = p.netForce;

    p.prevPosition = currentPos;
    Vector2D dv = Vector2D();
    if (p.mass != 0)
        {
            dv = currentForce * PHYSICS_TIME_STEP / p.mass;
        }
    currentVel += dv;
    currentPos += currentVel * PHYSICS_TIME_STEP;
    p.position = currentPos;
    p.velocity = currentVel;

    // Make call to update the bounding box for the p
    // updateBoundingBox(p);
}
void PhysicsSystem::integrateParticleVelocityVerlet(Particle& p)
{
    // 1. Calculate the acceleration *at the start of the frame* (F_n)
    // You are storing the net force from the *previous* frame in p.prevForce.
    // The force calculated *this* frame is in p.netForce (F_n).
    Vector2D acc_n = p.netForce / p.mass;

    // 2. Save F_n and update position (r_{n+1})
    p.prevForce = p.netForce; // Save F_n for the velocity update later

    p.prevPosition = p.position;
    p.position += p.velocity * PHYSICS_TIME_STEP + acc_n * 0.5 * PHYSICS_TIME_STEP * PHYSICS_TIME_STEP;

    // 3. The new force F_{n+1} must be calculated *after* this position update.
    // In your current structure, this happens in the *next* frame's force calculation loop.
    // The Velocity update (v_{n+1}) must be deferred until the next frame, or you need to
    // restructure the Update loop to calculate F_{n+1} for all particles immediately after
    // calculating r_{n+1}.

    // Since your UpdateSystemFrame is structured (Force Calc -> Integrate), you must
    // implement the half-step of the velocity update here.

    // Half-step velocity update (v_{n+1/2}):
    p.velocity += acc_n * 0.5 * PHYSICS_TIME_STEP;
    
    // 4. Reset netForce for next frame's force accumulation (which will be F_{n+1})
    p.netForce = Vector2D(0.0, 0.0); // Reset for the next frame's force calculation
}

// REPLACE the implementation of the existing integrateParticleLeapfrog
void PhysicsSystem::integrateParticleLeapfrog(Particle& p)
{
    if (p.mass == 0) return;
    double dt = PHYSICS_TIME_STEP;
    
    // F_n is currently in p.netForce, calculated using r_n
    Vector2D acc_n = p.netForce / p.mass; 
    
    // --- STEP 1: Half-Step Velocity (v_{n+1/2}) ---
    // Update velocity using the acceleration calculated this frame (a_n)
    // v_{n+1/2} = v_n + 0.5 * a_n * dt
    p.velocity += acc_n * 0.5 * dt;

    // --- STEP 2: Position Update (r_{n+1}) ---
    // Update position using the intermediate velocity v_{n+1/2}
    // r_{n+1} = r_n + v_{n+1/2} * dt
    p.prevPosition = p.position; // Save r_n
    p.position += p.velocity * dt; // Update to r_{n+1}
    
    // --- STEP 3: Store Force ---
    // Store F_n for the next frame's completion step (Step 1 of next frame)
    p.prevForce = p.netForce; 
    
    // p.netForce will be reset at the start of the next UpdateSystemFrame.
}
void PhysicsSystem::calculateGravForceBetweenParticles(Particle& particleA, Particle& particleB)
{
    double massA = particleA.mass;
    double massB = particleB.mass;
    Vector2D posA = particleA.position;
    Vector2D posB = particleB.position;
    Vector2D deltaDisp = posB - posA;
    double distSqr = deltaDisp.square_magnitude();

    // CRITICAL FIX: Use a small, constant softening parameter (EPSILON)
    // This is the core fix to prevent singularities (instability) when r -> 0.
    // A value of 0.01 to 1.0 is common, depending on coordinate scaling.
    const double EPSILON = 0.9; // TUNEABLE: Adjust this value
    const double EPSILON_SQR = EPSILON * EPSILON;
    
    double effective_dist_sqr = distSqr + EPSILON_SQR;
    
    if (effective_dist_sqr < 1e-12) return; 

    // Calculate the force direction (unit vector)
    Vector2D trigUnitVectors = deltaDisp.normalize();

    // Compute the gravitational force magnitude
    // F = G * m1 * m2 / (r^2 + epsilon^2)
    double force_magnitude = GRAVITATIONAL_CONSTANT * massA * massB / effective_dist_sqr;
    
    // F_B_on_A is the force vector on A due to B (points toward B, so along deltaDisp)
    Vector2D force_B_on_A = trigUnitVectors * force_magnitude; 
    
    Vector2D force_A_on_B = force_B_on_A * (-1.0);
    
    // Accumulate the net forces
    particleA.netForce += force_B_on_A;
    particleB.netForce += force_A_on_B;

    // The function MUST NOT include the 'if (distance > minDistance)' check.
}

void PhysicsSystem::applyInterClusterGravity(GameState& state) {
    auto& clusters = state.getClustersMutable();
    const double G = GRAVITATIONAL_CONSTANT;

    // 0) Recompute COMs for this frame
    for (size_t ci = 0; ci < clusters.size(); ++ci) {
        recomputeClusterProperties(state, ci);
    }

    // 1) Loop over all cluster pairs
    for (size_t i = 0; i < clusters.size(); ++i) {
        for (size_t j = i + 1; j < clusters.size(); ++j) {
            auto& clusterA = clusters[i];
            auto& clusterB = clusters[j];

            Vector2D delta = clusterB.centerOfMassPos - clusterA.centerOfMassPos;
            double distSqr = delta.square_magnitude();
            double dist = sqrt(distSqr);

            if (dist < 1e-6) continue; // Avoid divide by zero

            Vector2D direction = delta / dist; // normalized

            // Compute the gravitational force magnitude
            double forceMag = G * clusterA.totalMass * clusterB.totalMass / distSqr;
            Vector2D forceOnA = direction * forceMag;  // Force on cluster A
            Vector2D forceOnB = forceOnA * (-1.0);     // Newton's 3rd law

            // --- NEW: Apply at COM level ---
            Vector2D accA = forceOnA / clusterA.totalMass;
            Vector2D accB = forceOnB / clusterB.totalMass;

            for (int pIndex : clusterA.clusterParticleIndices) {
                Particle& p = state.getParticlesMutable()[pIndex];
                // Convert acceleration back to force F = ma, or just use forceOnA/totalMass * p.mass
                // Or simply: F = m*a
                Vector2D force_on_particle = accA * p.mass;
                p.netForce += force_on_particle; 
                
                // --- OLD: Remove this line entirely ---
                // p.velocity += accA * PHYSICS_TIME_STEP;
                // --- END OLD ---
            }
            for (int pIndex : clusterB.clusterParticleIndices) {
                Particle& p = state.getParticlesMutable()[pIndex];
                Vector2D force_on_particle = accB * p.mass;
                p.netForce += force_on_particle;

                // --- OLD: Remove this line entirely ---
                // p.velocity += accB * PHYSICS_TIME_STEP;
                // --- END OLD ---
            }
        }
    }
}
void PhysicsSystem::handleCollisions(GameState& state)
{   
    std::vector<Particle>& particles = state.getParticlesMutable();
    int size = state.getParticles().size();

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = i + 1; j < size; ++j) {

            Particle& particleA = particles[i];
            Particle& particleB = particles[j];
            handleElasticParticleCollision(particleA, particleB);
        }
    }
}


void PhysicsSystem::handleElasticParticleCollision(Particle& A, Particle& B)
{
    Vector2D x = B.position - A.position;
    double dist = x.magnitude();

    if (dist == 0) return;
    Vector2D n = x / dist; // collision normal

    Vector2D vA = A.velocity;
    Vector2D vB = B.velocity;

    double mA = A.mass;
    double mB = B.mass;

    // Project velocities onto normal
    double vA_n = vA.dot(n);
    double vB_n = vB.dot(n);

    // 1D elastic collision formula (normal direction only)
    double vA_n_new = (vA_n*(mA - mB) + 2*mB*vB_n) / (mA + mB);
    double vB_n_new = (vB_n*(mB - mA) + 2*mA*vA_n) / (mA + mB);

    // Change in normal components
    Vector2D vA_change = n * (vA_n_new - vA_n)*0.9;
    Vector2D vB_change = n * (vB_n_new - vB_n)*0.9;

    // Apply
    A.velocity = vA + vA_change;
    B.velocity = vB + vB_change;

    // --- Positional correction to prevent re-collision ---
    double overlap = A.radius + B.radius - dist;
    if (overlap > 0) {
        A.position = A.position - n * (overlap * 0.5);
        B.position = B.position + n * (overlap * 0.5);
    }
}



// void PhysicsSystem::createNewGravitationalCluster(UIState& UIState, GameState& state)
// {
//     // 1. Initial setup and validation
//     auto& inputState = UIState.getMutableInputState();
//     if (!inputState.isCreatingCluster) {
//         return;
//     }
    
//     // Create the cluster object and add it to the state
//     GravitationalCluster cluster;
//     cluster.centerOfMassPos = inputState.mouseCurrPosition;
//     cluster.totalMass = inputState.selectedMass;
//     cluster.initialRadius = inputState.selectedRadius; 

//     int clusterIndex = state.addCluster(cluster);
//     auto& newCluster = state.getClustersMutable()[clusterIndex];

//     int numParticles = 100; // Hardcoded particle count
//     double particleMass = cluster.totalMass / numParticles;

//     // 2. Initialize particles with positions and zero velocity
//     newCluster.clusterParticleIndices.reserve(numParticles);
//     for (int i = 0; i < numParticles; ++i) {
//         Particle p;
        
//         // Random position within a circle (uniform area distribution)
//         double angle = (double)rand() / RAND_MAX * 2.0 * PI;
//         double r = sqrt((double)rand() / RAND_MAX) * cluster.initialRadius;

//         p.position = newCluster.centerOfMassPos + Vector2D(r*cos(angle), r*sin(angle));
        
//         // Temporarily set velocity and previous position to zero/equal
//         p.prevPosition = p.position; 
//         p.velocity = Vector2D(0.0, 0.0); 
        
//         p.mass = particleMass;
//         p.clusterID = clusterIndex;
//         p.radius = r; // Particle radius currently equals distance from center (r), should be a small constant or calculated from mass/density
//         p.netForce = Vector2D(0.0, 0.0);
//         p.prevForce = Vector2D(0.0, 0.0); // Required for Velocity Verlet
        
//         // Add particle to state and track its index in the cluster
//         int particleIndex = state.addParticle(p);
//         newCluster.clusterParticleIndices.push_back(particleIndex);
//     }
    
//     // 3. Recompute COM (essential to get the exact COM and totalMass for the next step)
//     recomputeClusterProperties(state, clusterIndex);
    
//     // // 4. Apply Initial Keplerian Velocity Kick for stability
//     // for (int pIndex : newCluster.clusterParticleIndices) {
//     //     Particle& p = state.getParticlesMutable()[pIndex];
//     //     Vector2D r_vec = p.position - newCluster.centerOfMassPos;
//     //     double r = r_vec.magnitude();

//     //     // Keplerian velocity magnitude: V = sqrt(G*M/r)
//     //     double V_mag = 0.0;
//     //     if (r > 1e-6) {
//     //         // Using a square root term to mimic a smooth distribution rather than V=sqrt(M/r)
//     //         V_mag = sqrt(GRAVITATIONAL_CONSTANT * newCluster.totalMass / r); 
//     //     }
        
//     //     // V must be perpendicular to r_vec (tangential direction)
//     //     Vector2D r_norm = r_vec.normalize();
//     //     Vector2D tangent_norm = Vector2D(-r_norm.y_val, r_norm.x_val); 

//     //     // Apply rotational velocity. Scaling by 0.5 allows the system to contract slightly 
//     //     // toward a stable state, preventing immediate dispersal.
//     //     p.velocity = tangent_norm * V_mag * 0.5; 

//     //     // Set p.prevPosition to satisfy the first time step of Velocity Verlet/Leapfrog.
//     //     // r_{n-1} = r_n - v_n * dt
//     //     p.prevPosition = p.position - p.velocity * PHYSICS_TIME_STEP; 
//     // }
//     // 4. Apply Initial Velocity (Simplified)
//     for (int pIndex : newCluster.clusterParticleIndices) {
//         Particle& p = state.getParticlesMutable()[pIndex];
//         Vector2D r_vec = p.position - newCluster.centerOfMassPos;
//         double r = r_vec.magnitude();

//         // --- FIX: Simplify Initial Velocity ---
//         // A stable cluster needs rotation proportional to distance (solid body rotation)
//         // or zero velocity ("cold start"). A cold start will lead to collapse/oscillation.
//         // A mild solid-body rotation is more stable than Keplerian.
        
//         // TUNEABLE: Angular speed (e.g., 0.005)
//         const double OMEGA = 0.005; 
        
//         if (r > 1e-6) {
//             // v = omega * r. Velocity must be tangential.
//             Vector2D tangent_norm = Vector2D(-r_vec.y_val, r_vec.x_val).normalize(); 
//             p.velocity = tangent_norm * (r * OMEGA); 
//         } else {
//             p.velocity = Vector2D(0.0, 0.0);
//         }
        
//         // Set p.prevPosition to satisfy the first time step of Leapfrog/Velocity Verlet.
//         // r_{n-1} = r_n - v_n * dt
//         p.prevPosition = p.position - p.velocity * PHYSICS_TIME_STEP; 
//     }
//     // 5. Reset flags
//     inputState.isCreatingCluster = false;
//     inputState.dirty = false;
// }
void PhysicsSystem::createNewGravitationalCluster(UIState& UIState, GameState& state)
{
    // 1. Initial setup and validation (Unchanged)
    auto& inputState = UIState.getMutableInputState();
    if (!inputState.isCreatingCluster) {
        return;
    }
    
    // Create the cluster object and add it to the state
    GravitationalCluster cluster;
    cluster.centerOfMassPos = inputState.mouseCurrPosition;
    cluster.totalMass = inputState.selectedMass;
    cluster.initialRadius = inputState.selectedRadius; 

    int clusterIndex = state.addCluster(cluster);
    auto& newCluster = state.getClustersMutable()[clusterIndex];

    int numParticles = 100;
    double particleMass = cluster.totalMass / numParticles;

    // --- CRITICAL CONSTANT: Particle Radius Calculation ---
    // Assuming uniform density (D = M/V), particle radius should be small and constant
    // to avoid immediate overlap in the dense core, but large enough for collisions.
    // Use a fixed small radius, or one based on the total cluster volume.
    // For simplicity and stability, use a small fraction of the cluster's initial radius.
    const double PARTICLE_RADIUS = cluster.initialRadius * 0.10; // TUNEABLE: 5% of cluster radius

    // 2. Initialize particles with positions and initial velocity
    newCluster.clusterParticleIndices.reserve(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        Particle p;
        
        // Random position within a circle (uniform area distribution)
        double angle = (double)rand() / RAND_MAX * 2.0 * PI;
        double r = sqrt((double)rand() / RAND_MAX) * cluster.initialRadius;

        p.position = newCluster.centerOfMassPos + Vector2D(r*cos(angle), r*sin(angle));
        
        // Temporarily set velocity to zero/equal
        p.velocity = Vector2D(0.0, 0.0); 
        
        p.mass = particleMass;
        p.clusterID = clusterIndex;
        // FIX: Set a consistent particle radius for softening and collisions
        p.radius = PARTICLE_RADIUS; 
        
        p.netForce = Vector2D(0.0, 0.0);
        p.prevForce = Vector2D(0.0, 0.0); 
        
        // Add particle to state and track its index in the cluster
        int particleIndex = state.addParticle(p);
        newCluster.clusterParticleIndices.push_back(particleIndex);
    }
    
    // 3. Recompute COM (essential to get the exact COM and totalMass for the next step)
    recomputeClusterProperties(state, clusterIndex);
    
    // 4. Apply Simplified Rotational Kick (Solid Body Rotation) for stability
    for (int pIndex : newCluster.clusterParticleIndices) {
        Particle& p = state.getParticlesMutable()[pIndex];
        Vector2D r_vec = p.position - newCluster.centerOfMassPos;
        double r = r_vec.magnitude();

        // TUNEABLE: Angular speed (e.g., 0.005) - Controls rotation rate
        const double OMEGA = 0.005; 
        
        if (r > 1e-6) {
            // Calculate tangential velocity: V = omega * r
            Vector2D tangent_norm = Vector2D(-r_vec.y_val, r_vec.x_val).normalize(); 
            p.velocity = tangent_norm * (r * OMEGA); 
        } else {
            p.velocity = Vector2D(0.0, 0.0);
        }
        
        // CRITICAL LEAPFROG INIT STEP: Set p.prevPosition (r_{n-1})
        // The integrator starts by assuming p.position is r_n and p.velocity is v_n.
        // The first frame's Step 1 (velocity half-step completion) is skipped.
        // For the Position Leapfrog integrator to work, it needs $r_{n}$ and $v_{n}$ to calculate $r_{n+1}$.
        // If we set $r_{n-1}$ such that $r_n - r_{n-1} = v_n \Delta t$, the integrator has a clean start.
        // r_{n-1} = r_n - v_n * dt
        p.prevPosition = p.position - p.velocity * PHYSICS_TIME_STEP; 
    }

    // 5. Reset flags (Unchanged)
    inputState.isCreatingCluster = false;
    inputState.dirty = false;
}

void PhysicsSystem::recomputeClusterProperties(GameState& state, int clusterIndex)
{
    auto& cluster = state.getClustersMutable()[clusterIndex];
    const auto& particles = state.getParticles();

    Vector2D comPosition(0.0, 0.0);
    Vector2D comVelocity(0.0, 0.0);
    double totalMass = 0.0;

    for (int pIndex : cluster.clusterParticleIndices){
        const Particle& p = particles[pIndex];
        // Skip deleted or bad indices if they exist
        if (p.markedForDeletion) continue;
        comPosition += p.position * p.mass;
        comVelocity += p.velocity * p.mass;
        totalMass += p.mass;
    }

    if (totalMass > 0.0){
        comPosition = comPosition / totalMass;
        comVelocity = comVelocity / totalMass;
    }

    cluster.centerOfMassPos = comPosition;
    cluster.bulkVelocity = comVelocity;
    cluster.totalMass = totalMass;
}

