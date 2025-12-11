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
void PhysicsSystem::UpdateSystemFrame(GameState& state)
{
    // Update physics calculations for all entities in the game state

    // Extract body vector as modifiable
    std::vector<GravitationalBody>& bodies = state.getBodiesMutable();

    // 1) Reset net forces for everyone at frame start
    for (size_t i = 0; i < bodies.size(); ++i) {
        bodies[i].setNetForce(Vector2D(0.0, 0.0));
    }

    // Calculate Gravity pairwise
	for (size_t i = 0; i < bodies.size(); ++i) {
		for (size_t j = i + 1; j < bodies.size(); ++j) {
            calculateGravForceBetweenBodies(bodies[i], bodies[j]);
		}
	}


    // integrate the system forward after all the resultant forces are calculated 
    for (size_t i = 0; i < bodies.size(); ++i){
        integrateForwards(bodies[i]);
    }

    // if (bodies.size() >= 2){
    //     // Move cursor up 5 lines, clear those lines, print
    //     // std::cout << "\033[11A\033[J" << "\nBody 0 status\n" << bodies[0] << std::flush;
    //     // std::cout << "\033[1A\033[J" <<" Body 0 speed: " << bodies[0].getNetVelocity().magnitude() << std::flush;
    //     std::cout << "\r\033[K" << "Body 0: speed=" << bodies[0].getNetVelocity().magnitude() << std::flush;
    // }   
    // for (auto& b : bodies){
        //     b.updateGhostState();
    // }
        
    // handle collisions and update velocities of anything detected as colliding. later
    handleCollisions(state);

    manageLoad(state);
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

void PhysicsSystem::integrateForwards(GravitationalBody& body )
{
    
    // F = m a = m dv/dt. F*dt/m = dv. v_new = v + dv. v_new * dt = dx
    Vector2D currentVel = body.getNetVelocity();
    Vector2D currentPos = body.getPosition();
    Vector2D currentForce = body.getNetForce();

    body.setPrevPosition(currentPos);
    Vector2D dv = Vector2D();
    if (body.getMass() != 0)
        {
            dv = currentForce * PHYSICS_TIME_STEP / body.getMass();
        }
    currentVel += dv;
    currentPos += currentVel * PHYSICS_TIME_STEP;
    body.setPosition(currentPos);
    body.setNetVelocity(currentVel);

    // Make call to update the bounding box for the body
    updateBoundingBox(body);
}

void PhysicsSystem::calculateGravForceBetweenBodies(GravitationalBody& bodyA, GravitationalBody& bodyB)
{
    // Calculate gravitational force between bodyA and bodyB
    double massA = bodyA.getMass();
    double massB = bodyB.getMass();

    // Find distance between bodies
    Vector2D posA = bodyA.getPosition();
    Vector2D posB = bodyB.getPosition();
    Vector2D deltaDisp = posB - posA;
    double distance = deltaDisp.magnitude();
    double minDistance = bodyA.getRadius() + bodyB.getRadius();

    if (distance > minDistance){ 
        // sin(theta),cos(theta) unit vector
        Vector2D trigUnitVectors = deltaDisp.normalize();
        double force_magnitude = GRAVITATIONAL_CONSTANT * massA * massB / deltaDisp.square_magnitude();
        Vector2D force_B_on_A = Vector2D(force_magnitude * trigUnitVectors.x_val, force_magnitude * trigUnitVectors.y_val);
        // Newton's Third Law
        Vector2D force_A_on_B = force_B_on_A * (-1.0);
        bodyA.addNetForce(force_B_on_A);
        bodyB.addNetForce(force_A_on_B);
    }
    else{
        return;
    }
}


// MINE
void PhysicsSystem::handleCollisions(GameState& state)
{
    std::vector<GravitationalBody>& bodies = state.getBodiesMutable();
    int size = state.getBodies().size();

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = i + 1; j < size; ++j) {

            GravitationalBody& bodyA = bodies[i];
            GravitationalBody& bodyB = bodies[j];

            
            // Skip deleted bodies
            if (bodyA.getMarkedForDeletion() || bodyB.getMarkedForDeletion()){
                continue;
            }
            // // Skip if collision is disabled on either body
            if (!bodyA.getCollisionEnabled() || !bodyB.getCollisionEnabled())
                continue;

            // Check bounding boxes to reduce search space. If bounding boxes don't overlap, don't even bother. Bounding box method limited if extreme acceleration is in play.
            if (!bodyA.getBoundingBox().overlaps(bodyB.getBoundingBox()))
                continue; // skip


            // Actually deal with the collisions
            
            // Precise circle check (may include real overlap in case of clicking too close, which results in accretion or correction of positions)
            Vector2D direction_vector = bodyB.getPosition() - bodyA.getPosition();
            double dist = direction_vector.magnitude();
            double min_distance = bodyA.getRadius() + bodyB.getRadius();
            
            // Identify the smaller and larger body radius-wise
            // GravitationalBody* smallerR_body;
            // GravitationalBody* largerR_body;

            // // Identify the smaller and larger body mass-wise
            // GravitationalBody* smallerM_body;
            // GravitationalBody* largerM_body;

            // if (bodyA.getRadius() < bodyB.getRadius()) {
            //     smallerR_body = &bodyA;
            //     largerR_body  = &bodyB;
            // } else {
            //     smallerR_body = &bodyB;
            //     largerR_body  = &bodyA;
            // }
            
            // if (bodyA.getMass() < bodyB.getMass()) {
            //     smallerM_body = &bodyA;
            //     largerM_body  = &bodyB;
            // } else {
            //     smallerM_body = &bodyB;
            //     largerM_body  = &bodyA;
            // }
            // double mass_ratio = largerM_body->getMass() / smallerM_body->getMass();

            // if (dist < min_distance && dist > 0.0) {
            //     double relativeSpeed = (bodyA.getNetVelocity() - bodyB.getNetVelocity()).magnitude();
                
            //     handleDynamicCollision(bodyA, bodyB, state);
                
            //     // if (relativeSpeed < MAX_ELASTIC_COLLISION_SPEED){
            //         // handleElasticCollision(bodyA, bodyB);
            //     // }
            //     // else{
            //         // handleDynamicCollision(bodyA, bodyB, state);
            //         // handleElasticCollision(bodyA, bodyB);
            //     // }
            // }
            // else {
            //     // No collision
            //     continue;
            // }

            // If the mass ratio is small enough just accrete.

            // Figure out the mutual escape velocity to determine the threshold for accretion vs fragmentation.
    }
    bodies.erase(
        std::remove_if(bodies.begin(), bodies.end(),
                    [](const GravitationalBody& b){ return b.getMarkedForDeletion(); }),
        bodies.end()
    );

}

void PhysicsSystem::handleDynamicCollision(GravitationalBody& bodyA, GravitationalBody& bodyB, GameState& state)
{
    // Placeholder for dynamic collision handling logic
    // Currently, we just mark both bodies for deletion to avoid compilation errors.
    bodyA.setMarkedForDeletion(true);
    bodyB.setMarkedForDeletion(true);
}


// void PhysicsSystem::handleDynamicCollision(GravitationalBody& bodyA, GravitationalBody& bodyB, GameState& state)
// {
//     // --- 1. Extract Initial State and Define Target Momentum ---    
//     double mA = bodyA.getMass();
//     double mB = bodyB.getMass();
//     Vector2D vA = bodyA.getNetVelocity();
//     Vector2D vB = bodyB.getNetVelocity();
    
//     // Calculate total mass
//     double M = mA + mB;

//     // Calculate the total system momentum (MUST BE CONSERVED)
//     Vector2D P_initial = vA * mA + vB * mB;

//     // Calculate the velocity of the Center of Mass (CoM). 
//     // This is the target base velocity for all fragments.
//     Vector2D V_CoM = P_initial / M;
    
//     // Calculate the initial kinetic energy of the system for reference (used later to size the fragments/scatter)
//     // The relative velocity determines the intensity of the fragmentation.
//     Vector2D deltaV = vA - vB;
//     double relativeSpeed = deltaV.magnitude();


//     // --- 2. Determine Fragmentation Outcome and Generate Fragment Masses ---

//     // Define tuneable constants for fragmentation
//     const double FRAGMENTATION_RATIO = 5.0; // If mass ratio exceeds this, only the smaller body fragments.
//     const int MIN_FRAGMENTS = 8;
//     const int MAX_FRAGMENTS = 16;
    
//     // Pointers for clarity
//     GravitationalBody* smallerM_body = (mA < mB) ? &bodyA : &bodyB;
//     GravitationalBody* largerM_body = (mA < mB) ? &bodyB : &bodyA;
//     double mass_ratio = largerM_body->getMass() / smallerM_body->getMass();
    
//     // Vectors to hold fragment mass data
//     std::vector<double> fragmentsA_masses;
//     std::vector<double> fragmentsB_masses;
    
//     // --- Case 1: Both Bodies Break (Masses are similar) ---
//     if (mass_ratio <= FRAGMENTATION_RATIO) {
//         // Break both bodies
//         int numFragsA = MIN_FRAGMENTS + (rand() % (MAX_FRAGMENTS - MIN_FRAGMENTS + 1));
//         int numFragsB = MIN_FRAGMENTS + (rand() % (MAX_FRAGMENTS - MIN_FRAGMENTS + 1));
        
//         fragmentsA_masses = generateRandomFragmentMasses(mA, numFragsA);
//         fragmentsB_masses = generateRandomFragmentMasses(mB, numFragsB);
        
//         // Mark both original bodies for deletion
//         bodyA.setMarkedForDeletion(true);
//         bodyB.setMarkedForDeletion(true);
//     }
//     // --- Case 2: Only Smaller Body Breaks (Mass ratio is high) ---
//     else {
//         // Only the smaller body fragments, the larger body absorbs the impact (and survives)
//         double smallerMass = smallerM_body->getMass();
//         int numFrags = MIN_FRAGMENTS + (rand() % (MAX_FRAGMENTS - MIN_FRAGMENTS + 1));
        
//         // Determine which vector to populate based on which body is smaller
//         if (smallerM_body == &bodyA) {
//              fragmentsA_masses = generateRandomFragmentMasses(smallerMass, numFrags);
//         } else {
//              fragmentsB_masses = generateRandomFragmentMasses(smallerMass, numFrags);
//         }
        
//         // Mark only the smaller original body for deletion
//         smallerM_body->setMarkedForDeletion(true);
        
//         // The larger body (largerM_body) survives, but we must update its velocity 
//         // to reflect the momentum it absorbed from the smaller body.
//         // It now inherits the V_CoM of the entire system (the simplest approach).
//         largerM_body->setNetVelocity(V_CoM); 
//     }

//     // Helper: produce fragments for one body and return them
//     auto makeFragments = [&](GravitationalBody& sourceBody, GravitationalBody& nonSourceBody, std::vector<double>& masses)
//         {
//             std::vector<GravitationalBody> frags;
//             if (masses.empty()) return frags;

//             // Vector2D collisionCenter = (bodyA.getPosition() * mA + bodyB.getPosition() * mB) / M;
//             Vector2D direction = (sourceBody.getPosition() - nonSourceBody.getPosition()).normalize();
//             Vector2D scaledDirection = direction * nonSourceBody.getRadius();
//             Vector2D collisionPoint = nonSourceBody.getPosition() + scaledDirection;
//             Vector2D sourcePos = sourceBody.getPosition();
//             double sourceRadius = sourceBody.getRadius();

//             // Tuneable: The distance fragments scatter *from the collision center*
//             const double INITIAL_SCATTER_DISTANCE = sourceRadius * 0.75; 
            
//             for (double fMass : masses)
//             {
//                 GravitationalBody frag;
//                 frag.setMass(fMass);

//                 // Radius scales with M^(1/3)
//                 frag.setRadius(sourceRadius * std::pow(fMass / sourceBody.getMass(), 1.0/3.0));
                
//                 // Initial position: Scatter fragments from the collision center in a random direction
//                 // The position is CRUCIAL for stability; they must not overlap immediately.
//                 Vector2D scatterDir = randomDirectionVector();
//                 frag.setPosition(collisionPoint + scatterDir * INITIAL_SCATTER_DISTANCE);
//                 frag.setPrevPosition(frag.getPosition());

//                 // Base Velocity: All fragments start with the system's CoM velocity
//                 frag.setNetVelocity(V_CoM); 

//                 // Set flags
//                 frag.setIsFragment(true);

//                 frag.setCollisionEnabled(false);
                
//                 // Bounding box update is needed right away
//                 // We'll update the bounding box and then later apply the final scatter velocity
//                 // to ensure the BBox is large enough.
//                 // updateBoundingBox(frag); // Let's skip BBox update until after the final scatter velocity is applied

//                 frags.push_back(frag);
//             }
//             return frags;
//     };
//     // Instantiate all fragments
//     std::vector<GravitationalBody> allFragments;
    
//     if (!fragmentsA_masses.empty()) {
//         auto fragsA = makeFragments(bodyA, bodyB, fragmentsA_masses);
//         allFragments.insert(allFragments.end(), fragsA.begin(), fragsA.end());
//     }
//     if (!fragmentsB_masses.empty()) {
//         auto fragsB = makeFragments(bodyB, bodyA, fragmentsB_masses);
//         allFragments.insert(allFragments.end(), fragsB.begin(), fragsB.end());
//     }
//     // --- 4. Apply Impact-Driven Scatter Velocity ---
    
//     // TUNEABLE: Controls the intensity of the scatter (e.g., 0.1 to 0.5)
//     // Higher K means fragments fly apart faster after impact.
//     const double FRAGMENTATION_COEFFICIENT = 0.3; 
    
//     // Scatter magnitude is proportional to the relative impact speed
//     double scatterMagnitude = relativeSpeed * FRAGMENTATION_COEFFICIENT;
    
//     Vector2D collisionCenter = (bodyA.getPosition() * mA + bodyB.getPosition() * mB) / M;

//     for (auto& frag : allFragments)
//     {
//         // 1. Get the direction from the collision center (C) to the fragment's position (P_frag)
//         Vector2D direction_to_fragment = frag.getPosition() - collisionCenter;
        
//         // 2. Normalize the direction vector
//         Vector2D scatterDir = direction_to_fragment.normalize();
        
//         // 3. Calculate the scatter velocity vector
//         Vector2D scatterVel = scatterDir * scatterMagnitude;
        
//         // 4. Add the scatter velocity to the fragment's base velocity (V_CoM)
//         frag.setNetVelocity(frag.getNetVelocity() + scatterVel);
        
//         // NOTE: The total system momentum is now slightly incorrect due to randomization 
//         // in initial fragment positions and the fragmentation process. This is fixed next.
//     }
//     // --- 5. Final Momentum Correction and Stability Setup ---
    
//     // 5a. Calculate Raw Momentum of Fragments
//     Vector2D P_raw(0.0, 0.0);
//     double M_total_fragments = 0.0;
    
//     for (const auto& frag : allFragments)
//     {
//         P_raw += frag.getNetVelocity() * frag.getMass();
//         M_total_fragments += frag.getMass();
//     }
    
//     // 5b. Calculate Correction
//     // P_initial is the target momentum (from Step 1)
//     Vector2D P_diff = P_initial - P_raw;
    
//     // Apply correction to the entire mass of the fragments
//     Vector2D V_corr = P_diff / M_total_fragments;
    
//     // // 5c. Apply Correction and Final Setup
//     // const int GHOST_FRAMES = 5; // TUNEABLE: How many frames to disable collisions
    
//     for (auto& frag : allFragments)
//     {
//         // Apply the correction to ensure P_final == P_initial
//         frag.setNetVelocity(frag.getNetVelocity() + V_corr);
        
//     //     // Enable Ghosting/Cooldown to prevent immediate re-collision cascade
//     //     frag.setGhost(GHOST_FRAMES);
        
//     //     // Update the bounding box now that the final velocity is known
//     //     updateBoundingBox(frag);
//     }

//     // 5d. Insert New Fragments into Game State
//     // Since bodyA/bodyB were marked for deletion in Step 2, these are the replacements.
//     for (auto& f : allFragments)
//     {
//         state.addBody(f);
//     }
// } // End of handleDynamicCollision

void PhysicsSystem::handleElasticCollision(GravitationalBody& A, GravitationalBody& B)
{
    Vector2D x = B.getPosition() - A.getPosition();
    double dist = x.magnitude();

    if (dist == 0) return;
    Vector2D n = x / dist; // collision normal

    Vector2D vA = A.getNetVelocity();
    Vector2D vB = B.getNetVelocity();

    double mA = A.getMass();
    double mB = B.getMass();

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
    A.setNetVelocity(vA + vA_change);
    B.setNetVelocity(vB + vB_change);

    // --- Positional correction to prevent re-collision ---
    double overlap = A.getRadius() + B.getRadius() - dist;
    if (overlap > 0) {
        A.setPosition(A.getPosition() - n * (overlap * 0.5));
        B.setPosition(B.getPosition() + n * (overlap * 0.5));
    }
}


// /// MINE
// void PhysicsSystem::handleDynamicCollision(GravitationalBody& bodyA, GravitationalBody& bodyB, GameState& state)
// {
//     // if (!bodyA.getCollisionEnabled() || !bodyB.getCollisionEnabled()){
//     //     return;
//     // }
        
//     double mA = bodyA.getMass();
//     double mB = bodyB.getMass();
//     Vector2D vA = bodyA.getNetVelocity();
//     Vector2D vB = bodyB.getNetVelocity();

//     Vector2D P_initial = vA * mA + vB * mB; // total momentum

//     // Helper: produce fragments for one body and return them
//     auto makeFragmentsForBody = [&](GravitationalBody& source,
//                                     std::vector<double>& masses)
//     {
//         std::vector<GravitationalBody> frags;
//         frags.reserve(masses.size());

//         Vector2D pos = source.getPosition();

//         for (double fMass : masses)
//         {
//             GravitationalBody frag;
//             frag.setMass(fMass);

//             frag.setRadius(source.getRadius() * sqrt(fMass / source.getMass()));
//             frag.setPosition(pos + randomDirectionVector() * source.getRadius() * 0.3);
//             frag.setPrevPosition(frag.getPosition());
//             Vector2D scatter = randomDirectionVector() * (1.0 + rand() % 20);
//             frag.setNetVelocity(scatter);
//             frag.setIsFragment(true);
//             // instantiate as no collision (to prevent cascading failures)
//             frag.setCollisionEnabled(false);
            
//             updateBoundingBox(frag);


//             frags.push_back(frag);
//         }

//         return frags;
//     };
//     // Helper: momentum-correct fragment list
//     auto fixMomentum = [&](std::vector<GravitationalBody>& frags, Vector2D targetP)
//     {
//         double mTotal = 0.0;
//         Vector2D P_raw(0,0);

//         for (auto& f : frags) {
//             mTotal += f.getMass();
//             P_raw += f.getNetVelocity() * f.getMass();
//         }

//         Vector2D dP = targetP - P_raw;
//         Vector2D offsetVel = dP / mTotal;

//         for (auto& f : frags) {
//             f.setNetVelocity(f.getNetVelocity() + offsetVel);
//         }
//     };

//     // Vector to hold all fragments to add
//     std::vector<GravitationalBody> allFragments;

//     int nA = 8 + (rand( )% 16);
//     int nB = 8 + (rand() % 16);
//     // ------------------------------
//     // CASE 1: Both have equal mass → both break
//     // ------------------------------
//     if (mA == mB)
//     {

//         auto massesA = generateRandomFragmentMasses(mA, nA);
//         auto massesB = generateRandomFragmentMasses(mB, nB);

//         auto fragsA = makeFragmentsForBody(bodyA, massesA);
//         auto fragsB = makeFragmentsForBody(bodyB, massesB);

//         // Momentum-fix fragments individually if desired
//         fixMomentum(fragsA, vA * mA);
//         fixMomentum(fragsB, vB * mB);

//         allFragments.insert(allFragments.end(), fragsA.begin(), fragsA.end());
//         allFragments.insert(allFragments.end(), fragsB.begin(), fragsB.end());

//         bodyA.setMarkedForDeletion(true);
//         bodyB.setMarkedForDeletion(true);
//     }
//     // ------------------------------
//     // CASE 2: A heavier → B breaks
//     // ------------------------------
//     else if (mA > mB)
//     {
//         Vector2D vA_final = P_initial / (mA + mB);
//         bodyA.setNetVelocity(vA_final);

//         Vector2D P_frag_target = P_initial - vA_final * mA;

//         auto massesB = generateRandomFragmentMasses(mB, nB);
//         auto fragsB = makeFragmentsForBody(bodyB, massesB);

//         fixMomentum(fragsB, P_frag_target);

//         allFragments.insert(allFragments.end(), fragsB.begin(), fragsB.end());
//         bodyB.setMarkedForDeletion(true);
//     }
//     // ------------------------------
//     // CASE 3: B heavier → A breaks
//     // ------------------------------
//     else
//     {
//         Vector2D vB_final = P_initial / (mA + mB);
//         bodyB.setNetVelocity(vB_final);

//         Vector2D P_frag_target = P_initial - vB_final * mB;

//         auto massesA = generateRandomFragmentMasses(mA, nA);
//         auto fragsA = makeFragmentsForBody(bodyA, massesA);

//         fixMomentum(fragsA, P_frag_target);

//         allFragments.insert(allFragments.end(), fragsA.begin(), fragsA.end());
//         bodyA.setMarkedForDeletion(true);
//     }
//        // ----- Resolve initial overlaps among fragments with elastic separation -----
//     for (size_t i = 0; i < allFragments.size(); ++i) {
//         for (size_t j = i + 1; j < allFragments.size(); ++j) {
//             Vector2D delta = allFragments[j].getPosition() - allFragments[i].getPosition();
//             double dist = delta.magnitude();
//             double minDist = allFragments[i].getRadius() + allFragments[j].getRadius();
//             if (dist < minDist && dist > 0) {
//                 Vector2D dir = delta.normalize();
//                 double penetration = minDist - dist;

//                 // move positions slightly apart
//                 // allFragments[i].setPosition(allFragments[i].getPosition() - dir * penetration * 1.0);
//                 // allFragments[j].setPosition(allFragments[j].getPosition() + dir * penetration * 1.0);

//                 double separationFraction = 0.5; // 50% of penetration
//                 allFragments[i].setPosition(allFragments[i].getPosition() - dir * penetration * separationFraction);
//                 allFragments[j].setPosition(allFragments[j].getPosition() + dir * penetration * separationFraction);


//                 // apply elastic velocity push along separation
//                 // allFragments[i].setNetVelocity(allFragments[i].getNetVelocity() - dir * penetration / PHYSICS_TIME_STEP * 0.01);
//                 // allFragments[j].setNetVelocity(allFragments[j].getNetVelocity() + dir * penetration / PHYSICS_TIME_STEP * 0.01);
//                 double velocityPushFactor = 0.20; // tweakable
//                 Vector2D push = dir * (penetration / PHYSICS_TIME_STEP * velocityPushFactor * separationFraction);
//                 allFragments[i].setNetVelocity(allFragments[i].getNetVelocity() - push);
//                 allFragments[j].setNetVelocity(allFragments[j].getNetVelocity() + push);

                


//             }
//         }
//     }
//     // //     // ----- Enable collisions immediately -----
//     for (auto& f : allFragments){
//         // f.setCollisionEnabled(true);
//         f.setCollisionEnabled(false);
//         // f.setGhost(3);
//     }
//     // ---- ADD ALL FRAGMENTS AT ONCE ----
//     for (auto& f : allFragments)
//         state.addBody(f);
// }

// void PhysicsSystem::handleAccretion(GravitationalBody& smallerM, GravitationalBody& largerM)
// {
//     // 1. Conservation of Momentum
//     double m1 = smallerM.getMass();
//     double m2 = largerM.getMass();
//     Vector2D v1 = smallerM.getNetVelocity();
//     Vector2D v2 = largerM.getNetVelocity();

//     double newMass = m1 + m2;
//     // P_final = P1_initial + P2_initial
//     Vector2D newVelocity = (v1 * m1 + v2 * m2) / newMass;

//     // 2. Update Collector Body
//     largerM.setMass(newMass);
//     largerM.setNetVelocity(newVelocity);
    
//     // Update radius: Scale by mass ratio, assuming constant density (R ∝ M^(1/3))
//     double newRadius = largerM.getRadius() * std::pow(newMass / m1, 1.0/3.0);
//     largerM.setRadius(newRadius);
//     updateBoundingBox(largerM);

//     // 3. Mark the Absorbed Fragment for Deletion
//     smallerM.setMarkedForDeletion(true);
// }
std::vector<double> PhysicsSystem::generateRandomFragmentMasses(double totalMass, int numFragments)
{
    std::vector<double> masses;
    masses.resize(numFragments);

    std::vector<double> weights(numFragments);

    double W = 0.0;
    for (int i = 0; i < numFragments; i++)
    {
        // Wide variation
        // double w = -log((double)rand()/RAND_MAX); 

        // Mild variation
        // double x = randNormal();  // via Box-Muller
        // double w = fabs(x);

        // Explosion of tiny objects
        double w = pow((double)rand()/RAND_MAX, 3.0); // higher exponent = more tiny fragments

        // double w = (double)rand() / RAND_MAX; 
        if (w < 1e-6) w = 1e-6;  // avoid zero mass fragments
        weights[i] = w;
        W += w;
    }

    for (int i = 0; i < numFragments; i++)
        masses[i] = totalMass * (weights[i] / W);

    return masses;
}

Vector2D PhysicsSystem::randomDirectionVector()
{
    double angle = ((double)rand() / RAND_MAX) * 2.0 * PI;
    return Vector2D(cos(angle), sin(angle));
}
// Return true if `body` overlaps ANY other body in `bodies` (ignoring itself and deleted bodies)
// bool PhysicsSystem::isOverlappingAny(const GravitationalBody& body, const GameState& state) {
//     for (const auto& other : state.getBodies()) {
//         if (&other == &body) continue; // skip self
//         if (other.getMarkedForDeletion()) continue;
//         if (other.getIsGhost()) continue;
//         // Quick AABB reject/accept first
//         // if (!body.getBoundingBox().overlaps(other.getBoundingBox())) continue;
//         // precise radius check
//         Vector2D d = other.getPosition() - body.getPosition();
//         double dist = d.magnitude();
//         double minDist = body.getRadius() + other.getRadius();
//         if (dist <= minDist) return true;
//     }
//     return false;
// }
void PhysicsSystem::updateBoundingBox(GravitationalBody& body) {
    // next pos propagated but with a fudge factor due to potential extreme acceleration in between frames. Probably need to do a derivation for what the maximum is to most optimize the size of the box
    // but this is proof of concept.
    Vector2D expectedNextPos = body.getPosition() + body.getNetVelocity()*PHYSICS_TIME_STEP;
    Vector2D currentPos = body.getPosition();
    float radius = body.getRadius();

    double minX = std::min(currentPos.x_val, expectedNextPos.x_val) - radius;
    double maxX = std::max(currentPos.x_val, expectedNextPos.x_val) + radius;
    double minY = std::min(currentPos.y_val, expectedNextPos.y_val) - radius;
    double maxY = std::max(currentPos.y_val, expectedNextPos.y_val) + radius;

    // Expand by 50% for safety margin
    const double expansionFactor = 0.05;
    double expandX = (maxX - minX) * expansionFactor * 0.5;
    double expandY = (maxY - minY) * expansionFactor * 0.5;

    BoundingBox box;
    box.min_X = minX - expandX;
    box.max_X = maxX + expandX;
    box.min_Y = minY - expandY;
    box.max_Y = maxY + expandY;

    body.setBoundingBox(box);
}