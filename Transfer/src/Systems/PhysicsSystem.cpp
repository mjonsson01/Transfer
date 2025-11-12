#include "PhysicsSystem.h"
#include "iostream"
PhysicsSystem::PhysicsSystem()
{
    // Initialize physics system variables if needed
}

PhysicsSystem::~PhysicsSystem()
{
    // Cleanup if necessary
}

void PhysicsSystem::UpdateSystemFrame(GameState& state)
{
    // Update physics calculations for all entities in the game state

    // Extract body vector as modifiable
    std::vector<GravitationalBody>& bodies = state.getBodiesMutable();


	for (size_t i = 0; i < bodies.size(); ++i) {
		for (size_t j = i + 1; j < bodies.size(); ++j) {
            calculateGravForceBetweenBodies(bodies[i], bodies[j]);
		}
	}


    // integrate the system forward after all the resultant forces are calculated 
    for (size_t i = 0; i < bodies.size(); ++i){
        updateBodyVectors(bodies[i]);
    }

    // handle collisions and update velocities of anything detected as colliding. later
    handleCollisions(state);
}

void PhysicsSystem::updateBodyVectors(GravitationalBody& body )
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

    // Extract values of interest for legibility
    double massA = bodyA.getMass();
    double massB = bodyB.getMass();

    // Find distance between bodies
    Vector2D posA = bodyA.getPosition();
    Vector2D posB = bodyB.getPosition();
    Vector2D deltaDisp = posB - posA;
    double distance = deltaDisp.magnitude();

    if (distance > bodyA.getRadius() + bodyB.getRadius()){ 
        // sin(theta),cos(theta) unit vector
        Vector2D trigUnitVectors = deltaDisp.normalize();
        double force_magnitude = GRAVITATIONAL_CONSTANT * massA * massB / deltaDisp.square_magnitude();
        Vector2D force_B_on_A = Vector2D(force_magnitude * trigUnitVectors.x_val, force_magnitude * trigUnitVectors.y_val);
        // Newton's Third Law
        Vector2D force_A_on_B = force_B_on_A * (-1.0);
        bodyA.setNetForce(force_B_on_A);
        bodyB.setNetForce(force_A_on_B);
    }
    else{
        // collisions handled elsewhere.
        return;
    }
}

void PhysicsSystem::handleCollisions(GameState& state)
{
    int size = state.getBodies().size();
    std::vector<GravitationalBody>& bodies = state.getBodiesMutable();
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = i + 1; j < size; ++j) {

            GravitationalBody& bodyA = bodies[i];
            GravitationalBody& bodyB = bodies[j];
            // Check bounding boxes to reduce search space. Bounding box method limited if extreme speed is in play.
            if (!bodyA.getBoundingBox().overlaps(bodyB.getBoundingBox()))
                continue; // skip

            // Actually deal with the collisions
            Vector2D direction_vector = bodyB.getPosition() - bodyA.getPosition();
            double dist = direction_vector.magnitude();
            double min_distance = bodyA.getRadius() + bodyB.getRadius();
            if (dist < min_distance && dist > 0.0) {
                Vector2D unit_dir_vector = direction_vector.normalize();
                handleElasticCollision(bodyA, bodyB, unit_dir_vector); // if collision is hard enough, break the mass :D TBI
            }
        }
    }
}

void PhysicsSystem::handleElasticCollision(GravitationalBody& bodyA, GravitationalBody& bodyB, Vector2D& unit_dir_vector)
{
    // unit dir vector points from A to B
    // They will need to transfer momentum to one another and resolve based on their masses. 
    // conservation of momentum says m1v1i+m2v2i = m1v1f+m2v2f
    // And conservation of energy (assuming elastic) says 1/2 m1 v1i^2 + 1/2 m2 v2i^2 = 1/2m1 vf^2 + 1/2 m2 vf^2 
    double m_a = bodyA.getMass();
    double m_b = bodyB.getMass();
    Vector2D v_a_initial = bodyA.getNetVelocity();
    Vector2D v_b_initial = bodyB.getNetVelocity();
    Vector2D bodyA_final_vel = v_a_initial * (m_a - m_b)/(m_a + m_b) + v_b_initial * (2*m_b/(m_a+m_b));
    Vector2D bodyB_final_vel = v_b_initial * (m_a - m_b)/(m_a + m_b) + v_a_initial * (2*m_a/(m_a+m_b));
    bodyA.setNetVelocity(bodyA_final_vel);
    bodyB.setNetVelocity(bodyB_final_vel);  
}

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
    const double expansionFactor = 0.5;
    double expandX = (maxX - minX) * expansionFactor * 0.5;
    double expandY = (maxY - minY) * expansionFactor * 0.5;

    BoundingBox box;
    box.min_X = minX - expandX;
    box.max_X = maxX + expandX;
    box.min_Y = minY - expandY;
    box.max_Y = maxY + expandY;

    body.setBoundingBox(box);
}