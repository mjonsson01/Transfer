// File: Transfer/src/Systems/Physics/VelocityVerletIntegration.cpp

#include "Systems/Physics/VelocityVerletIntegration.hpp"
#include "Utilities/Math/CustomMathUtilities.hpp"

void integrateVelocityVerletFirstHalfStep(UniverseState& universe_state)
{
    std::vector<GravitationalBody>& particles = universe_state.particles;
    std::vector<GravitationalBody>& macro_bodies = universe_state.macro_bodies;

    for (auto& particle : particles)
    {
        applyVelocityVerletFirstHalfStep(particle);
    }
    for (auto& macro_body : macro_bodies)
    {
        applyVelocityVerletFirstHalfStep(macro_body);
    }
    universe_state.ship.applyVelocityVerletFirstHalfStep();
}

void applyVelocityVerletFirstHalfStep(GravitationalBody& grav_body)
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

void applyVelocityVerletSecondHalfStep(GravitationalBody& grav_body)
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

void integrateVelocityVerletSecondHalfStep(UniverseState& universe_state)
{
    std::vector<GravitationalBody>& particles = universe_state.particles;
    std::vector<GravitationalBody>& macro_bodies = universe_state.macro_bodies;

    for (auto& particle : particles)
    {
        applyVelocityVerletSecondHalfStep(particle);
    }
    for (auto& macro_body : macro_bodies)
    {
        applyVelocityVerletSecondHalfStep(macro_body);
    }
    universe_state.ship.applyVelocityVerletSecondHalfStep();
}
