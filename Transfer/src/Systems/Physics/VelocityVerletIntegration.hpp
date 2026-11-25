// File: Transfer/src/Systems/Physics/VelocityVerletIntegration.hpp

#pragma once

// Custom Imports
#include "Systems/Physics/UniverseState.hpp"

// Velocity Verlet ("kick, drift, kick"). One physics tick is:
//   1. integrateVelocityVerletFirstHalfStep   half kick with last tick's forces, then move
//   2. the forces are recomputed at the new positions (PhysicsSystem::recomputeAllForces)
//   3. integrateVelocityVerletSecondHalfStep  the other half kick with the new forces
// The order matters: step 1 zeroes every netForce, so step 2 must refill it before step 3 reads it.
// PhysicsSystem::integrateUniversePhysics runs these three; the real tick and the trajectory prediction both use it.

// One body, first half: half the velocity change from last tick's net force, then a full position step.
// Saves previousPosition (for render interpolation) and prevForce, then zeroes netForce for this tick's forces.
// Massless and static bodies skip the kick but still move with their velocity.
void applyVelocityVerletFirstHalfStep(GravitationalBody& grav_body);
// One body, second half: the other half of the velocity change, from the net force at the new position.
// Massless and static bodies are left alone.
void applyVelocityVerletSecondHalfStep(GravitationalBody& grav_body);

// The first half step for every particle, macro body and the ship
void integrateVelocityVerletFirstHalfStep(UniverseState& universe_state);
// The second half step for every particle, macro body and the ship (after the forces were recomputed)
void integrateVelocityVerletSecondHalfStep(UniverseState& universe_state);