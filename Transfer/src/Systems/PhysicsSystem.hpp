// File: Transfer/src/Systems/PhysicsSystem.hpp

#pragma once

// Custom Imports
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "Entities/Physics/GravitationalBody.hpp"
#include "Entities/Physics/GravitationalBodyPair.hpp"
#include "Utilities/Constants/EngineConstants.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"
#include "Utilities/Math/CustomMathUtilities.hpp"
#include "Utilities/Math/Vector2D.hpp"
#include "Utilities/Physics/UniformParticleGrid.hpp"
#include "Utilities/Rendering/CameraTransform.hpp"

// Standard Library Imports
#include <algorithm>
#include <cmath>

#include <iostream>

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

class PhysicsSystem
{
  public:
    // Constructor and Destructor
    PhysicsSystem();
    ~PhysicsSystem();

    // Method to update Physics System. Handles all physics interactions and
    // body instantiations for one physics frame
    void UpdateSystemFrame(GameState& gameState, UIState& uiState);
    // Helper method called in the destructor to clear up physics-related
    // contents
    void CleanUp();
    void UpdateGravBodyInstantiations(GameState& gameState, UIState& uiState);

  private:
    // --- Collision Handling ---
    // Top-level collision handler. Makes decisions about the kinds of
    // collisions encountered and dispatches to the subhandlers
    void handleCollisions(GameState& gameState);
    void handleMacroMacroCollisions(GameState& gameState);
    void handleMacroParticleCollisions(GameState& gameState);
    void handleParticleParticleCollisions(GameState& gameState);
    void handleDynamicCollision(GravitationalBodyPair& gravBodyPair, const CollisionInfo& collisionInfo,
                                GameState& gameState);
    // Handles a 'bouncy' (elastic) collision between two bodies, when the collision
    // satisfies Engine-Constant-defined constraints
    void handleElasticCollisions(GravitationalBody& smallerBody, GravitationalBody& largerBody);
    void handleAccretion(GravitationalBodyPair& gravBodyPair);
    void promoteOversizedParticles(GameState& gameState); // TODO: Prune? currently uncalled, see UpdateSystemFrame
    void substituteWithParticles(GravitationalBody& originalBody, GameState& gameState, uint32_t targetFragmentCount);
    void substituteWithParticlesFromImpact(GravitationalBody& originalBody, GameState& gameState,
                                           uint32_t targetFragmentCount, const Vector2D& impactPoint);

    // --- Gravity ---
    void updateAllForces(GameState& gameState); // Gravity calculation dispatch helper
    void updateGravityForSystem(GameState& gameState);
    void calculateGravity(GravitationalBody& firstBody,
                          GravitationalBody& secondBody); // Calculate and apply gravity between two
                                                          // gravitational bodies

    // --- Integration (Velocity Verlet) ---
    void integrateForwardsVelocityVerletPhase1(GameState& gameState);
    void applyVelocityVerletPhase1(GravitationalBody& gravBody);
    void integrateForwardsVelocityVerletPhase2(GameState& gameState);
    void applyVelocityVerletPhase2(GravitationalBody& gravBody);

    // --- Gravitational Body Creation Mechanisms ---
    void createMacroBody(GameState& gameState,
                         InputState& inputState); // Creates a Macro Gravitational Body
                                                  // with the user-defined attributes
    void createParticle(GameState& gameState,
                        InputState& inputState); // TODO: Prune? declared, never defined or called
    void createParticleCluster(GameState& gameState,
                               InputState& inputState); // TODO: Prune? declared, never defined or called

    // --- Utility ---
    void calculateTotalEnergy(GameState& gameState); // TODO: Prune? currently uncalled, see UpdateSystemFrame
                                                     //  Calculates total energy of all Macro Bodies and
                                                     //  Particles on screen.

    // --- Cleanup ---
    void cleanupParticles(GameState& gameState);   // Clears any Particles from the screen flagged
                                                   // as marked for deletion
    void cleanupMacroBodies(GameState& gameState); // Clears any Macro Bodies from the screen
                                                   // flagged as marked for deletion

    // --- Data Members ---
    UniformParticleGrid particleGrid;
};