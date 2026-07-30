// File: Transfer/src/Systems/PhysicsSystem.h

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
#include <numeric>
#include <random>

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

  public:
    // Method to update Physics System. Handles all physics interactions and
    // body instantiations for one physics frame
    void UpdateSystemFrame(GameState& gameState, UIState& UIState);
    // Helper method called in the destructor to clear up physics-related
    // contents
    void CleanUp();
    void UpdateGravBodyInstantiations(GameState& gameState, UIState& UIState);

  private:
    void updateAllForces(GameState& gameState);
    void handleMacroMacroCollisions(GameState& gameState);
    void handleMacroParticleCollisions(GameState& gameState);
    void handleParticleParticleCollisions(GameState& gameState);
    void handleDynamicCollision(GravitationalBodyPair& gravBodyPair, const CollisionInfo& collisionInfo,
                                GameState& gameState);
    void handleAccretion(GravitationalBodyPair& gravBodyPair);
    void promoteOversizedParticles(GameState& gameState);
    void substituteWithParticles(GravitationalBody& originalBody, GameState& gameState, uint32_t targetFragmentCount);
    void substituteWithParticlesFromImpact(GravitationalBody& originalBody, GameState& gameState,
                                           uint32_t targetFragmantCount, const Vector2D& impactPoint);

    void handleCollisions(GameState& gameState);

    // Sub-level collision handlers
    void handleElasticCollisions(GravitationalBody& smallerBody,
                                 GravitationalBody& largerBody); // Handles 'bouncy' collision if the collision
                                                                 // satisfies Engine-Constant-defined constraints
                                                                 // the pieces if the collision satisfies

    UniformParticleGrid particleGrid;

    // Gravity Methods
    void updateGravityForSystem(GameState& gameState); // Gravity calculation dispatch helper
    void calculateGravity(GravitationalBody& body1,
                          GravitationalBody& body2); // Calculate and apply gravity between two
                                                     // gravitational bodies

    // Integration system
    void integrateForwardsVelocityVerletPhase1(GameState& gameState);
    void applyVelocityVerletPhase1(GravitationalBody& gravBody);
    void integrateForwardsVelocityVerletPhase2(GameState& gameState);
    void applyVelocityVerletPhase2(GravitationalBody& gravBody);
    // Top-level collision handler. Makes decisions about the kinds of
    // collisions encountered and dispatches to the subhandlers

    // Gravitational Body Creation Mechanisms
    void createMacroBody(GameState& gameState,
                         InputState& inputState); // Creates a Macro Gravitational Body
                                                  // with the user-defined attributes
    void createParticle(GameState& gameState,
                        InputState& inputState); // Creates a Particle Gravitational Body with
                                                 // the user-defined attriubutes
    void createParticleCluster(GameState& gameState,
                               InputState& inputState); // Creates a cluster of Particle Gravitational
                                                        // Bodies with the user-defined attributes

    // Utility Functions
    void calculateTotalEnergy(GameState& gameState); // Calculates total energy of all Macro Bodies
                                                     // and Particles on Screen.
    // Cleanup Functions
    void cleanupParticles(GameState& gameState);   // Clears any Particles from the screen flagged
                                                   // as marked for deletion
    void cleanupMacroBodies(GameState& gameState); // Clears any Macro Bodies from the screen
                                                   // flagged as marked for deletion
};