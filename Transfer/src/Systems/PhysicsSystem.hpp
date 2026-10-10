// File: Transfer/src/Systems/PhysicsSystem.hpp

#pragma once

// Custom Imports
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "DynamoEngine/Math/Vector2.hpp"
#include "Entities/Physics/GravitationalBody.hpp"
#include "Entities/Physics/GravitationalBodyPair.hpp"
#include "Utilities/Constants/EngineConstants.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"
#include "Utilities/Constants/PhysicsConstants.hpp"
#include "Utilities/Math/CustomMathUtilities.hpp"
#include "Utilities/Physics/UniformParticleGrid.hpp"
#include "Utilities/Rendering/CameraTransform.hpp"

// Standard Library Imports
#include <algorithm>
#include <cmath>

#include <iostream>

struct CollisionInfo
{
    double distance;
    DynamoEngine::Vector2D unitNormalVector;       // unit normal (from bodyA to bodyB)
    DynamoEngine::Vector2D relativeVelocityVector; // vB - vA
    double normalSpeed;                            // signed speed along normal vector
    double absNormalSpeed;                         // abs value of signed speed along normal vector
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
    void UpdateSystemFrame(GameState& game_state, UIState& uiState);
    // Helper method called in the destructor to clear up physics-related
    // contents
    void CleanUp();
    void UpdateGravBodyInstantiations(GameState& game_state, UIState& uiState);

  private:
    // --- Collision Handling ---
    // Top-level collision handler. Makes decisions about the kinds of
    // collisions encountered and dispatches to the subhandlers
    void handleCollisions(GameState& game_state);
    void handleMacroMacroCollisions(GameState& game_state);
    void handleMacroParticleCollisions(GameState& game_state);
    void handleParticleParticleCollisions(GameState& game_state);
    void handleDynamicCollision(GravitationalBodyPair& gravBodyPair, const CollisionInfo& collisionInfo,
                                GameState& game_state);
    // Handles a 'bouncy' (elastic) collision between two bodies, when the collision
    // satisfies Engine-Constant-defined constraints
    void handleElasticCollisions(GravitationalBody& smallerBody, GravitationalBody& largerBody);
    void handleAccretion(GravitationalBodyPair& gravBodyPair);
    void promoteOversizedParticles(GameState& game_state); // TODO: Prune? currently uncalled, see UpdateSystemFrame
    // Both append the new fragments to `fragments_out`. During collisions that is m_pending_fragments, never the
    // particles vector itself: the collision loops are still walking over (and holding references into) particles.
    void substituteWithParticles(GravitationalBody& original_body, std::vector<GravitationalBody>& fragments_out,
                                 uint32_t targetFragmentCount);
    void substituteWithParticlesFromImpact(GravitationalBody& original_body,
                                           std::vector<GravitationalBody>& fragments_out, uint32_t targetFragmentCount,
                                           const DynamoEngine::Vector2D& impactPoint);
    // Particles alive now plus fragments waiting to join them (what MAX_LIVE_PARTICLES limits)
    size_t liveParticleCount(const GameState& game_state) const;
    // Live particles plus the worst case still to come: every macro body shattering into DEFAULT_FRAGMENT_COUNT
    size_t potentialParticleCount(const GameState& game_state) const;
    // Checks if there is room for a given number of new particles within the MAX_LIVE_PARTICLES limit
    bool hasRoomForParticles(const GameState& game_state, size_t new_particle_count) const;

    // --- Gravity ---
    void updateAllForces(GameState& game_state); // Gravity calculation dispatch helper
    void updateGravityForSystem(GameState& game_state);
    void calculateGravity(GravitationalBody& firstBody,
                          GravitationalBody& secondBody); // Calculate and apply gravity between two
                                                          // gravitational bodies

    // --- Integration (Velocity Verlet) ---
    void integrateForwardsVelocityVerletPhase1(GameState& game_state);
    void applyVelocityVerletPhase1(GravitationalBody& gravBody);
    void integrateForwardsVelocityVerletPhase2(GameState& game_state);
    void applyVelocityVerletPhase2(GravitationalBody& gravBody);

    // --- Gravitational Body Creation Mechanisms ---
    void createMacroBody(GameState& game_state,
                         DEPRECATED_InputState& inputState); // Creates a Macro Gravitational Body
                                                             // with the user-defined attributes
    void createParticle(GameState& game_state,
                        DEPRECATED_InputState& inputState); // TODO: Prune? declared, never defined or called
    void createParticleCluster(GameState& game_state,
                               DEPRECATED_InputState& inputState); // TODO: Prune? declared, never defined or called

    // --- Utility ---
    void calculateTotalEnergy(GameState& game_state); // TODO: Prune? currently uncalled, see UpdateSystemFrame
                                                      //  Calculates total energy of all Macro Bodies and
                                                      //  Particles on screen.

    // --- Player Physics ---
    void updatePlayerPhysics(GameState& game_state, UIState& uiState);
    // --- Cleanup ---
    void cleanupParticles(GameState& game_state);   // Clears any Particles from the screen flagged
                                                    // as marked for deletion
    void cleanupMacroBodies(GameState& game_state); // Clears any Macro Bodies from the screen
                                                    // flagged as marked for deletion

    uint32_t survivableFragmentCount(const GravitationalBody& body, uint32_t maxCount);
    // --- Data Members ---
    UniformParticleGrid particleGrid;

    // Fragments created during this tick's collision passes. Appended to the particles vector only after all three
    // passes are done, so nothing is ever added to particles while a loop is iterating over it.
    std::vector<GravitationalBody> m_pending_fragments;
};