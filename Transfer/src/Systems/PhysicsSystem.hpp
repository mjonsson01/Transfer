// File: Transfer/src/Systems/PhysicsSystem.hpp

#pragma once

// Custom Imports
#include "Core/GameState.hpp"
#include "Core/PlayerInput.hpp"
#include "Core/SpawnSettings.hpp"
#include "Core/Trajectories.hpp"
#include "Core/UIState.hpp"
#include "Core/UniverseState.hpp"
#include "DynamoEngine/Math/Vector2.hpp"
#include "DynamoEngine/Physics/Collision2D.hpp"
#include "Entities/Physics/GravitationalBody.hpp"
#include "Entities/Physics/GravitationalBodyPair.hpp"
#include "Systems/Physics/VelocityVerletIntegration.hpp"
#include "Utilities/Constants/EngineConstants.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"
#include "Utilities/Constants/PhysicsConstants.hpp"
#include "Utilities/Math/CustomMathUtilities.hpp"
#include "Utilities/Physics/FragmentLayout.hpp"
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
    void advanceUniverseState(UniverseState& universe_state, const ShipControls& ship_controls);
    // Helper method called in the destructor to clear up physics-related
    // contents
    void cleanUp();
    // This frame's spawn, undo and clear requests (PlayerInput::commands), with the spawn panel's settings
    void updateGravBodyInstantiations(UniverseState& universe_state, const PlayerInput& player_input,
                                      const SpawnSettings& spawn_settings);
    // Planning mode stops time, but W/S still set up the burn: moves the ship's throttle in REAL time (frame_seconds)
    void adjustThrottleWhilePlanning(UniverseState& universe_state, const ShipControls& ship_controls,
                                     float frame_seconds);
    // Once per FRAME (not per tick). Planning mode: handles the I request, wipes the lines when the universe was
    // edited, and grows them. Outside planning mode: clears them. Defined in Systems/Physics/TrajectoryPrediction.cpp
    // `player_input` gives the I request and the spawn being aimed (for the placement preview)
    void updateTrajectories(float frame_seconds, const UniverseState& universe_state, const PlayerInput& player_input,
                            Trajectories& trajectories, const UIState& ui_state);

  private:
    void integrateUniversePhysics(UniverseState& universe_state);
    // Removes the most recently spawned planet or cluster that still exists, with all of its debris (plain Delete)
    void removeNewestSpawn(UniverseState& universe_state);

    // --- Collision Handling ---
    // Top-level collision handler. Makes decisions about the kinds of
    // collisions encountered and dispatches to the subhandlers
    void handleCollisionEvents(UniverseState& universe_state);
    void handleMacroMacroCollisions(UniverseState& universe_state);
    void handleMacroParticleCollisions(UniverseState& universe_state);
    void handleParticleParticleCollisions(UniverseState& universe_state);
    // The player's ship against planets and debris: never overlapping, bouncing with SHIP_RESTITUTION
    void handleShipCollisions(UniverseState& universe_state);
    // Pushes the ship and one body apart along the contact normal, then bounces them if they're still closing.
    // Returns the kinetic energy the bounce absorbed (0 if they weren't closing).
    double resolveShipContact(Starship& ship, GravitationalBody& body, const DynamoEngine::CircleContact& contact);

    void handleDynamicCollision(GravitationalBodyPair& grav_body_pair, const CollisionInfo& collisionInfo,
                                UniverseState& universe_state);
    // Handles a 'bouncy' (elastic) collision between two bodies, when the collision
    // satisfies Engine-Constant-defined constraints
    void handleElasticCollisions(GravitationalBody& smallerBody, GravitationalBody& largerBody);
    void handleAccretion(GravitationalBodyPair& grav_body_pair);
    // TODO: Prune? currently uncalled, see updateSystemFrame
    void promoteOversizedParticles(UniverseState& universe_state);
    // Both append the new fragments to `fragments_out`. During collisions that is m_pending_fragments, never the
    // particles vector itself: the collision loops are still walking over (and holding references into) particles.
    void substituteWithParticles(GravitationalBody& original_body, std::vector<GravitationalBody>& fragments_out,
                                 uint32_t target_fragment_count);
    void substituteWithParticlesFromImpact(GravitationalBody& original_body,
                                           std::vector<GravitationalBody>& fragments_out,
                                           uint32_t target_fragment_count, const DynamoEngine::Vector2D& impactPoint);
    // Particles alive now plus fragments waiting to join them (what MAX_LIVE_PARTICLES limits)
    size_t liveParticleCount(const UniverseState& universe_state) const;
    // Live particles plus the worst case still to come: every macro body shattering into DEFAULT_FRAGMENT_COUNT
    size_t potentialParticleCount(const UniverseState& universe_state) const;
    // Checks if there is room for a given number of new particles within the MAX_LIVE_PARTICLES limit
    bool hasRoomForParticles(const UniverseState& universe_state, size_t new_particle_count) const;

    // --- Gravity ---
    void recomputeAllForces(UniverseState& universe_state); // Gravity calculation dispatch helper
    void addGravityBetweenGravBodies(UniverseState& universe_state);
    void addGravityBetweenGravBodyPair(GravitationalBody& firstBody,
                                       GravitationalBody& secondBody); // Calculate and apply gravity between two
                                                                       // gravitational bodies
    // Gravity between the player's ship and every planet, both ways (the ship pulls the gravitational bodies too)
    void addGravityBetweenShipAndBodies(UniverseState& universe_state);

    // --- Gravitational Body Creation Mechanisms ---
    void createMacroBody(UniverseState& universe_state, const SpawnProperties& spawn_aim,
                         const SpawnSettings& spawn_settings); // Creates a Macro Gravitational Body
                                                               // with the user-defined attributes
    void createParticleCluster(UniverseState& universe_state, const SpawnProperties& spawn_aim,
                               const SpawnSettings& spawn_settings);

    // --- Player Physics ---
    void applyShipInputPhysics(UniverseState& universe_state, const ShipControls& ship_controls);
    // --- Cleanup ---
    void cleanupGravBodies(UniverseState& universe_state);  // calls the subcleaners in succession
    void cleanupParticles(UniverseState& universe_state);   // Clears any Particles from the screen flagged
                                                            // as marked for deletion
    void cleanupMacroBodies(UniverseState& universe_state); // Clears any Macro Bodies from the screen
                                                            // flagged as marked for deletion
    // --- Trajectory prediction (planning mode), defined in Systems/Physics/TrajectoryPrediction.cpp ---
    // Copies the universe into the prediction and starts the paths (one per macro body, plus the ship)
    void startTrajectories(const UniverseState& universe_state, Trajectories& trajectories);
    // Drops the paths and the prediction: back to Idle
    void clearTrajectories(Trajectories& trajectories);
    // True when the real universe no longer matches the one the trace started from (a spawn, a delete, ...)
    bool universeChangedSinceTraceStart(const UniverseState& universe_state) const;
    // The integrator only (no collisions): steps the prediction for this frame's share of ticks, recording points
    void advanceIntegratedTrajectories(float frame_seconds, Trajectories& trajectories);
    // Appends the prediction's current positions to every path
    void recordTrajectoryPoints(Trajectories& trajectories);
    // While a planet is being placed (left button held on the world): predicts where it would go. Restarts whenever
    // the planet changes (the mouse moved), then grows as fast as MAX_PREVIEW_MILLISECONDS_PER_FRAME allows.
    // The preview universe has NO debris (it would make every restart slow): it ignores particle clouds, the traced
    // lines don't.
    void updateSpawnPreviewTrajectory(const UniverseState& universe_state, const PlayerInput& player_input,
                                      const SpawnSettings& spawn_settings, Trajectories& trajectories);
    void clearSpawnPreviewTrajectory(Trajectories& trajectories);
    // The planet a left-button release would spawn right now (defined in PhysicsSystem.cpp next to createMacroBody,
    // built by the same function, so the preview can never disagree with the real spawn). False if none would spawn.
    bool makeSpawnPreviewBody(const SpawnProperties& spawn_aim, const SpawnSettings& spawn_settings,
                              GravitationalBody& body_out) const;

    // --- Data Members ---
    UniformParticleGrid particleGrid;

    // Fragments created during this tick's collision passes. Appended to the particles vector only after all three
    // passes are done, so nothing is ever added to particles while a loop is iterating over it.
    std::vector<GravitationalBody> m_pending_fragments;

    // Trajectory prediction
    UniverseState m_prediction_universe_state;  // the copy being stepped forward
    UniverseState m_trace_start_universe_state; // the real universe when the trace started, to notice edits
    float m_trace_time_accumulator = 0.0f;      // game seconds the prediction still owes (same idea as Game::Run's)
    int m_trace_ticks_done = 0;                 // tracing stops at the horizon

    // Placement preview (its own little prediction, restarted whenever the planet being placed changes)
    UniverseState m_preview_universe_state; // the planets + ship + the planet being placed (always the last body)
    GravitationalBody m_preview_spawn_body; // the planet the preview was started for, to notice when it changes
    int m_preview_ticks_done = 0;
    bool m_has_spawn_preview = false;
};