// File: Transfer/src/Core/GameState.h

#pragma once

// Custom Imports
#include "Camera/CameraState.hpp"
#include "Systems/Physics/Trajectories.hpp"
#include "Systems/Physics/UniverseState.hpp"
#include "Systems/Audio/MusicModeEnum.hpp"
#include "Player/Player.hpp"
#include "Constants/GameSystemConstants.hpp"
// Standard Library Imports
#include <functional>
#include <utility>
#include <vector>

class GameState
{
  public:
    // Constructor and Destructor
    GameState();
    ~GameState();

  public:
    // Getters and Setters for Game State
    bool IsPlaying() const { return isPlaying; }
    void SetPlaying(bool playing) { isPlaying = playing; }

    bool getIsShuttingDownAudioSystem() const { return isShuttingDownAudioSystem; }
    void setIsShuttingDownAudioSystem(bool shutdown) { isShuttingDownAudioSystem = shutdown; }

    // Getter and setter for the alpha rendering variable
    float getAlpha() const { return alpha; }
    void setAlpha(float alpha_in) { alpha = alpha_in; }

    CameraState& getCameraStateMutable() { return cameraState; }
    const CameraState& getCameraState() const { return cameraState; }

    Player& getPlayerMutable() { return player; }
    const Player& getPlayer() const { return player; }

    // In the running game this is the main thread's COPY of the simulation (Core/PhysicsThread.hpp): what is drawn,
    // read by the UI and the trajectory lines, and replaced by a newer copy every frame. Change the simulation with
    // editUniverse, never through getUniverseStateMutable (the next copy would undo it). Tools that step a universe
    // themselves, without the physics thread (the determinism check, harnesses), use it directly.
    UniverseState& getUniverseStateMutable() { return m_universe_state; }
    const UniverseState& getUniverseState() const { return m_universe_state; }

    // Queues a change to the simulation's universe (e.g. a scene adding its collision box). Game hands the queued
    // changes to the physics thread every frame, which applies them between two ticks; the copy here shows them as
    // soon as the next copy arrives.
    void editUniverse(std::function<void(UniverseState&)> edit) { m_pending_universe_edits.push_back(std::move(edit)); }
    // The queued changes (Game takes them every frame)
    std::vector<std::function<void(UniverseState&)>>& pendingUniverseEdits() { return m_pending_universe_edits; }

    Trajectories& getTrajectoriesMutable() { return m_trajectories; }
    const Trajectories& getTrajectories() const { return m_trajectories; }

  private:
    // State variables
    bool isPlaying = false;
    bool isShuttingDownAudioSystem = false;
    // Frame helper vars
    float alpha = 0.0f;

    // Contains all Physics integratable objects
    UniverseState m_universe_state;
    std::vector<std::function<void(UniverseState&)>> m_pending_universe_edits;
    Trajectories m_trajectories;

    Player player;

    CameraState cameraState;
};
