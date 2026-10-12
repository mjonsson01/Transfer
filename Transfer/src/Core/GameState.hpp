// File: Transfer/src/Core/GameState.h

#pragma once

// Custom Imports
#include "Core/CameraState.hpp"
#include "Core/UniverseState.hpp"
#include "Entities/Sound/MusicModeEnum.hpp"
#include "Player/Player.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"
// Standard Library Imports
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

    UniverseState& getUniverseStateMutable() { return m_universe_state; }
    const UniverseState& getUniverseState() const { return m_universe_state; }

  private:
    // State variables
    bool isPlaying = false;
    bool isShuttingDownAudioSystem = false;
    // Frame helper vars
    float alpha = 0.0f;

    // Contains all Physics integratable objects
    UniverseState m_universe_state;

    Player player;

    CameraState cameraState;
};
