// File: Transfer/src/Core/GameState.h

#pragma once

// Custom Imports
#include "Entities/Physics/GravitationalBody.h"
#include "Entities/Sound/MusicModeEnum.h"
#include "Utilities/Constants/GameSystemConstants.h"

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

    // Getters for Particles with a mutable and nonmutable version to respect
    // System hierarchies
    const std::vector<GravitationalBody>& getParticles() const { return particles; }
    std::vector<GravitationalBody>& getParticlesMutable() { return particles; }

    // Getters for Bodies with a mutable and nonmutable version to respect
    // System hierarchies
    const std::vector<GravitationalBody>& getMacroBodies() const { return macroBodies; }
    std::vector<GravitationalBody>& getMacroBodiesMutable() { return macroBodies; }

    // Getter and setter for the alpha rendering variable
    float getAlpha() const { return alpha; }
    void setAlpha(float alphaIn) { alpha = alphaIn; }

    void incrementMaxIDInstantiated() { maxIDInstantiated += 1; }
    int getMaxIDInstantiated() const { return maxIDInstantiated; }

  private:
    // State variables
    bool isPlaying = false;
    bool isShuttingDownAudioSystem = false;
    // Frame helper vars
    float alpha = 0.0f;

    // Database for all the Macro Bodies
    std::vector<GravitationalBody> macroBodies;
    int maxIDInstantiated = 0;
    // Database for all the Particle Bodies
    std::vector<GravitationalBody> particles;
};
