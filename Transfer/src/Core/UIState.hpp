// File: Transfer/src/Core/UIState.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>

// Custom Imports
#include "Core/DEPRECATED_InputState.hpp"
#include "Entities/Sound/MusicModeEnum.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"

// Standard Library Imports
#include <queue>
#include <string>
#include <vector>

class UIState
{
  public:
    UIState();
    ~UIState();
    DEPRECATED_InputState& getMutableDEPRECATED_InputState() { return inputState; }
    const DEPRECATED_InputState& getDEPRECATED_InputState() const { return inputState; }
    float getFPS() { return framesPerSecond; }
    void setFPS(float fps) { framesPerSecond = fps; }
    bool getRenderDebug() { return renderDebug; }
    void setRenderDebug(bool rd) { renderDebug = rd; }
    float getTimeScaleFactor() const { return static_cast<float>(inputState.selectedSimSpeedScale); }

    void QueueSoundEffect(const std::string& soundName) { pendingSoundEffects.push(soundName); };
    bool HasPendingSoundEffects() const { return !pendingSoundEffects.empty(); };
    std::string PopNextSoundEffect()
    {
        if (pendingSoundEffects.empty())
            return "";

        std::string sound = pendingSoundEffects.front();
        pendingSoundEffects.pop();

        return sound;
    }
    MusicMode getRequestedMusicMode() const { return requestedMusicMode; }
    void setRequestedMusicMode(MusicMode musicMode) { requestedMusicMode = musicMode; }
    bool getPlayMusic() const { return playMusic; }
    void setPlayMusic(bool spm) { playMusic = spm; }
    bool getPlaySoundEffects() const { return playSoundEffects; }
    void setPlaySoundEffects(bool pse) { playSoundEffects = pse; }

  private:
    DEPRECATED_InputState inputState;
    float framesPerSecond = TARGET_FPS;

    bool renderDebug = VIEW_DEBUG; // Toggles rendering of debug elements like
    // collision boxes, spawn areas, etc.

    // SOUND STUFF
    bool playMusic = false;
    bool playSoundEffects = false;
    MusicMode requestedMusicMode = MusicMode::NONE;
    std::queue<std::string> pendingSoundEffects;
};
