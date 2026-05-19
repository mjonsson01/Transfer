// File: Transfer/src/Core/UIState.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>

// Custom Imports
#include "Core/InputState.hpp"
#include "Entities/Sound/MusicModeEnum.hpp"
#include "Entities/UIElements/UIElementIdentifierEnum.hpp"
#include "Scenes/SceneIdentifierEnum.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"

// Standard Library Imports
#include <vector>

class UIState
{
  public:
    UIState();
    ~UIState();
    InputState& getMutableInputState() { return inputState; }
    const InputState& getInputState() const { return inputState; }
    float getFPS() { return framesPerSecond; }
    void setFPS(float fps) { framesPerSecond = fps; }
    bool getAllUIVisibility() { return allUIElementsVisible; }
    void invertUIElementsVisibility() { allUIElementsVisible = !allUIElementsVisible; }
    bool getRenderDebug() { return renderDebug; }
    void setRenderDebug(bool rd) { renderDebug = rd; }
    float getTimeScaleFactor() const { return static_cast<float>(inputState.selectedSimSpeedScale); }
    SceneIdentifier getCurrentSceneID() const { return currentScene; }
    void setCurrentScene(SceneIdentifier scene_desired) { currentScene = scene_desired; }
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
    InputState inputState;
    float framesPerSecond = TARGET_FPS;

    bool renderDebug = VIEW_DEBUG; // Toggles rendering of debug elements like
    // collision boxes, spawn areas, etc.
    bool allUIElementsVisible = true; // Default to true because we want all elements visible.

    SceneIdentifier currentScene = SceneIdentifier::NONE;
    // SOUND STUFF
    bool playMusic = false;
    bool playSoundEffects = false;
    MusicMode requestedMusicMode = MusicMode::NONE;
    std::queue<std::string> pendingSoundEffects;
};
