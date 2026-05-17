// File: Transfer/src/Core/UIState.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>

// Custom Imports
#include "Core/InputState.h"
#include "Entities/UIElements/UIElementIdentifierEnum.h"
#include "Scenes/SceneIdentifierEnum.h"
#include "Utilities/Constants/GameSystemConstants.h"

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

  private:
    InputState inputState;
    float framesPerSecond = TARGET_FPS;

    bool renderDebug = VIEW_DEBUG;    // Toggles rendering of debug elements like
                                      // collision boxes, spawn areas, etc.
    bool allUIElementsVisible = true; // Default to true because we want all elements visible.

    SceneIdentifier currentScene = SceneIdentifier::NO_SCENE;
};
