// File: Transfer/src/Core/UIState.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>

// Custom Imports
#include "Core/InputState.h"
#include "Utilities/GameSystemConstants.h"

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

    bool getGameScene() { return levelScene; }
    void setGameScene(bool ls) { levelScene = ls; }
    bool getLevelEditorScene() { return levelEditorScene; }
    void setLevelEditorScene(bool les) { levelEditorScene = les; }
    bool getLevelSelectScene() { return levelSelectScene; }
    void setLevelSelectScene(bool lss) { levelEditorScene = lss; }
    bool getPauseMenuActive() { return pauseMenuActive; }
    void setPauseMenuActive(bool pma) { pauseMenuActive = pma; }
    bool getStartMenuActive() { return startMenuActive; }
    void setStartMenuActive(bool sma) { startMenuActive = sma; }
    bool getRenderDebug() { return renderDebug; }
    void setRenderDebug(bool rd) { renderDebug = rd; }

  private:
    InputState inputState;
    float framesPerSecond = TARGET_FPS;

    // Menus
    bool startMenuActive = false;
    bool pauseMenuActive = false;

    // Scenes
    bool levelScene = false;
    bool levelEditorScene = false;
    bool levelSelectScene = false;

    bool renderDebug = false;         // Toggles rendering of debug elements like
                                      // collision boxes, spawn areas, etc.
    bool allUIElementsVisible = true; // Default to true because we want all elements visible.
};
