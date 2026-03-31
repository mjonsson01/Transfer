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
        float getFPS() { return framesPerSecond;}
        void setFPS(float fps) { framesPerSecond = fps;}
        bool getAllUIVisibility() {return allUIElementsVisible;}
        void invertUIElementsVisibility() { allUIElementsVisible = !allUIElementsVisible; }

        bool getLevelScene() {return levelScene;}
        void setLevelScene(bool ls) {levelScene = ls;}
        bool getPauseMenuActive() {return pauseMenuActive;}
        void setPauseMenuActive(bool pma) {pauseMenuActive = pma;}
        bool getStartMenuActive() {return startMenuActive;}
        void setStartMenuActive(bool sma) {startMenuActive = sma;}
        bool getLevelEditorScene() {return levelEditorScene;}
        void setLevelEditorScene(bool les) {levelEditorScene = les;}
        bool getRenderDebug() {return renderDebug;}
        void setRenderDebug(bool rd) {renderDebug = rd;}
    private:
        InputState inputState;
        float framesPerSecond = TARGET_FPS;
        bool startMenuActive = false;
        bool pauseMenuActive = false;
        bool levelScene = false;
        bool levelEditorScene = false;
        bool renderDebug = false; // Toggles rendering of debug elements like collision boxes, spawn areas, etc.
        bool allUIElementsVisible = true; // Default to true because we want all elements visible.
};
