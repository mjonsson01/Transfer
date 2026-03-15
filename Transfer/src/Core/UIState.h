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
        bool getAllVisibility() {return allUIElementsVisible;}
        // bool getFPSVisibility() { return allUIElementsVisibility.FPSVisibile;}

    private:
        InputState inputState;
        float framesPerSecond = TARGET_FPS;
        bool allUIElementsVisible = true; // Default to true because we want all elements visible.
};
