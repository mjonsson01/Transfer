// File: Transfer/src/Systems/InputSystem.h

#pragma once

// SDL3 Imports
#include "SDL3/SDL.h"

// Custom Imports
#include "Core/UIState.h"
#include "Core/GameState.h"
#include "Utilities/GameSystemConstants.h"
#include "Utilities/EngineConstants.h"

// Standard Library Imports
#include <iostream>
class InputSystem
{
    public:
        // Constructor and Destructor
        InputSystem();
        ~InputSystem();

    public:
        // Main method to process input
        void ProcessSystemInputFrame(GameState& gameState, UIState& UIState);

        // Clean up helper
        void CleanUp();

    private: 
        // void handleMassSliderInput(SDL_Event& event, UIState& UIState); 
};