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
        // Helper method to handle Keyboard Inputs
        void handleKeyboardInput(SDL_Event& event, GameState& gameState, UIState& UIState);
        // Helper method to handle Mouse Inputs
        void handleMouseInput(SDL_Event& event, GameState& gameState, UIState& UIState);



        // void handleMassSliderInput(SDL_Event& event, UIState& UIState); 
};