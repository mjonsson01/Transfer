// File: Transfer/src/Systems/InputSystem.h

#pragma once

// SDL3 Imports
#include "SDL3/SDL.h"

// Custom Imports
#include "Core/UIState.h"
#include "Core/GameState.h"
#include "Utilities/GameSystemConstants.h"
#include "Utilities/EngineConstants.h"

class InputSystem
{
    public:
        // Constructor and Destructor
        InputSystem();
        ~InputSystem();

    public:
        // Methods to process input
        void ProcessSystemInputFrame(GameState& state, UIState& UIState);

        void CleanUp();

    private:
        // Internal state variables for input handling can be added here
    
    private: 
        // subordinate rendering functions.
        void createNewBody(SDL_Event& event, GameState& state);

};