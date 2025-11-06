#pragma once

#include "Core/GameState.h"
#include "SDL3/SDL.h"
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
        void ProcessSystemInputFrame(GameState& state);

    private:
        // Internal state variables for input handling can be added here
    
    private: 
        // subordinate rendering functions.
        void createNewBody(SDL_Event& event, GameState& state);

};