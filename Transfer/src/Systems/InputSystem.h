// File: Transfer/src/Systems/InputSystem.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>

// Custom Imports
#include "Core/GameState.h"
#include "Core/UIState.h"
#include "Utilities/EngineConstants.h"
#include "Utilities/GameSystemConstants.h"
#include "Utilities/TransferInputs.h"

// Standard Library Imports
#include <iostream>

// struct LevelEditorInputs. TBI

class InputSystem {
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
    void routeSDL_EventInputInGame(SDL_Event* event);
    void translateAndPassTransferInputsOff(UIState& UIState);
    TransferInputs transferInputs; // in game inputs
};