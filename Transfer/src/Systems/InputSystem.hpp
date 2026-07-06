// File: Transfer/src/Systems/InputSystem.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>

// Custom Imports
#include "Core/CameraState.hpp"
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "Scenes/SceneIdentifierEnum.hpp"
#include "Utilities/Constants/EngineConstants.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"
#include "Utilities/Math/CustomMathUtilities.hpp"
#include "Utilities/Rendering/CameraTransform.hpp"
#include "Utilities/UserInput/TransferInputs.hpp"
// Standard Library Imports
#include <algorithm>
#include <cmath>
#include <iostream>
// struct LevelEditorInputs. TBI

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
    void routeSDL_EventInputInGame(SDL_Event* event);
    void routeSDL_EventInputInMenu(SDL_Event* event);
    void translateAndPassTransferInputsOff(UIState& UIState);
    void translateAndPassMenuInputsOff(UIState& UIState);
    TransferInputs transferInputs; // in game inputs
};