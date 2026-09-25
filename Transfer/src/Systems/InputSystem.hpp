// File: Transfer/src/Systems/InputSystem.hpp

#pragma once

// Custom Imports
#include "Core/CameraState.hpp"
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "Scenes/SceneIdentifierEnum.hpp"
#include "Utilities/Constants/EngineConstants.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"
#include "Utilities/Math/CustomMathUtilities.hpp"
#include "Utilities/Rendering/CameraTransform.hpp"

// Engine Imports
#include "DynamoEngine/Input/InputEvent.hpp"
#include "DynamoEngine/Input/InputState.hpp"
#include "DynamoEngine/Input/SDLInputIntake.hpp"
#include "DynamoEngine/Math/Vector2.hpp"

// Standard Library Imports
#include <algorithm>
#include <cmath>
#include <iostream>

// Game-Side Input: reads engine's InputState and turns it into Transfer's meaning
// (camera controls, spawn requests, scene changes), written into DEPRECATED_InputState for now.
class InputSystem
{
  public:
    // Constructor and Destructor
    InputSystem();
    ~InputSystem();

  public:
    // Main method to process input
    void processSystemInputFrame(GameState& game_state, UIState& ui_state);

    // Clean up helper
    void cleanUp();

  private:
    // Game rule: where a creation drag started (shift re-anchors it). Runs BEFORE the event is applied,
    // because it needs the button state from before this event.
    void trackDragAnchor(const DynamoEngine::InputEvent& event);
    void updateCamera(GameState& game_state); // zoom around cursor, middle-mouse pan, star-field clamp
    void translateGameInputs(UIState& ui_state);
    void translateMenuInputs(UIState& ui_state);
    void copySharedPointerState(DEPRECATED_InputState& legacy_state); // fields both translators pass on

  private:
    DynamoEngine::SDLInputIntake m_intake; // direct SDL events translated through the intake
    DynamoEngine::InputState m_input;
    std::vector<DynamoEngine::InputEvent> m_frame_events; // Events piled on a frame by frame basis
    DynamoEngine::Vector2F m_mouse_drag_anchor;           // Transfer's drag start (pressPosition + shift re-anchor)
};