// File: Transfer/src/Systems/InputSystem.hpp

#pragma once

// Custom Imports
#include "Core/CameraState.hpp"
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "Scenes/TransferScenes.hpp"
#include "Utilities/Constants/EngineConstants.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"
#include "Utilities/Math/CustomMathUtilities.hpp"
#include "Utilities/Rendering/CameraTransform.hpp"

// Engine Imports
#include "DynamoEngine/Input/InputEvent.hpp"
#include "DynamoEngine/Input/InputState.hpp"
#include "DynamoEngine/Input/SDLInputIntake.hpp"
#include "DynamoEngine/Math/Vector2.hpp"
#include "DynamoEngine/Scenes/SceneManager.hpp"
#include "DynamoEngine/UI/UIInputResult.hpp"
#include "DynamoEngine/UI/UIRoot.hpp"

// Standard Library Imports
#include <algorithm>
#include <cmath>
#include <iostream>

// Game-Side Input: reads engine's InputState, lets the current scene's UI take its share first,
// then turns what's left into Transfer's meaning (camera controls, spawn requests, scene changes),
// written into DEPRECATED_InputState for now.
class InputSystem
{
  public:
    // Constructor and Destructor
    InputSystem();
    ~InputSystem();

  public:
    // Main method to process input. `frame_seconds` is how long the last frame took (for UI updates).
    void processSystemInputFrame(GameState& game_state, UIState& ui_state, DynamoEngine::SceneManager& scenes,
                                 float frame_seconds);

    // Clean up helper
    void cleanUp();

  private:
    // Game rule: where a creation drag started (shift re-anchors it). Runs BEFORE the event is applied,
    // because it needs the button state from before this event.
    void trackDragAnchor(const DynamoEngine::InputEvent& event);
    // Lays out and updates the scene's UI for the current window, then gives it this frame's mouse input
    DynamoEngine::UIInputResult updateSceneUI(DynamoEngine::UIRoot& ui, const CameraState& camera_state,
                                              float frame_seconds);
    void updateCamera(GameState& game_state); // zoom around cursor, middle-mouse pan, star-field clamp
    void translateGameInputs(GameState& game_state, UIState& ui_state, DynamoEngine::SceneManager& scenes);
    void translateMenuInputs(UIState& ui_state, DynamoEngine::SceneManager& scenes);
    void copySharedPointerState(DEPRECATED_InputState& legacy_state); // fields both translators pass on

    void updateVisor(CameraState& camera_state); // tab cycles the visor view;

  private:
    DynamoEngine::SDLInputIntake m_intake; // direct SDL events translated through the intake
    DynamoEngine::InputState m_input;
    std::vector<DynamoEngine::InputEvent> m_frame_events; // Events piled on a frame by frame basis
    DynamoEngine::Vector2F m_mouse_drag_anchor;           // Transfer's drag start (pressPosition + shift re-anchor)
};