// File: Transfer/src/Systems/Input/InputSystem.hpp

#pragma once

// Custom Imports
#include "Camera/CameraLimits.hpp"
#include "Camera/CameraState.hpp"
#include "Core/GameState.hpp"
#include "Systems/Input/PlayerInput.hpp"
#include "Core/UIState.hpp"
#include "Scenes/TransferScenes.hpp"
#include "Constants/EngineConstants.hpp"
#include "Constants/GameSystemConstants.hpp"
#include "Utilities/CustomMathUtilities.hpp"
#include "Camera/CameraTransform.hpp"

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
// then turns what's left into Transfer's meaning (camera controls, scene changes, and PlayerInput: ship controls,
// spawn aiming and one-frame commands).
class InputSystem
{
  public:
    // Constructor and Destructor
    InputSystem();
    ~InputSystem();

  public:
    // Main method to process input. `frame_seconds` is how long the last frame took (for UI updates).
    // Rewrites `player_input` for this frame (its one-frame commands are cleared first).
    void processSystemInputFrame(GameState& game_state, UIState& ui_state, DynamoEngine::SceneManager& scenes,
                                 PlayerInput& player_input, float frame_seconds);

    // Clean up helper
    void cleanUp();

  private:
    // Game rule: where a creation drag started (shift re-anchors it). Runs BEFORE the event is applied,
    // because it needs the button state from before this event.
    void trackDragAnchor(const DynamoEngine::InputEvent& event);
    // Lays out and updates the scene's UI for the current window, then gives it this frame's mouse input
    DynamoEngine::UIInputResult updateSceneUI(DynamoEngine::UIRoot& ui, const CameraState& camera_state,
                                              float frame_seconds);
    // Zoom around the cursor, pan (middle-mouse drag or the arrow keys), then keep the camera inside the world
    void updateCamera(GameState& game_state, float frame_seconds);
    // `pointer_captured`: the UI took the mouse this frame, so a button release there doesn't spawn anything
    void translateGameInputs(GameState& game_state, UIState& ui_state, DynamoEngine::SceneManager& scenes,
                             PlayerInput& player_input, bool pointer_captured);
    // Esc: the pause menu returns to the paused scene, the start menu starts the game
    void translateMenuInputs(const UIState& ui_state, DynamoEngine::SceneManager& scenes);
    // Where the next body would spawn right now, in world space: at the cursor, or (with Shift) at the drag's start,
    // moving toward the cursor. The one place this is worked out: the spawn and both previews all use it.
    SpawnProperties aimSpawn(bool is_cluster, bool with_velocity, const CameraState& camera_state) const;

    void updateVisor(CameraState& camera_state); // tab cycles the visor view;

  private:
    DynamoEngine::SDLInputIntake m_intake; // direct SDL events translated through the intake
    DynamoEngine::InputState m_input;
    std::vector<DynamoEngine::InputEvent> m_frame_events; // Events piled on a frame by frame basis
    DynamoEngine::Vector2F m_mouse_drag_anchor;           // Transfer's drag start (pressPosition + shift re-anchor)
};