// File: Transfer/src/Scenes/TransferScenes.hpp

#pragma once

// Custom Imports
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "DynamoEngine/Scenes/Scene.hpp"
#include "DynamoEngine/Scenes/SceneManager.hpp"

// Every scene in Transfer. A plain enum (not an enum class) inside a namespace: the names are written
// TransferScene::Game, and they are already SceneIds, so no casts are needed.
namespace TransferScene
{
enum : DynamoEngine::SceneId
{
    StartMenu,
    Game,
    Pause,
    TestVisual,
};
} // namespace TransferScene

// Builds every Transfer scene and registers it with `scenes`.
// The UI writes the player's choices (slider values, sounds) into `ui_state`, so `ui_state` must outlive the scenes.
// The game must outlive the scenes for the same reason, the lambdas are full of references (primarily to the camera
// state, in game state)
void addTransferScenes(DynamoEngine::SceneManager& scenes, UIState& ui_state, GameState& game_state);
