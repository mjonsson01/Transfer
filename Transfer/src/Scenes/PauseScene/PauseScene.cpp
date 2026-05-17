// File: Transfer/src/Scenes/PauseScene/PauseScene.cpp

#include "Scenes/PauseScene/PauseScene.h"

PauseScene::PauseScene() : Scene(SceneIdentifier::PAUSE_SCENE) {}

void PauseScene::populateMe()
{
    // PlayGameButton* play_game_button = new PlayGameButton();
    for (auto& [UI_element_ID, UI_element_ptr] : sceneUIElements)
    {
        switch (UI_element_ID)
        {
        default:
            break;
        }
    }
}