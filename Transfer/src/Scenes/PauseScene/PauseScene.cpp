// File: Transfer/src/Scenes/PauseScene/PauseScene.cpp

#include "Scenes/PauseScene/PauseScene.h"

PauseScene::PauseScene() : Scene(SceneIdentifier::PAUSE_SCENE)
{
    sceneUIElements.insert({UIElementIdentifier::PLAY_GAME_BUTTON_INDEX, nullptr});
}

void PauseScene::populateMe()
{
    // std::cout << "gamescenepopulatecalled" << std::endl;
    PlayGameButton* play_game_button = new PlayGameButton();
    for (auto& [UI_element_ID, UI_element_ptr] : sceneUIElements)
    {
        switch (UI_element_ID)
        {
        case UIElementIdentifier::PLAY_GAME_BUTTON_INDEX:
            UI_element_ptr = play_game_button;
            break;
        default:
            break;
        }
    }
}