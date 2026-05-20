// File: Transfer/src/Scenes/StartMenuScene/StartMenuScene.cpp

#include "StartMenuScene.hpp"

StartMenuScene::StartMenuScene() : Scene(SceneIdentifier::GAME_SCENE)
{
    sceneUIElements.insert({UIElementIdentifier::PLAY_GAME_BUTTON_INDEX, nullptr});
}

void StartMenuScene::populateMe()
{
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