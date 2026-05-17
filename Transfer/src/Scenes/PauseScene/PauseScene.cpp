// File: Transfer/src/Scenes/PauseScene/PauseScene.cpp

#include "Scenes/PauseScene/PauseScene.h"

PauseScene::PauseScene() : Scene(SceneIdentifier::PAUSE_SCENE)
{
    sceneUIElements.insert({UIElementIdentifier::RESUME_BUTTON_INDEX, nullptr});
}

void PauseScene::populateMe()
{
    // PlayGameButton* play_game_button = new PlayGameButton();
    for (auto& [UI_element_ID, UI_element_ptr] : sceneUIElements)
    {
        ResumeButton* resume_button = new ResumeButton();
        switch (UI_element_ID)
        {
        case UIElementIdentifier::RESUME_BUTTON_INDEX:
            UI_element_ptr = resume_button;
            break;
        default:
            break;
        }
    }
}