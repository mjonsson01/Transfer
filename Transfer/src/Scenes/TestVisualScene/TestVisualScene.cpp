// File: Transfer/src/Scenes/TestVisualScene/TestVisualScene.cpp

#include "Scenes/TestVisualScene/TestVisualScene.hpp"

TestVisualScene::TestVisualScene() : Scene(SceneIdentifier::PAUSE_SCENE)
{
    // sceneUIElements.insert({UIElementIdentifier::RESUME_BUTTON_INDEX, nullptr});
    sceneUIElements.insert({UIElementIdentifier::VISOR_MENU_SELECTION_INDEX, nullptr});
    // No ui elements for now
}

void TestVisualScene::populateMe()
{
    DropDownMenu* visor_menu = new DropDownMenu();
    for (auto& [UI_element_ID, UI_element_ptr] : sceneUIElements)
    {
        switch (UI_element_ID)
        {
        case UIElementIdentifier::VISOR_MENU_SELECTION_INDEX:
            UI_element_ptr = visor_menu;
            break;
        default:
            break;
        }
    }
}