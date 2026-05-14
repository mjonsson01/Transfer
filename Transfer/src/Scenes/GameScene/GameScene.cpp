// File: Transfer/src/Scenes/GameScene/GameScene.cpp

#include "GameScene.h"

GameScene::GameScene() : Scene(SceneIdentifier::GAME_SCENE)
{
    sceneUIElements.insert({UIElementIdentifier::FPS_COUNTER_INDEX, nullptr});
    sceneUIElements.insert({UIElementIdentifier::MASS_SLIDER_INDEX, nullptr});
    sceneUIElements.insert({UIElementIdentifier::RADIUS_SLIDER_INDEX, nullptr});
    sceneUIElements.insert({UIElementIdentifier::SIMULATION_SPEED_SLIDER_INDEX, nullptr});
}

void GameScene::populateMe()
{
    // std::cout << "gamescenepopulatecalled" << std::endl;
    FPSCounter* fps_counter = new FPSCounter();
    MassSlider* mass_slider = new MassSlider();
    SimulationSpeedSlider* simulation_speed_slider = new SimulationSpeedSlider();
    RadiusSlider* radius_slider = new RadiusSlider();
    for (auto& [UI_element_ID, UI_element_ptr] : sceneUIElements)
    {
        switch (UI_element_ID)
        {
        case UIElementIdentifier::FPS_COUNTER_INDEX:
            UI_element_ptr = fps_counter;
            break;
        case UIElementIdentifier::MASS_SLIDER_INDEX:
            UI_element_ptr = mass_slider;
            break;
        case UIElementIdentifier::RADIUS_SLIDER_INDEX:
            UI_element_ptr = radius_slider;
            break;
        case UIElementIdentifier::SIMULATION_SPEED_SLIDER_INDEX:
            UI_element_ptr = simulation_speed_slider;
            break;
        default:
            break;
        }
    }
}