// File: Transfer/src/Scenes/GameScene/GameScene.h

#pragma once

#include "Entities/UIElements/Buttons/PlayGameButton/PlayGameButton.hpp"
#include "Entities/UIElements/Overlay/FPSCounter.hpp"
#include "Entities/UIElements/Sliders/MassSlider.hpp"
#include "Entities/UIElements/Sliders/RadiusSlider.hpp"
#include "Entities/UIElements/Sliders/SimulationSpeedSlider.hpp"
#include "Entities/UIElements/UIElement.hpp"
#include "Entities/UIElements/UIElementIdentifierEnum.hpp"
#include "Scenes/Scene.hpp"
#include "Scenes/SceneIdentifierEnum.hpp"
#include <iostream>

class GameScene : public Scene
{
  public:
    GameScene();
    ~GameScene() = default;
    void populateMe() override;
};