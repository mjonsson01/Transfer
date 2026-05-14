// File: Transfer/src/Scenes/GameScene/GameScene.h

#pragma once

#include "Entities/UIElements/Overlay/FPSCounter.h"
#include "Entities/UIElements/Sliders/MassSlider.h"
#include "Entities/UIElements/Sliders/RadiusSlider.h"
#include "Entities/UIElements/Sliders/SimulationSpeedSlider.h"
#include "Entities/UIElements/UIElement.h"
#include "Entities/UIElements/UIElementIdentifierEnum.h"
#include "Scenes/Scene.h"
#include "Scenes/SceneIdentifierEnum.h"
#include <iostream>

class GameScene : public Scene
{
  public:
    GameScene();
    ~GameScene() = default;
    void populateMe() override;
};