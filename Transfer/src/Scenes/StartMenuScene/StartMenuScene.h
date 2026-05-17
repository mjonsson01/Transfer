// File: Transfer/src/Scenes/StartMenuScene/StartMenuScene.h

#pragma once

#include "Entities/UIElements/Buttons/PlayGameButton/PlayGameButton.h"
#include "Entities/UIElements/Overlay/FPSCounter.h"
#include "Entities/UIElements/Sliders/MassSlider.h"
#include "Entities/UIElements/Sliders/RadiusSlider.h"
#include "Entities/UIElements/Sliders/SimulationSpeedSlider.h"
#include "Entities/UIElements/UIElement.h"
#include "Entities/UIElements/UIElementIdentifierEnum.h"
#include "Scenes/Scene.h"
#include "Scenes/SceneIdentifierEnum.h"
#include <iostream>

class StartMenuScene : public Scene
{
  public:
    StartMenuScene();
    ~StartMenuScene() = default;
    void populateMe() override;
};