// File: Transfer/src/Scenes/PauseScene/PauseScene.h

#pragma once

#include "Entities/UIElements/Buttons/PlayGameButton/PlayGameButton.h"
#include "Entities/UIElements/UIElement.h"
#include "Entities/UIElements/UIElementIdentifierEnum.h"
#include "Scenes/Scene.h"
#include "Scenes/SceneIdentifierEnum.h"
#include <iostream>

class PauseScene : public Scene
{
  public:
    PauseScene();
    ~PauseScene() = default;
    void populateMe() override;
};