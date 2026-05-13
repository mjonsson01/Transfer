// File: Transfer/src/Scenes/GameScene/GameScene.h

#pragma once

#include "Entities/UIElements/UIElement.h"
#include "Entities/UIElements/UIElementIdentifierEnum.h"
#include "Scenes/Scene.h"
#include "Scenes/SceneIdentifierEnum.h"

class GameScene : public Scene
{
  public:
    GameScene();
    ~GameScene() = default;
    void populateMe() override;
};