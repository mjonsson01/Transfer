// File: Transfer/src/Scenes/TestVisualScene/TestVisualScene.h

#pragma once

#include "Entities/UIElements/Buttons/PlayGameButton/PlayGameButton.h"
#include "Entities/UIElements/Buttons/ResumeButton/ResumeButton.h"
#include "Entities/UIElements/UIElement.h"
#include "Entities/UIElements/UIElementIdentifierEnum.h"
#include "Scenes/Scene.h"
#include "Scenes/SceneIdentifierEnum.h"
#include <iostream>

class TestVisualScene : public Scene
{
  public:
    TestVisualScene();
    ~TestVisualScene() = default;
    void populateMe() override;
};