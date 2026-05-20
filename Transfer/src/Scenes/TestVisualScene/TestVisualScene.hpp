// File: Transfer/src/Scenes/TestVisualScene/TestVisualScene.h

#pragma once

#include "Entities/UIElements/Buttons/PlayGameButton/PlayGameButton.hpp"
#include "Entities/UIElements/Buttons/ResumeButton/ResumeButton.hpp"
#include "Entities/UIElements/UIElement.hpp"
#include "Entities/UIElements/UIElementIdentifierEnum.hpp"
#include "Scenes/Scene.hpp"
#include "Scenes/SceneIdentifierEnum.hpp"
#include <iostream>

class TestVisualScene : public Scene
{
  public:
    TestVisualScene();
    ~TestVisualScene() = default;
    void populateMe() override;
};