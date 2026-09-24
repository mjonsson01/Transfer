// File: Transfer/src/Scenes/TestVisualScene/TestVisualScene.h

#pragma once


#include "Entities/UIElements/UIElement.hpp"
#include "Scenes/Scene.hpp"
#include "Scenes/SceneIdentifierEnum.hpp"
#include "Entities/UIElements/DropDownMenu/DropDownMenu.hpp"
#include <iostream>

class TestVisualScene : public Scene
{
  public:
    TestVisualScene();
    ~TestVisualScene() = default;
    void populateMe() override;
};