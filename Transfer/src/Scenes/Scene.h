// File: Transfer/src/Scenes/Scene.h

#pragma once

#include "Entities/UIElements/UIElement.h"
#include "Entities/UIElements/UIElementIdentifierEnum.h"
#include "Scenes/SceneIdentifierEnum.h"

#include <unordered_map>

class Scene
{
  public:
    Scene(SceneIdentifier sceneID);
    virtual ~Scene();
    void CleanUpSceneElements();

  private:
    SceneIdentifier sceneIdentifier;
    std::unordered_map<UIElementIdentifier, UIElement*> sceneUIElements;
};