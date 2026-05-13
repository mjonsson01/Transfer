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
    virtual void populateMe() {}; // should be overridden by each child class.
    void CleanUpSceneElements();

  protected:
    SceneIdentifier sceneIdentifier;
    std::unordered_map<UIElementIdentifier, UIElement*> sceneUIElements;
};