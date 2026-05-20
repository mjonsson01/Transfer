// File: Transfer/src/Scenes/Scene.h

#pragma once

#include "Entities/UIElements/UIElement.hpp"
#include "Entities/UIElements/UIElementIdentifierEnum.hpp"
#include "Scenes/SceneIdentifierEnum.hpp"

#include <unordered_map>

class Scene
{
  public:
    Scene(SceneIdentifier sceneID);
    virtual ~Scene();
    virtual void populateMe() {}; // should be overridden by each child class.
    void CleanUpSceneElements();
    std::unordered_map<UIElementIdentifier, UIElement*>& getSceneElements() { return sceneUIElements; }

  protected:
    SceneIdentifier sceneIdentifier;
    std::unordered_map<UIElementIdentifier, UIElement*> sceneUIElements;
};