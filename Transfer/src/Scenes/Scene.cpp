// File: Transfer/src/Scenes/Scene.cpp
#include "Scene.h"

Scene::Scene(SceneIdentifier sceneID) : sceneIdentifier(sceneID) {}

Scene::~Scene() {};

void Scene::CleanUpSceneElements()
{
    for (auto& [UI_element_ID, UI_element_ptr] : sceneUIElements)
    {
        if (UI_element_ptr)
        {
            // std::cout << "UIElement Deleted" << std::endl;
            delete UI_element_ptr;
            UI_element_ptr = nullptr;
        }
    }
    sceneUIElements.clear();
}