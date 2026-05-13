// File: Transfer/src/Scenes/Scene.cpp
#include "Scene.h"

Scene::Scene(SceneIdentifier sceneID) : sceneIdentifier(sceneID) {}

Scene::~Scene() {};

void Scene::CleanUpSceneElements()
{
    for (auto& [id, UI_element_ptr] : sceneUIElements)
    {
        if (UI_element_ptr)
        {
            delete UI_element_ptr;
            UI_element_ptr = nullptr;
        }
    }
    sceneUIElements.clear();
}