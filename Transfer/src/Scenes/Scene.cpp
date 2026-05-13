// File: Transfer/src/Scenes/Scene.cpp
#include "Scene.h"

Scene::Scene(SceneIdentifier sceneID) : sceneIdentifier(sceneID) {}

Scene::~Scene() {};
void Scene::CleanUpSceneElements()
{
    for (auto& scene_element : sceneUIElements)
    {
        delete scene_element.second;
    }
    sceneUIElements.clear();
}