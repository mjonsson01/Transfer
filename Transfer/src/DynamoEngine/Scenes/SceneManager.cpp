// File: Transfer/src/DynamoEngine/Scenes/SceneManager.cpp

#include "DynamoEngine/Scenes/SceneManager.hpp"

// SDL Imports
#include <SDL3/SDL_log.h>

// Standard Library Imports
#include <cassert>
#include <utility>

namespace DynamoEngine
{
void SceneManager::addScene(SceneId id, std::unique_ptr<Scene> scene)
{
    if (hasScene(id))
    {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SceneManager: scene id %u was registered twice", id);
        assert(false && "addScene: duplicate scene id");
        return;
    }
    m_scenes[id] = std::move(scene);
}

void SceneManager::requestSwitch(SceneId id)
{
    if (!hasScene(id))
    {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SceneManager: no scene is registered with id %u", id);
        assert(false && "requestSwitch: unknown scene id");
        return;
    }
    m_pending_scene_id = id;
}

void SceneManager::applyPendingSwitch()
{
    if (!m_pending_scene_id.has_value())
    {
        return; // nothing requested this frame
    }
    const SceneId next_scene_id = m_pending_scene_id.value();
    m_pending_scene_id.reset();

    if (m_current_scene_id == next_scene_id)
    {
        return; // already there
    }

    if (hasCurrentScene())
    {
        currentScene().onExit();
    }
    m_current_scene_id = next_scene_id;
    currentScene().onEnter();
}

SceneId SceneManager::currentSceneId() const
{
    assert(hasCurrentScene() && "currentSceneId: no scene is active yet");
    return m_current_scene_id.value();
}

Scene& SceneManager::currentScene()
{
    assert(hasCurrentScene() && "currentScene: no scene is active yet");
    return *m_scenes.at(m_current_scene_id.value());
}

const Scene& SceneManager::currentScene() const
{
    assert(hasCurrentScene() && "currentScene: no scene is active yet");
    return *m_scenes.at(m_current_scene_id.value());
}
} // namespace DynamoEngine