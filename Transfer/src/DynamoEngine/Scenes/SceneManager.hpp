// File: Transfer/src/DynamoEngine/Scenes/SceneManager.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/Scenes/Scene.hpp"

// Standard Library Imports
#include <memory>
#include <optional>
#include <unordered_map>

namespace DynamoEngine
{
// Owns every scene and knows which one is active.
//
// Switching is deferred: requestSwitch() only remembers the request, and applyPendingSwitch() (called once at the
// end of each frame) performs it. That way a button can switch scenes from inside its own click handler without
// the scene it belongs to being torn down mid-click.
class SceneManager
{
  public:
    // Registers a scene under `id` and takes ownership of it. Registering the same id twice is a mistake:
    // it is logged, stops Debug builds, and is ignored in Release builds.
    void addScene(SceneId id, std::unique_ptr<Scene> scene);

    bool hasScene(SceneId id) const { return m_scenes.contains(id); }

    // Asks to switch to scene `id` at the end of this frame. The latest request in a frame wins.
    // An unknown id is logged, stops Debug builds, and is ignored in Release builds (the game stays where it is).
    void requestSwitch(SceneId id);

    // Performs the requested switch, if there is one: the old scene's UI lets go of any press or hover,
    // then onExit() on the old scene, then onEnter() on the new one.
    // Call once per frame, after everything else. Requesting the scene that is already active does nothing.
    void applyPendingSwitch();

    // --- The active scene --- //
    bool hasCurrentScene() const { return m_current_scene_id.has_value(); }
    SceneId currentSceneId() const; // only valid when hasCurrentScene()
    Scene& currentScene();          // only valid when hasCurrentScene()
    const Scene& currentScene() const;

  private:
    std::unordered_map<SceneId, std::unique_ptr<Scene>> m_scenes;
    std::optional<SceneId> m_current_scene_id; // empty until the first switch
    std::optional<SceneId> m_pending_scene_id; // empty when no switch has been requested this frame
};
} // namespace DynamoEngine