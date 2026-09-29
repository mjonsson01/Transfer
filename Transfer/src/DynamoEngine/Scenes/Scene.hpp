// File: Transfer/src/DynamoEngine/Scenes/Scene.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/UI/UIRoot.hpp"

// Standard Library Imports
#include <cstdint>

namespace DynamoEngine
{
// A number that names a scene. The game defines its own list of names, e.g.
//     namespace TransferScene { enum : DynamoEngine::SceneId { StartMenu, Game, Pause }; }
using SceneId = uint32_t;

// What kind of scene this is: which engine systems run while it is active
struct SceneSettings
{
    bool runs_simulation = false; // physics ticks while this scene is active
    bool draws_world = false;     // bodies, stars, and ship are rendered behind the UI

    static SceneSettings menu() { return {.runs_simulation = false, .draws_world = false}; }
    static SceneSettings simulation() { return {.runs_simulation = true, .draws_world = true}; }
};

// One screen of the game: its settings plus the UI it shows.
// Simple scenes need no subclass; override onEnter/onExit/update only when a scene needs its own behavior.
class Scene
{
  public:
    explicit Scene(SceneSettings settings) : m_settings(settings) {}
    virtual ~Scene() = default;

    const SceneSettings& settings() const { return m_settings; }

    UIRoot& ui() { return m_ui; }
    const UIRoot& ui() const { return m_ui; }

    // --- Override these --- //
    virtual void onEnter() {}                   // this scene just became active
    virtual void onExit() {}                    // this scene is about to stop being active
    virtual void update(float delta_seconds) {} // once per frame while active

  private:
    SceneSettings m_settings;
    UIRoot m_ui;
};
} // namespace DynamoEngine