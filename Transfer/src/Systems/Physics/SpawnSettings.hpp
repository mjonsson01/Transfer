// File: Transfer/src/Core/SpawnSettings.hpp

#pragma once

// The spawn panel's on/off settings, by name (for SpawnSettings::setFlag)
enum class SpawnFlag
{
    Static,
    Shatterable,
    Accretable,
    Bounce,
    Ghost,
};

// What the next planet or cluster the player spawns will be like. The spawn panel (sliders and checkboxes) edits it and
// PhysicsSystem's create functions read it. Nothing resets it after a spawn: a setting stays until the player
// changes it.
struct SpawnSettings
{
    // The radius and mass sliders (0 until the player moves them: nothing spawns with radius <= 1 or zero mass)
    double radius = 0.0;
    double mass = 0.0;

    // Force-static (GravitationalBody::isForceStatic): never pulled by gravity, an immovable wall in bounces, keeps its
    // velocity forever, and only ever absorbs much lighter bodies. Clusters ignore this setting
    bool is_force_static = false;
    // May shatter into debris when hit hard enough (GravitationalBody::isShatterable). Planets only: particles never
    // shatter.
    bool is_shatterable = true;

    // May be absorbed by a much heavier body (GravitationalBody::isAccretable). For a cluster: every particle.
    bool is_accretable = true;

    // Only ever bounces: never shatters, absorbs or gets absorbed, whatever it hits (GravitationalBody::isBounce).
    // For a cluster: every particle.
    bool is_bounce = false;

    // A ghost passes through everything (planets, debris, the ship) but still pulls and is pulled by gravity
    // (GravitationalBody::isCollidable = false). For a cluster: every particle.
    bool is_ghost = false;

    // Switches one on/off setting (the spawn panel's checkboxes all go through here). Switching one ON switches off
    // the settings it makes meaningless, so the panel never shows a combination that doesn't do what it says:
    //   Ghost (never collides)    -> off: Bounce, Accretable, Shatterable (they only act in collisions)
    //   Bounce (only ever bounces) -> off: Ghost, Accretable, Shatterable (a bounce never merges or shatters)
    //   Accretable, Shatterable    -> off: Ghost, Bounce
    //   Static                     -> works with everything
    // Switching a setting OFF never switches anything back on.
    void setFlag(SpawnFlag flag, bool on);
    // Whether that setting is on (what its checkbox shows)
    bool flag(SpawnFlag flag) const;
};