// File: Transfer/src/Core/SpawnSettings.hpp

#pragma once

// What the next planet or cluster the player spawns will be like. The spawn panel's checkboxes edit it and
// PhysicsSystem's create functions read it. Unlike DEPRECATED_InputState's isCreating... flags, nothing resets it
// after a spawn: a setting stays until the player changes it.
struct SpawnSettings
{
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
};