// File: Transfer/src/Core/PlayerInput.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/Math/Vector2.hpp"

// What the player is doing this frame, in Transfer's terms. InputSystem writes it every frame (from the engine's
// InputState, after the UI took its share); the other systems only read it. Game owns it and passes it explicitly.

// The ship's controls (W/S, A/D)
struct ShipControls
{
    // +1 = W (throttle up), -1 = S (throttle down), 0 = neither, or both. The throttle itself is the ship's state:
    // it moves while a key is held and stays where it is when released (Starship::applyThrust)
    int thrust = 0;
    // +1 = A (counterclockwise on screen), -1 = D (clockwise), 0 = neither, or both
    int rotation = 0;
};

// Where and how the next body would spawn, already in WORLD space: computed once, by InputSystem, so the spawn
// itself, the see-through preview and the placement trajectory can never disagree
struct SpawnAim
{
    bool is_cluster = false;    // right button (left = a planet)
    bool with_velocity = false; // Shift: the body starts where the drag began, moving toward the cursor
    DynamoEngine::Vector2D world_position = {0.0, 0.0};
    DynamoEngine::Vector2D world_velocity = {0.0, 0.0}; // 0 without Shift
};

// One-frame requests. InputSystem clears them at the start of every frame, so they last exactly one frame and
// nobody else ever has to reset them.
struct PlayerCommands
{
    // A spawn drag was released on the world this frame: spawn what spawn_aim describes
    bool spawn = false;
    // Only meaningful when spawn is true
    SpawnAim spawn_aim;
    // Delete
    bool undo_newest_spawn = false;
    // Shift + Delete
    bool clear_all_bodies = false;
    // I, while planning
    bool toggle_trajectories = false;
};

struct PlayerInput
{
    // A spawn button held down on the world (not on the UI): what would spawn on release, for the previews
    bool is_aiming_spawn = false;
    SpawnAim spawn_aim; // only meaningful while is_aiming_spawn is true

    PlayerCommands commands;
    ShipControls ship;
};
