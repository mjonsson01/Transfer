// File: Transfer/src/Systems/Physics/UniverseState.hpp

#pragma once

// Custom Imports
#include "Player/Starship.hpp"
#include "Systems/Physics/GravitationalBody.hpp"

// Standard Library Imports
#include <vector>

// Contains all Physics integratable objects. The state is portable to be copyable
struct UniverseState
{
    std::vector<GravitationalBody> macro_bodies;
    std::vector<GravitationalBody> particles;
    Starship ship;
    // Whether the ship is in this universe at all. Without it, nothing touches `ship`: no gravity, collisions,
    // controls, prediction line or drawing. True by default; the game sets it from RUN_GAME_WITH_PLAYER_SHIP.
    bool has_ship = true;
    int newest_spawn_id = 0; // most recent GravBody spawn ID

    // The world's edges, in world pixels. The camera stays inside them; with reflective_walls on, every body (and the
    // ship) bounces off them. The game sets them at start-up (the star field's extent).
    DynamoEngine::Vector2D world_min_corner = {-1.0e9, -1.0e9};
    DynamoEngine::Vector2D world_max_corner = {1.0e9, 1.0e9};
    bool reflective_walls = false;

    // An optional collision box (the test scene's): its four edges are walls. A body whose center was inside stays
    // inside; a body coming from outside bounces off the outside. Drawn as lines, independent of the world's edges.
    bool has_collision_box = false;
    DynamoEngine::Vector2D collision_box_min_corner = {0.0, 0.0};
    DynamoEngine::Vector2D collision_box_max_corner = {0.0, 0.0};
};