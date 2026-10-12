// File: Transfer/src/Core/UniverseState.hpp

#pragma once

// Custom Imports
#include "Entities/Physics/GravitationalBody.hpp"
#include "Player/Starship.hpp"

// Standard Library Imports
#include <vector>

// Contains all Physics integratable objects. The state is portable to be copyable
struct UniverseState
{
    std::vector<GravitationalBody> macro_bodies;
    std::vector<GravitationalBody> particles;
    Starship ship;
    int newest_spawn_id = 0; // most recent GravBody spawn ID
};