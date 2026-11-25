// File: Transfer/src/Enties/Physics/GravitationalBodyPair.hpp

#pragma once

#include "Systems/Physics/GravitationalBody.hpp"

struct GravitationalBodyPair
{
    GravitationalBody* heavierBody;
    GravitationalBody* lighterBody;
    double ratio; // heavy/light (>= 1)
    bool equal;
};