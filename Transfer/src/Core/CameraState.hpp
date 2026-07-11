// File: Transfer/src/Core/CameraState.hpp

#pragma once

// Custom imports
#include "Utilities/Constants/GameSystemConstants.hpp"
#include "Utilities/Math/Vector2D.hpp"

struct CameraState
{
    double zoom = STARTUP_ZOOM_VALUE;
    Vector2D offset = {0.0, 0.0}; // pan offset
    Vector2D twinklingStarOffset = {0.0, 0.0};

    float windowWidth = static_cast<float>(SCREEN_WIDTH);
    float windowHeight = static_cast<float>(SCREEN_HEIGHT);

    float maxDisplayWidth = static_cast<float>(SCREEN_WIDTH);
    float maxDisplayHeight = static_cast<float>(SCREEN_HEIGHT);

    float renderAlpha;
};