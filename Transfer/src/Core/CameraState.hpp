// File: Transfer/src/Core/CameraState.hpp

#pragma once

// Custom imports
#include "DynamoEngine/Math/Vector2.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"

struct CameraState
{
    double zoom = STARTUP_ZOOM_VALUE;
    DynamoEngine::Vector2D offset = {0.0, 0.0}; // pan offset
    DynamoEngine::Vector2D twinkling_star_offset = {0.0, 0.0};

    float window_width = static_cast<float>(SCREEN_WIDTH);
    float window_height = static_cast<float>(SCREEN_HEIGHT);

    float maxDisplayWidth = static_cast<float>(SCREEN_WIDTH);
    float maxDisplayHeight = static_cast<float>(SCREEN_HEIGHT);

    float renderAlpha;
};