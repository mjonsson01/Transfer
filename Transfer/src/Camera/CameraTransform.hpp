// File: Transfer/src/Utilities/Rendering/CameraTransform.hpp

#pragma once

#include "Core/CameraState.hpp"
#include "DynamoEngine/Math/Vector2.hpp"

inline DynamoEngine::Vector2D ScreenToWorldCoordinates(const DynamoEngine::Vector2D& screenPoint,
                                                       const CameraState& camera_state)
{
    return (screenPoint / camera_state.zoom) - camera_state.offset;
}