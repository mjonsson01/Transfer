// File: Transfer/src/Utilities/Rendering/CameraTransform.hpp

#pragma once

#include "Core/CameraState.hpp"
#include "Utilities/Math/Vector2D.hpp"

inline Vector2D ScreenToWorldCoordinates(const Vector2D& screenPoint, const CameraState& camera_state)
{
    return (screenPoint / camera_state.zoom) - camera_state.offset;
}