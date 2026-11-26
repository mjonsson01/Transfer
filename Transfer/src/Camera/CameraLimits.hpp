// File: Transfer/src/Camera/CameraLimits.hpp

#pragma once

// Custom Imports
#include "Camera/CameraState.hpp"
#include "DynamoEngine/Math/Vector2.hpp"

// Keeps the camera inside the world (world_min_corner .. world_max_corner, world units): the view's center can't move
// so far that the view leaves the world, and when the view is bigger than the world it is centered on it. The star
// layer's own (parallax) pan is kept inside the star field the same way.
// Returns how far the camera was pushed back (world units): (0, 0) when the view was already inside.
DynamoEngine::Vector2D keepViewInsideWorld(CameraState& camera_state, const DynamoEngine::Vector2D& world_min_corner,
                                           const DynamoEngine::Vector2D& world_max_corner);
