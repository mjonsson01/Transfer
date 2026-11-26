// File: Transfer/src/DynamoEngine/Constants/GlobalConstants.hpp

#pragma once

namespace DynamoEngine
{
// Near-zero comparison epsilon.
constexpr double EPSILON = 1e-8;

// The window size the UI is designed at (UI points). The UI scales up and down from here so that this whole
// 1280 x 720 area always fits inside the window, whatever its shape (see UIRoot::uiScale).
constexpr float UI_REFERENCE_WIDTH = 1280.0f;
constexpr float UI_REFERENCE_HEIGHT = 720.0f;

// Hardcoded Pi
constexpr double PI = 3.141592653589793;

// Two pi
constexpr double TWO_PI = 2.0 * PI;

} // namespace DynamoEngine
