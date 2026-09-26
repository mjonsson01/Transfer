// File: Transfer/src/DynamoEngine/UI/UIPlacement.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/Math/Vector2.hpp"

// Standard Library Imports
#include <cstdint>

namespace DynamoEngine
{
// Which spot of the parent an element attaches to. The element's matching spot sits on it:
// BottomCenter = "my bottom-center on my parent's bottom-center". Fill = cover the whole parent.
enum class UIAlign : uint8_t
{
    TopLeft,
    TopCenter,
    TopRight,
    CenterLeft,
    Center,
    CenterRight,
    BottomLeft,
    BottomCenter,
    BottomRight,
    Fill,
};

struct UIPlacement
{
    UIAlign align = UIAlign::TopLeft;
    Vector2F size = {0.0f, 0.0f}; // in UI points (ignored by Fill)
    float margin = 0.0f;          // gap between the element and the parent's edge
};
} // namespace DynamoEngine