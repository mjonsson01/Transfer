// File: Transfer/src/DynamoEngine/UI/UIPlacement.cpp

#include "DynamoEngine/UI/UIPlacement.hpp"

namespace
{
using DynamoEngine::UIAlign;
using DynamoEngine::Vector2F;

// The spot on a rectangle that an alignment refers to, as fractions of its width and height:
// (0, 0) = top-left corner, (0.5, 0.5) = center, (1, 1) = bottom-right corner
Vector2F alignmentSpot(UIAlign align)
{
    switch (align)
    {
    case UIAlign::TopLeft:
        return {0.0f, 0.0f};
    case UIAlign::TopCenter:
        return {0.5f, 0.0f};
    case UIAlign::TopRight:
        return {1.0f, 0.0f};
    case UIAlign::CenterLeft:
        return {0.0f, 0.5f};
    case UIAlign::Center:
        return {0.5f, 0.5f};
    case UIAlign::CenterRight:
        return {1.0f, 0.5f};
    case UIAlign::BottomLeft:
        return {0.0f, 1.0f};
    case UIAlign::BottomCenter:
        return {0.5f, 1.0f};
    case UIAlign::BottomRight:
        return {1.0f, 1.0f};
    case UIAlign::Fill:
        return {0.0f, 0.0f}; // Fill is handled separately in placeInside
    }
    return {0.0f, 0.0f};
}

} // namespace

namespace DynamoEngine
{
SDL_FRect UIPlacement::placeInside(const SDL_FRect& parent_rect) const
{
    if (align == UIAlign::Fill)
    {
        return {parent_rect.x + margin, parent_rect.y + margin, parent_rect.w - 2.0f * margin,
                parent_rect.h - 2.0f * margin};
    }

    const Vector2F spot = alignmentSpot(align);

    // 1. Find the spot on the parent (e.g. the parent's bottom-center)
    const float parent_spot_x = parent_rect.x + spot.x_val * parent_rect.w;
    const float parent_spot_y = parent_rect.y + spot.y_val * parent_rect.h;

    // 2. Put this element's matching spot on it (e.g. my bottom-center on the parent's bottom-center)
    float x = parent_spot_x - spot.x_val * size.x_val;
    float y = parent_spot_y - spot.y_val * size.y_val;

    // 3. Step inward from the edge by the margin: attached to the left edge -> move right, right edge -> move left,
    //    centered -> don't move. (1 - 2 * spot) is +1, 0, or -1 for spots 0, 0.5, and 1.
    x += (1.0f - 2.0f * spot.x_val) * margin;
    y += (1.0f - 2.0f * spot.y_val) * margin;

    return {x, y, size.x_val, size.y_val};
}
} // namespace DynamoEngine
