// File: Transfer/src/DynamoEngine/UI/UIElement.cpp

#include "DynamoEngine/UI/UIElement.hpp"

// Standard Library Imports
#include <cassert>
#include <utility>

namespace
{
using DynamoEngine::UIAlign;
using DynamoEngine::UIPlacement;
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

// The rectangle `placement` describes inside `parent_rect`
SDL_FRect placeInside(const UIPlacement& placement, const SDL_FRect& parent_rect)
{
    const float margin = placement.margin;

    if (placement.align == UIAlign::Fill)
    {
        return {parent_rect.x + margin, parent_rect.y + margin, parent_rect.w - 2.0f * margin,
                parent_rect.h - 2.0f * margin};
    }

    const Vector2F spot = alignmentSpot(placement.align);
    const Vector2F size = placement.size;

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
} // namespace

namespace DynamoEngine
{
UIElement& UIElement::addChild(std::unique_ptr<UIElement> child)
{
    assert(child != nullptr && "addChild: child must not be null");

    child->m_parent = this;
    m_children.push_back(std::move(child));
    return *m_children.back();
}

void UIElement::updateLayout(const SDL_FRect& parent_rect)
{
    m_rect = placeInside(m_placement, parent_rect);

    for (const std::unique_ptr<UIElement>& child : m_children)
    {
        child->updateLayout(m_rect); // children are placed inside this element
    }
}

bool UIElement::containsPoint(Vector2F point) const
{
    // Left and top edges count as inside, right and bottom edges don't -- so two elements that touch
    // never both claim the same point
    bool inside_x = point.x_val >= m_rect.x && point.x_val < m_rect.x + m_rect.w;
    bool inside_y = point.y_val >= m_rect.y && point.y_val < m_rect.y + m_rect.h;
    return inside_x && inside_y;
}

UILayer UIElement::layer() const
{
    if (m_layer.has_value())
    {
        return m_layer.value(); // set explicitly on this element
    }
    if (m_parent != nullptr)
    {
        return m_parent->layer(); // inherit (keeps asking up the tree until someone has one set)
    }
    return UILayer::HUD; // the root, with nothing set
}
} // namespace DynamoEngine