// File: Transfer/src/DynamoEngine/UI/UIElement.cpp

#include "DynamoEngine/UI/UIElement.hpp"

// Standard Library Imports
#include <cassert>
#include <utility>

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
    m_rect = m_placement.placeInside(parent_rect);

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

void UIElement::requestSound(UISound sound)
{
    if (m_parent != nullptr)
    {
        m_parent->requestSound(sound); // pass it up; the root decides what actually happens
    }
}

} // namespace DynamoEngine