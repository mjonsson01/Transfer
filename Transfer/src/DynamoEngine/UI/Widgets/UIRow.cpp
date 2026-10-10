// File: Transfer/src/DynamoEngine/UI/Widgets/UIRow.cpp

#include "DynamoEngine/UI/Widgets/UIRow.hpp"

// Standard Library Imports
#include <algorithm>
#include <memory>

namespace DynamoEngine
{
void UIRow::updateSize()
{
    // Children first: a row or column nested inside this one must know its size THIS frame, not last frame's
    for (const std::unique_ptr<UIElement>& child : children())
    {
        child->updateSize();
    }

    // Fit the children: widths added up (plus spacing), height of the tallest child
    float total_width = 0.0f;
    float tallest = 0.0f;
    for (const std::unique_ptr<UIElement>& child : children())
    {
        total_width += child->placement().size.x_val;
        tallest = std::max(tallest, child->placement().size.y_val);
    }
    if (!children().empty())
    {
        total_width += m_spacing * static_cast<float>(children().size() - 1);
    }
    m_placement.size = {total_width, tallest};
}

void UIRow::updateLayout(const SDL_FRect& parent_rect)
{
    // 1. Size the row (and, through updateSize, every row and column inside it) to fit its children
    updateSize();

    // 2. Place the row itself inside its parent
    m_rect = m_placement.placeInside(parent_rect);

    // 3. Give each child its slot, left to right
    float slot_x = m_rect.x;
    for (const std::unique_ptr<UIElement>& child : children())
    {
        const float slot_width = child->placement().size.x_val;
        child->updateLayout(SDL_FRect{slot_x, m_rect.y, slot_width, m_rect.h});
        slot_x += slot_width + m_spacing;
    }
}

void UIColumn::updateSize()
{
    // Children first: a row or column nested inside this one must know its size THIS frame, not last frame's
    for (const std::unique_ptr<UIElement>& child : children())
    {
        child->updateSize();
    }

    // Fit the children: heights added up (plus spacing), width of the widest child
    float total_height = 0.0f;
    float widest = 0.0f;
    for (const std::unique_ptr<UIElement>& child : children())
    {
        total_height += child->placement().size.y_val;
        widest = std::max(widest, child->placement().size.x_val);
    }
    if (!children().empty())
    {
        total_height += m_spacing * static_cast<float>(children().size() - 1);
    }
    m_placement.size = {widest, total_height};
}

void UIColumn::updateLayout(const SDL_FRect& parent_rect)
{
    // 1. Size the column (and, through updateSize, every row and column inside it) to fit its children
    updateSize();

    // 2. Place the column itself inside its parent
    m_rect = m_placement.placeInside(parent_rect);

    // 3. Give each child its slot, top to bottom
    float slot_y = m_rect.y;
    for (const std::unique_ptr<UIElement>& child : children())
    {
        const float slot_height = child->placement().size.y_val;
        child->updateLayout(SDL_FRect{m_rect.x, slot_y, m_rect.w, slot_height});
        slot_y += slot_height + m_spacing;
    }
}
} // namespace DynamoEngine
