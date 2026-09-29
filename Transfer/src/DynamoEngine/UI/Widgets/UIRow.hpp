// File: Transfer/src/DynamoEngine/UI/Widgets/UIRow.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/UI/UIElement.hpp"

namespace DynamoEngine
{
// Lines its children up left to right with `spacing` points between them.
// A row sizes itself to fit its children, so only its alignment and margin matter when placing it.
// Each child is given a slot as wide as the child and as tall as the row; the child's own alignment
// decides where it sits inside that slot (e.g. CenterLeft = vertically centered).
class UIRow : public UIElement
{
  public:
    explicit UIRow(float spacing) : m_spacing(spacing) {}
    void updateLayout(const SDL_FRect& parent_rect) override;

  private:
    float m_spacing;
};

// Stacks its children top to bottom with `spacing` points between them. The vertical twin of UIRow.
class UIColumn : public UIElement
{
  public:
    explicit UIColumn(float spacing) : m_spacing(spacing) {}
    void updateLayout(const SDL_FRect& parent_rect) override;

  private:
    float m_spacing;
};
} // namespace DynamoEngine
