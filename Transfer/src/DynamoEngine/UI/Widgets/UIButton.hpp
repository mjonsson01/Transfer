// File: Transfer/src/DynamoEngine/UI/Widgets/UIButton.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/UI/UIAction.hpp"
#include "DynamoEngine/UI/UIElement.hpp"

// Standard Library Imports
#include <string>

namespace DynamoEngine
{
// A clickable rectangle with centered text. It darkens while hovered and darkens further while held down.
// The click happens on release, and only if the mouse is still over the button (so you can "back out" of a click).
class UIButton : public UIElement
{
  public:
    explicit UIButton(std::string text) : m_text(std::move(text)) {}

    void setOnClick(UIAction action) { m_on_click = std::move(action); }
    void setText(std::string text) { m_text = std::move(text); }
    const std::string& text() const { return m_text; }

    void draw(UIGeometryBuilder& builder) const override;
    bool onMousePressed(Vector2F mouse_position) override;
    void onMouseReleased(Vector2F mouse_position, bool released_inside) override;
    void onMouseEntered() override { m_is_hovered = true; }
    void onMouseExited() override { m_is_hovered = false; }

  private:
    std::string m_text;
    UIAction m_on_click; // empty = the button does nothing when clicked
    bool m_is_hovered = false;
    bool m_is_held_down = false;
};
} // namespace DynamoEngine
