// File: Transfer/src/DynamoEngine/UI/Widgets/UICheckbox.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/UI/UIElement.hpp"

// Standard Library Imports
#include <functional>
#include <string>

namespace DynamoEngine
{
// A labelled on/off switch: a small square on the left (filled in while checked) and the label to its right, on the
// same gray background as a button. Clicking anywhere on it flips it. Like a button, the click happens on release,
// and only if the mouse is still over it.
class UICheckbox : public UIElement
{
  public:
    using ToggledAction = std::function<void(bool is_checked)>;
    using CheckedSource = std::function<bool()>;

    UICheckbox(std::string label, bool starts_checked);

    // Called with the NEW state every time the player flips it
    void setOnToggled(ToggledAction action) { m_on_toggled = std::move(action); }
    // Where the checked state comes from every frame (e.g. a spawn setting), so the box never shows a stale value
    void setCheckedSource(CheckedSource source) { m_checked_source = std::move(source); }

    void setChecked(bool is_checked) { m_is_checked = is_checked; } // changes the state WITHOUT calling the action
    bool isChecked() const { return m_is_checked; }
    const std::string& label() const { return m_label; }

    // UIElement overrides
    void update(float delta_seconds) override;
    void draw(UIGeometryBuilder& builder) const override;
    bool onMousePressed(Vector2F mouse_position) override;
    void onMouseReleased(Vector2F mouse_position, bool released_inside) override;
    void onMouseEntered() override { m_is_hovered = true; }
    void onMouseExited() override { m_is_hovered = false; }

  private:
    SDL_FRect boxRect() const; // the square on the left, computed from m_rect every time (never stored)

    std::string m_label;
    ToggledAction m_on_toggled;     // empty = nobody is listening
    CheckedSource m_checked_source; // empty = only the player and setChecked() change it
    bool m_is_checked = false;
    bool m_is_hovered = false;
    bool m_is_held_down = false;
};
} // namespace DynamoEngine