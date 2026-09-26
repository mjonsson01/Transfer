// File: Transfer/src/DynamoEngine/UI/Widgets/UISlider.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/UI/UIElement.hpp"

// Standard Library Imports
#include <functional>
#include <string>

namespace DynamoEngine
{
// Converts between where the knob is (0 = far left, 1 = far right) and the value it stands for.
// Needs both directions: dragging uses to_value, placing the knob for a given value uses to_position.
struct SliderMapping
{
    std::function<double(double position)> to_value;
    std::function<double(double value)> to_position;

    // A straight line from min_value (far left) to max_value (far right)
    static SliderMapping linear(double min_value, double max_value);
};

// A horizontal slider: a track, a knob, and a "Label: value" line underneath.
// The whole rect is clickable; pressing anywhere jumps the knob there, then it follows the mouse until release.
class UISlider : public UIElement
{
  public:
    using ValueChangedAction = std::function<void(double new_value)>;

    UISlider(std::string label, SliderMapping mapping, double starting_value);

    void setOnValueChanged(ValueChangedAction action) { m_on_value_changed = std::move(action); }
    void setValue(double value); // moves the knob without calling the value-changed action
    double value() const { return m_value; }

    void draw(UIGeometryBuilder& builder) const override;
    bool onMousePressed(Vector2F mouse_position) override;
    void onMouseDragged(Vector2F mouse_position) override;
    void onMouseEntered() override { m_is_hovered = true; }
    void onMouseExited() override { m_is_hovered = false; }

  private:
    void moveKnobTo(float mouse_x);
    SDL_FRect trackRect() const;
    SDL_FRect knobRect() const;

    std::string m_label;
    SliderMapping m_mapping;
    ValueChangedAction m_on_value_changed; // empty = nobody is listening
    double m_value = 0.0;
    double m_position = 0.0; // 0..1, where the knob is along the track
    int m_last_tick = -1;    // the tick mark the knob was last at (for the tick sound)
    bool m_is_hovered = false;
};
} // namespace DynamoEngine
