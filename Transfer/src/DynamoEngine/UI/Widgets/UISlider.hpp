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
// The knob darkens slightly while the cursor is over the knob itself, and while it is being dragged.
class UISlider : public UIElement
{
  public:
    using ValueChangedAction = std::function<void(double new_value)>;

    UISlider(std::string label, SliderMapping mapping, double starting_value);

    void setOnValueChanged(ValueChangedAction action) { m_on_value_changed = std::move(action); }
    void setValue(double value); // moves the knob without calling the value-changed action
    double value() const { return m_value; }
    // Makes the value snap to whole multiples of `step` (e.g. 0.1: 0.1, 0.2, 0.3, ...), shown with `decimal_places`
    // decimals. The knob jumps from step to step. Without it the value is continuous, shown by std::to_string.
    // (Pick a range whose ends are multiples of the step.)
    void setStep(double step, int decimal_places);
    // The "Label: value" line drawn under the track
    std::string text() const;

    void draw(UIGeometryBuilder& builder) const override;
    bool onMousePressed(Vector2F mouse_position) override;
    void onMouseDragged(Vector2F mouse_position) override;
    void onMouseReleased(Vector2F mouse_position, bool released_inside) override;
    void onMouseHover(Vector2F mouse_position) override;
    void onMouseExited() override { m_is_knob_hovered = false; }

  private:
    void moveKnobTo(float mouse_x);
    // `value` rounded to the nearest whole multiple of m_step (unchanged without a step)
    double snapToStep(double value) const;
    SDL_FRect trackRect() const;
    SDL_FRect knobRect() const;

    std::string m_label;
    SliderMapping m_mapping;
    ValueChangedAction m_on_value_changed; // empty = nobody is listening
    double m_value = 0.0;
    double m_position = 0.0;   // 0..1, where the knob is along the track
    double m_step = 0.0;       // 0 = continuous (setStep)
    int m_decimal_places = -1; // -1 = std::to_string's format
    int m_last_tick = -1;    // the tick mark the knob was last at (for the tick sound)

    // The knob is highlighted while the cursor is on the knob itself, or while it's being dragged
    bool m_is_knob_hovered = false;
    bool m_is_dragging = false;
};
} // namespace DynamoEngine
