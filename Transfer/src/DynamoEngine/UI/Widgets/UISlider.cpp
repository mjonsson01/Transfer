// File: Transfer/src/DynamoEngine/UI/Widgets/UISlider.cpp

#include "DynamoEngine/UI/Widgets/UISlider.hpp"

// Standard Library Imports
#include <algorithm>
#include <cmath>
#include <utility>

namespace
{
constexpr float TRACK_HEIGHT = 12.0f;
constexpr float KNOB_WIDTH = 20.0f;
constexpr float KNOB_HEIGHT = 30.0f;
constexpr int TICK_COUNT = 30; // a tick sound plays each time the knob crosses one of these marks

constexpr SDL_Color TRACK_COLOR = {128, 128, 128, 255};        // gray
constexpr SDL_Color KNOB_COLOR = {255, 255, 255, 255};         // white
constexpr SDL_Color KNOB_HOVERED_COLOR = {215, 215, 215, 255}; // slightly darker while hovered
} // namespace

namespace DynamoEngine
{
SliderMapping SliderMapping::linear(double min_value, double max_value)
{
    SliderMapping mapping;
    mapping.to_value = [min_value, max_value](double position) { return min_value + position * (max_value - min_value); };
    mapping.to_position = [min_value, max_value](double value) { return (value - min_value) / (max_value - min_value); };
    return mapping;
}

UISlider::UISlider(std::string label, SliderMapping mapping, double starting_value)
    : m_label(std::move(label)), m_mapping(std::move(mapping))
{
    setValue(starting_value);
}

void UISlider::setValue(double value)
{
    m_value = value;
    m_position = std::clamp(m_mapping.to_position(value), 0.0, 1.0);
    m_last_tick = static_cast<int>(std::round(m_position * TICK_COUNT));
}

// --- Geometry: the knob sits at the top of the rect, the track is centered on it, the text goes underneath --- //

SDL_FRect UISlider::trackRect() const
{
    const float track_y = m_rect.y + (KNOB_HEIGHT - TRACK_HEIGHT) / 2.0f;
    return {m_rect.x, track_y, m_rect.w, TRACK_HEIGHT};
}

SDL_FRect UISlider::knobRect() const
{
    const float travel = m_rect.w - KNOB_WIDTH; // how far the knob can move
    const float knob_x = m_rect.x + static_cast<float>(m_position) * travel;
    return {knob_x, m_rect.y, KNOB_WIDTH, KNOB_HEIGHT};
}

void UISlider::draw(UIGeometryBuilder& builder) const
{
    builder.addRect(trackRect(), TRACK_COLOR);
    const bool is_knob_highlighted = m_is_knob_hovered || m_is_dragging;
    builder.addRect(knobRect(), is_knob_highlighted ? KNOB_HOVERED_COLOR : KNOB_COLOR);
    builder.addText(m_label + ": " + std::to_string(m_value), {m_rect.x, m_rect.y + KNOB_HEIGHT});
}

// --- Input --- //

bool UISlider::onMousePressed(Vector2F mouse_position)
{
    m_is_dragging = true;
    moveKnobTo(mouse_position.x_val); // clicking anywhere on the slider jumps the knob there
    return true;
}

void UISlider::onMouseDragged(Vector2F mouse_position)
{
    moveKnobTo(mouse_position.x_val);
}

void UISlider::onMouseReleased(Vector2F mouse_position, bool released_inside)
{
    m_is_dragging = false;
}

void UISlider::onMouseHover(Vector2F mouse_position)
{
    // Same edge rule as UIElement::containsPoint: left and top edges are inside, right and bottom are not
    const SDL_FRect knob = knobRect();
    const bool inside_x = mouse_position.x_val >= knob.x && mouse_position.x_val < knob.x + knob.w;
    const bool inside_y = mouse_position.y_val >= knob.y && mouse_position.y_val < knob.y + knob.h;
    m_is_knob_hovered = inside_x && inside_y;
}

void UISlider::moveKnobTo(float mouse_x)
{
    // Center the knob on the mouse, then convert "how far along the track" into 0..1
    const float travel = m_rect.w - KNOB_WIDTH;
    const float knob_left = mouse_x - KNOB_WIDTH / 2.0f;
    m_position = (travel > 0.0f) ? std::clamp((knob_left - m_rect.x) / travel, 0.0f, 1.0f) : 0.0;

    const double new_value = m_mapping.to_value(m_position);
    const bool value_changed = (new_value != m_value);
    m_value = new_value;

    // A tick sound each time the knob reaches a different tick mark
    const int tick = static_cast<int>(std::round(m_position * TICK_COUNT));
    if (tick != m_last_tick)
    {
        requestSound(UISound::Tick);
        m_last_tick = tick;
    }

    if (value_changed && m_on_value_changed)
    {
        m_on_value_changed(m_value);
    }
}
} // namespace DynamoEngine
