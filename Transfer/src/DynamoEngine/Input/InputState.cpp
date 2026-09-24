// File: Transfer/src/DynamoEngine/Input/InputState.cpp

#include "DynamoEngine/Input/InputState.hpp"

namespace DynamoEngine
{
void InputState::beginInputFrame()
{
    m_keys_pressed_this_frame.reset();
    m_keys_released_this_frame.reset();
    m_mouse_buttons_pressed_this_frame.reset();
    m_mouse_buttons_released_this_frame.reset();
    m_mouse_position_delta = Vector2F();
    m_mouse_scroll_delta = 0.0f;
}
void InputState::applyInputEvent(const InputEvent& event)
{
    switch (event.type)
    {
    case InputEventType::KeyDown:
        if (!event.is_repeat) // auto repeat driven by the OS should not be counted as a new press
        {
            m_keys_pressed_this_frame.set(index(event.key));
        }
        m_keys_down.set(index(event.key));
        break;

    case InputEventType::KeyUp:
        m_keys_down.reset(index(event.key));
        m_keys_released_this_frame.set(index(event.key));
        break;

    case InputEventType::MouseButtonDown:
        m_mouse_position = event.mouse_position;
        m_mouse_buttons_down.set(index(event.mouse_button));
        m_mouse_button_pressed_positions[index(event.mouse_button)] = event.mouse_position;
        m_mouse_buttons_pressed_this_frame.set(index(event.mouse_button));
        break;

    case InputEventType::MouseButtonUp:
        m_mouse_position = event.mouse_position;
        m_mouse_buttons_down.reset(index(event.mouse_button));
        m_mouse_buttons_released_this_frame.set(index(event.mouse_button));
        break;

    case InputEventType::MouseMove:
        m_mouse_position = event.mouse_position;
        m_mouse_position_delta += event.mouse_position_delta;
        break;

    case InputEventType::MouseWheel:
        m_mouse_scroll_delta += event.scroll_delta;
        break;

    case InputEventType::FocusLost:
        m_keys_down.reset();
        m_mouse_buttons_down.reset();
        break;

    case InputEventType::Quit:
        m_quit_requested = true;
        break;

    default:
        break;
    }
}

} // namespace DynamoEngine
