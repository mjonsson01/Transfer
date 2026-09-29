// File: Transfer/src/DynamoEngine/Input/SDLInputIntake.cpp

#include "DynamoEngine/Input/SDLInputIntake.hpp"

// SDL Imports
#include <SDL3/SDL_mouse.h>
#include <SDL3/SDL_stdinc.h>

// Standard Library Imports
#include <utility>

// Empty namespace for internal linkage only
namespace
{
DynamoEngine::MouseButton translateMouseButton(Uint8 sdl_button)
{
    switch (sdl_button)
    {
    case SDL_BUTTON_LEFT:
        return DynamoEngine::MouseButton::Left;
    case SDL_BUTTON_RIGHT:
        return DynamoEngine::MouseButton::Right;
    case SDL_BUTTON_MIDDLE:
        return DynamoEngine::MouseButton::Middle;
    default:
        return DynamoEngine::MouseButton::None;
    }
}
} // namespace

namespace DynamoEngine
{
bool translateSDLEvent(const SDL_Event& sdl_event, InputEvent& output_event)
{
    InputEvent event;
    switch (sdl_event.type)
    {
    case SDL_EVENT_QUIT:
        event.type = InputEventType::Quit;
        break;

    case SDL_EVENT_WINDOW_RESIZED:
        event.type = InputEventType::WindowResize;
        event.window_width = sdl_event.window.data1;
        event.window_height = sdl_event.window.data2;
        break;

    case SDL_EVENT_WINDOW_FOCUS_LOST:
        event.type = InputEventType::FocusLost;
        break;

    case SDL_EVENT_KEY_DOWN:
    case SDL_EVENT_KEY_UP:
        event.type = (sdl_event.type == SDL_EVENT_KEY_DOWN) ? InputEventType::KeyDown : InputEventType::KeyUp;
        event.key = sdl_event.key.scancode;
        event.is_repeat = sdl_event.key.repeat;
        event.key_modifiers = sdl_event.key.mod;
        break;

    case SDL_EVENT_MOUSE_BUTTON_DOWN:
    case SDL_EVENT_MOUSE_BUTTON_UP:
        event.type = (sdl_event.type == SDL_EVENT_MOUSE_BUTTON_DOWN) ? InputEventType::MouseButtonDown
                                                                     : InputEventType::MouseButtonUp;
        event.mouse_button = translateMouseButton(sdl_event.button.button);
        if (event.mouse_button == MouseButton::None)
        {
            return false; // Button we don't track
        }
        event.type = (sdl_event.type == SDL_EVENT_MOUSE_BUTTON_DOWN) ? InputEventType::MouseButtonDown
                                                                     : InputEventType::MouseButtonUp;
        event.mouse_position = Vector2F(sdl_event.button.x, sdl_event.button.y);
        break;

    case SDL_EVENT_MOUSE_MOTION:
        event.type = InputEventType::MouseMove;
        event.mouse_position = Vector2F(sdl_event.motion.x, sdl_event.motion.y);
        event.mouse_position_delta = Vector2F(sdl_event.motion.xrel, sdl_event.motion.yrel);
        break;

    case SDL_EVENT_MOUSE_WHEEL:
        event.type = InputEventType::MouseWheel;
        event.scroll_delta = sdl_event.wheel.y;
        if (sdl_event.wheel.direction == SDL_MOUSEWHEEL_FLIPPED)
        {
            event.scroll_delta *= -1.0f;
        }
        event.mouse_position = Vector2F(sdl_event.wheel.mouse_x, sdl_event.wheel.mouse_y);
        break;
    case SDL_EVENT_TEXT_INPUT:
        event.type = InputEventType::TextInput;
        event.text_input = sdl_event.text.text;
        break;
    default:
        return false; // Not an input event we care about
        break;
    }
    output_event = std::move(event);
    return true;
}

void SDLInputIntake::pollEvents(std::vector<InputEvent>& output_events)
{
    SDL_Event sdl_event;
    while (SDL_PollEvent(&sdl_event))
    {
        InputEvent event;
        if (translateSDLEvent(sdl_event, event))
        {
            output_events.push_back(std::move(event));
        }
    }
}
} // namespace DynamoEngine