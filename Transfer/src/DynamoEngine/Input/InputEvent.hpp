// File: Transfer/src/DynamoEngine/Input/InputEvent.hpp

#pragma once

// SDL3 Imports
#include <SDL3/SDL_keycode.h>
#include <SDL3/SDL_scancode.h>

// Custom Imports
#include "DynamoEngine/Math/Vector2.hpp"

// Standard Library Imports
#include <cstdint>
#include <string>

namespace DynamoEngine
{
using Key = SDL_Scancode;
enum class MouseButton : uint8_t
{
    None = 0,
    Left = 1,
    Right = 2,
    Middle = 3,
    Count = 4
};

enum class InputEventType : uint8_t
{
    None = 0,
    KeyDown = 1,
    KeyUp = 2,
    MouseButtonDown = 3,
    MouseButtonUp = 4,
    MouseMove = 5,
    MouseWheel = 6,
    TextInput = 7,
    WindowResize = 8,
    FocusLost = 9,
    Quit = 10
};

struct InputEvent
{
    InputEventType type = InputEventType::None;

    // Keyboard Input
    Key key = SDL_Scancode::SDL_SCANCODE_UNKNOWN;
    bool is_repeat = false;       // Track whether we are holding the key (not a new press)
    SDL_Keymod key_modifiers = 0; // holding shift, alt, control

    // Mouse Input
    MouseButton mouse_button = MouseButton::None;
    Vector2F mouse_position;       // Current mouse position
    Vector2F mouse_position_delta; // The movement delta since the previous MouseMove Event
    float scrollDelta = 0.0f;      // positive is scroll 'up', negative is scroll 'down'

    // Text Input
    std::string text_input;

    // Window Resize Vars, in screen-space points
    int window_width = 0;
    int window_height = 0;
};
} // namespace DynamoEngine