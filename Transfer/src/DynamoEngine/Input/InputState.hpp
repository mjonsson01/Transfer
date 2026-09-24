// File: Transfer/src/DynamoEngine/Input/InputState.hpp

#pragma once

// SDL Imports
#include <SDL3/SDL_scancode.h> // IWYU pragma: export

// Custom Imports
#include "DynamoEngine/Input/InputEvent.hpp" // IWYU pragma: export
#include "DynamoEngine/Math/Vector2.hpp"     // IWYU pragma: export

// Standard Library Imports
#include <array>
#include <bitset>
#include <cstddef>

namespace DynamoEngine
{
// Current State of Keyboard and Mouse

class InputState
{
  public:
    InputState() = default;
    ~InputState() = default;

  public:
    // --- Top Level Methods --- //
    void beginInputFrame(); // cleans up all old frame vars, called before applying a new event on the following frame
    void applyInputEvent(const InputEvent& event); // fold a single input event into the current state

    // --- Keyboard --- //

    // Checks to see current status of the key this frame
    bool isKeyDown(Key key) const { return m_keys_down[index(key)]; }

    // Checks to see if the key was pressed this frame
    bool wasKeyPressed(Key key) const { return m_keys_pressed_this_frame[index(key)]; }

    // Checks to see if the key was released this frame
    bool wasKeyReleased(Key key) const { return m_keys_released_this_frame[index(key)]; }

    // --- Helpers for KeyBoard Modifiers --- //

    // Is either shift key pressed
    bool isShiftDown() const { return isKeyDown(SDL_SCANCODE_LSHIFT) || isKeyDown(SDL_SCANCODE_RSHIFT); }

    // Is either control key pressed
    bool isCtrlDown() const { return isKeyDown(SDL_SCANCODE_RCTRL) || isKeyDown(SDL_SCANCODE_LCTRL); }

    // Is either alt key pressed
    bool isAltDown() const { return isKeyDown(SDL_SCANCODE_LALT) || isKeyDown(SDL_SCANCODE_RALT); }

    // Is Cmd (on mac) or Win key (on windows) pressed
    bool isSuperDown() const { return isKeyDown(SDL_SCANCODE_LGUI) || isKeyDown(SDL_SCANCODE_RGUI); }

    // --- Mouse Buttons --- //

    // Checks to see current status of the mouse button this frame
    bool isMouseButtonDown(MouseButton button) const { return m_mouse_buttons_down[index(button)]; }

    // Checks to see if the mouse was pressed this frame
    bool wasMouseButtonPressed(MouseButton button) const { return m_mouse_buttons_pressed_this_frame[index(button)]; }

    // Checks to see if the mouse was released this frame
    bool wasMouseButtonReleased(MouseButton button) const { return m_mouse_buttons_released_this_frame[index(button)]; }

    // Checks if any mouse buttons are pressed this frame
    bool isAnyMouseButtonDown() const { return m_mouse_buttons_down.any(); }

    // --- Mouse position/motion in screen space --- //

    // Returns the current mouse position this frame
    const Vector2F& mousePosition() const { return m_mouse_position; }

    // Returns the position this button was most recently pressed
    const Vector2F& pressPosition(MouseButton button) const { return m_mouse_button_pressed_positions[index(button)]; }

    // Returns the net position delta the mouse has accumulated this frame
    const Vector2F& mousePositionDeltaThisFrame() const { return m_mouse_position_delta; }

    // Returns the net scroll delta the mouse has accumulated this frame
    float mouseScrollDeltaThisFrame() const { return m_mouse_scroll_delta; }

    // --- Application Level Control --- //

    // Fetches the quit requested status
    bool quitRequested() const { return m_quit_requested; }

  private:
    // index() functions to map types into size_t values to iterate through the key and mouse button maps
    static constexpr std::size_t index(MouseButton button) { return static_cast<std::size_t>(button); }

    // index() functions to map types into size_t values to iterate through the key and mouse button maps
    static constexpr std::size_t index(Key key) { return static_cast<std::size_t>(key); }

    // Easy lookup for Mouse Button Count
    static constexpr std::size_t MOUSE_BUTTON_COUNT = static_cast<std::size_t>(MouseButton::Count);

    // Lookup structures for keys and buttons
    std::bitset<SDL_SCANCODE_COUNT> m_keys_down;                        // Keys currently down this input frame
    std::bitset<SDL_SCANCODE_COUNT> m_keys_pressed_this_frame;          // Keys pressed this specific frame
    std::bitset<SDL_SCANCODE_COUNT> m_keys_released_this_frame;         // Keys released this specific frame
    std::bitset<MOUSE_BUTTON_COUNT> m_mouse_buttons_down;               // Mouse buttons currently down this input frame
    std::bitset<MOUSE_BUTTON_COUNT> m_mouse_buttons_pressed_this_frame; // Mouse buttons pressed this specific frame
    std::bitset<MOUSE_BUTTON_COUNT> m_mouse_buttons_released_this_frame; // Mouse buttons released this specific frame

    // Current mouse position this frame
    Vector2F m_mouse_position;

    // Most recent mouse position when a mouse button was pressed.
    std::array<Vector2F, MOUSE_BUTTON_COUNT> m_mouse_button_pressed_positions;

    // Accumulated position delta of the mouse position
    Vector2F m_mouse_position_delta;

    // Accumulated scroll delta of the mouse
    float m_mouse_scroll_delta = 0.0f;

    bool m_quit_requested = false; // Flag to request a shutdown
};
} // namespace DynamoEngine