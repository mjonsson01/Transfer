// File: Transfer/src/DynamoEngine/UI/UIInputResult.hpp

#pragma once

namespace DynamoEngine
{
// What the UI took from this frame's input. The game checks this before acting on the mouse or keyboard,
// so a click on a button never also spawns a planet behind it.
struct UIInputResult
{
    bool pointer_captured =
        false; // the UI owns the current press (including the frame it's released): game ignores mouse buttons
    bool pointer_over_ui = false; // the cursor is over a UI element: game skips hover-picking and wheel zoom
    bool keyboard_captured =
        false; // a focused element owns key presses: game ignores keys (used once keyboard focus exists)
};
} // namespace DynamoEngine