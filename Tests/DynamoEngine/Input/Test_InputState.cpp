// File: Tests/DynamoEngine/Input/Test_InputState.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/Input/InputState.hpp"
#include "TestingUtilities/VectorsNear.hpp"

using namespace DynamoEngine;

// --------- FIXTURE --------- //
// Every TEST_F gets a brand new InputState, already inside its first frame (beginInputFrame called).
class InputStateTest : public ::testing::Test
{
  protected:
    void SetUp() override { state.beginInputFrame(); }

    // --- Event builders: keep each test focused on behavior, not on filling structs --- //
    static InputEvent keyEvent(InputEventType type, Key key, bool is_repeat = false)
    {
        InputEvent event;
        event.type = type;
        event.key = key;
        event.is_repeat = is_repeat;
        return event;
    }

    static InputEvent mouseButtonEvent(InputEventType type, MouseButton button, Vector2F position)
    {
        InputEvent event;
        event.type = type;
        event.mouse_button = button;
        event.mouse_position = position;
        return event;
    }

    static InputEvent mouseMoveEvent(Vector2F position, Vector2F delta)
    {
        InputEvent event;
        event.type = InputEventType::MouseMove;
        event.mouse_position = position;
        event.mouse_position_delta = delta;
        return event;
    }

    static InputEvent wheelEvent(float scroll_delta)
    {
        InputEvent event;
        event.type = InputEventType::MouseWheel;
        event.scroll_delta = scroll_delta;
        return event;
    }

    static InputEvent eventOfType(InputEventType type)
    {
        InputEvent event;
        event.type = type;
        return event;
    }

    static constexpr float TOLERANCE = 1e-6f;

    InputState state;
};

// --------- DEFAULTS --------- //

TEST_F(InputStateTest, StartsWithNothingHeldAndNoQuit)
{
    EXPECT_FALSE(state.isKeyDown(SDL_SCANCODE_W));
    EXPECT_FALSE(state.wasKeyPressed(SDL_SCANCODE_W));
    EXPECT_FALSE(state.wasKeyReleased(SDL_SCANCODE_W));
    EXPECT_FALSE(state.isAnyMouseButtonDown());
    EXPECT_FALSE(state.wasMouseButtonPressed(MouseButton::Left));
    EXPECT_TRUE(VectorsNear(state.mousePosition(), {0.0f, 0.0f}, TOLERANCE));
    EXPECT_TRUE(VectorsNear(state.mousePositionDeltaThisFrame(), {0.0f, 0.0f}, TOLERANCE));
    EXPECT_FLOAT_EQ(state.mouseScrollDeltaThisFrame(), 0.0f);
    EXPECT_FALSE(state.quitRequested());
}

// --------- KEYBOARD --------- //

TEST_F(InputStateTest, KeyPressIsAnEdgeButKeyDownPersists)
{
    state.applyInputEvent(keyEvent(InputEventType::KeyDown, SDL_SCANCODE_W));
    EXPECT_TRUE(state.isKeyDown(SDL_SCANCODE_W));
    EXPECT_TRUE(state.wasKeyPressed(SDL_SCANCODE_W));

    state.beginInputFrame(); // next frame: key still physically held, but no new press
    EXPECT_TRUE(state.isKeyDown(SDL_SCANCODE_W));
    EXPECT_FALSE(state.wasKeyPressed(SDL_SCANCODE_W));
}

TEST_F(InputStateTest, KeyAutoRepeatIsNotANewPress)
{
    state.applyInputEvent(keyEvent(InputEventType::KeyDown, SDL_SCANCODE_ESCAPE));
    state.beginInputFrame();

    state.applyInputEvent(keyEvent(InputEventType::KeyDown, SDL_SCANCODE_ESCAPE, /*is_repeat=*/true));
    EXPECT_TRUE(state.isKeyDown(SDL_SCANCODE_ESCAPE));
    EXPECT_FALSE(state.wasKeyPressed(SDL_SCANCODE_ESCAPE)); // holding Esc must not re-toggle the pause menu
}

TEST_F(InputStateTest, KeyReleaseIsAnEdgeAndClearsDown)
{
    state.applyInputEvent(keyEvent(InputEventType::KeyDown, SDL_SCANCODE_W));
    state.beginInputFrame();

    state.applyInputEvent(keyEvent(InputEventType::KeyUp, SDL_SCANCODE_W));
    EXPECT_FALSE(state.isKeyDown(SDL_SCANCODE_W));
    EXPECT_TRUE(state.wasKeyReleased(SDL_SCANCODE_W));

    state.beginInputFrame();
    EXPECT_FALSE(state.wasKeyReleased(SDL_SCANCODE_W));
}

TEST_F(InputStateTest, KeyPressedAndReleasedInSameFrameKeepsBothEdges)
{
    state.applyInputEvent(keyEvent(InputEventType::KeyDown, SDL_SCANCODE_SPACE));
    state.applyInputEvent(keyEvent(InputEventType::KeyUp, SDL_SCANCODE_SPACE));

    EXPECT_TRUE(state.wasKeyPressed(SDL_SCANCODE_SPACE)); // a very fast tap still counts as a press
    EXPECT_TRUE(state.wasKeyReleased(SDL_SCANCODE_SPACE));
    EXPECT_FALSE(state.isKeyDown(SDL_SCANCODE_SPACE));
}

TEST_F(InputStateTest, KeysAreTrackedIndependently)
{
    state.applyInputEvent(keyEvent(InputEventType::KeyDown, SDL_SCANCODE_W));
    EXPECT_TRUE(state.isKeyDown(SDL_SCANCODE_W));
    EXPECT_FALSE(state.isKeyDown(SDL_SCANCODE_S));
}

TEST_F(InputStateTest, ModifierHelpersAcceptEitherSide)
{
    state.applyInputEvent(keyEvent(InputEventType::KeyDown, SDL_SCANCODE_RSHIFT));
    EXPECT_TRUE(state.isShiftDown()); // right Shift counts, not just left

    state.applyInputEvent(keyEvent(InputEventType::KeyDown, SDL_SCANCODE_LCTRL));
    EXPECT_TRUE(state.isCtrlDown());

    state.applyInputEvent(keyEvent(InputEventType::KeyDown, SDL_SCANCODE_RALT));
    EXPECT_TRUE(state.isAltDown());

    state.applyInputEvent(keyEvent(InputEventType::KeyDown, SDL_SCANCODE_LGUI));
    EXPECT_TRUE(state.isSuperDown());

    state.applyInputEvent(keyEvent(InputEventType::KeyUp, SDL_SCANCODE_RSHIFT));
    EXPECT_FALSE(state.isShiftDown());
}

// --------- MOUSE BUTTONS --------- //

TEST_F(InputStateTest, MouseButtonPressRecordsDownEdgeAndPosition)
{
    state.applyInputEvent(mouseButtonEvent(InputEventType::MouseButtonDown, MouseButton::Left, {100.0f, 200.0f}));

    EXPECT_TRUE(state.isMouseButtonDown(MouseButton::Left));
    EXPECT_TRUE(state.wasMouseButtonPressed(MouseButton::Left));
    EXPECT_TRUE(state.isAnyMouseButtonDown());
    EXPECT_TRUE(VectorsNear(state.pressPosition(MouseButton::Left), {100.0f, 200.0f}, TOLERANCE));
    EXPECT_TRUE(VectorsNear(state.mousePosition(), {100.0f, 200.0f}, TOLERANCE));
}

TEST_F(InputStateTest, MouseMoveUpdatesPositionButNotPressPosition)
{
    state.applyInputEvent(mouseButtonEvent(InputEventType::MouseButtonDown, MouseButton::Left, {100.0f, 200.0f}));
    state.applyInputEvent(mouseMoveEvent({150.0f, 260.0f}, {50.0f, 60.0f}));

    EXPECT_TRUE(VectorsNear(state.mousePosition(), {150.0f, 260.0f}, TOLERANCE));
    EXPECT_TRUE(
        VectorsNear(state.pressPosition(MouseButton::Left), {100.0f, 200.0f}, TOLERANCE)); // drag start is fixed
}

TEST_F(InputStateTest, MouseButtonReleaseIsAnEdgeAndClearsDown)
{
    state.applyInputEvent(mouseButtonEvent(InputEventType::MouseButtonDown, MouseButton::Right, {10.0f, 10.0f}));
    state.beginInputFrame();

    state.applyInputEvent(mouseButtonEvent(InputEventType::MouseButtonUp, MouseButton::Right, {20.0f, 30.0f}));
    EXPECT_FALSE(state.isMouseButtonDown(MouseButton::Right));
    EXPECT_TRUE(state.wasMouseButtonReleased(MouseButton::Right));
    EXPECT_FALSE(state.isAnyMouseButtonDown());
    EXPECT_TRUE(VectorsNear(state.mousePosition(), {20.0f, 30.0f}, TOLERANCE)); // release position is tracked too
}

TEST_F(InputStateTest, FastClickKeepsBothEdgesInOneFrame)
{
    state.applyInputEvent(mouseButtonEvent(InputEventType::MouseButtonDown, MouseButton::Left, {5.0f, 5.0f}));
    state.applyInputEvent(mouseButtonEvent(InputEventType::MouseButtonUp, MouseButton::Left, {5.0f, 5.0f}));

    EXPECT_TRUE(state.wasMouseButtonPressed(MouseButton::Left));
    EXPECT_TRUE(state.wasMouseButtonReleased(MouseButton::Left));
    EXPECT_FALSE(state.isMouseButtonDown(MouseButton::Left));
}

TEST_F(InputStateTest, MouseButtonsAreTrackedIndependently)
{
    state.applyInputEvent(mouseButtonEvent(InputEventType::MouseButtonDown, MouseButton::Middle, {0.0f, 0.0f}));
    EXPECT_TRUE(state.isMouseButtonDown(MouseButton::Middle));
    EXPECT_FALSE(state.isMouseButtonDown(MouseButton::Left));
    EXPECT_FALSE(state.isMouseButtonDown(MouseButton::Right));
}

// --------- PER-FRAME DELTAS --------- //

TEST_F(InputStateTest, MouseDeltasAccumulateWithinAFrameAndResetNextFrame)
{
    state.applyInputEvent(mouseMoveEvent({10.0f, 0.0f}, {10.0f, 0.0f}));
    state.applyInputEvent(mouseMoveEvent({15.0f, 5.0f}, {5.0f, 5.0f}));
    EXPECT_TRUE(VectorsNear(state.mousePositionDeltaThisFrame(), {15.0f, 5.0f}, TOLERANCE));

    state.beginInputFrame();
    EXPECT_TRUE(VectorsNear(state.mousePositionDeltaThisFrame(), {0.0f, 0.0f}, TOLERANCE));
    EXPECT_TRUE(VectorsNear(state.mousePosition(), {15.0f, 5.0f}, TOLERANCE)); // position itself persists
}

TEST_F(InputStateTest, ScrollAccumulatesWithinAFrameAndResetsNextFrame)
{
    state.applyInputEvent(wheelEvent(1.0f));
    state.applyInputEvent(wheelEvent(2.0f));
    EXPECT_FLOAT_EQ(state.mouseScrollDeltaThisFrame(), 3.0f);

    state.beginInputFrame();
    EXPECT_FLOAT_EQ(state.mouseScrollDeltaThisFrame(), 0.0f);
}

// --------- WINDOW / APPLICATION --------- //

TEST_F(InputStateTest, FocusLostReleasesEverythingHeldWithoutQuitting)
{
    state.applyInputEvent(keyEvent(InputEventType::KeyDown, SDL_SCANCODE_W));
    state.applyInputEvent(mouseButtonEvent(InputEventType::MouseButtonDown, MouseButton::Left, {0.0f, 0.0f}));

    state.applyInputEvent(eventOfType(InputEventType::FocusLost)); // e.g. Alt-Tab / Cmd-Tab away

    EXPECT_FALSE(state.isKeyDown(SDL_SCANCODE_W));
    EXPECT_FALSE(state.isAnyMouseButtonDown());
    EXPECT_FALSE(state.quitRequested()); // regression: FocusLost must not fall through into Quit
}

TEST_F(InputStateTest, QuitRequestPersistsAcrossFrames)
{
    state.applyInputEvent(eventOfType(InputEventType::Quit));
    EXPECT_TRUE(state.quitRequested());

    state.beginInputFrame(); // quitting is not a per-frame edge
    EXPECT_TRUE(state.quitRequested());
}

TEST_F(InputStateTest, EventsWithoutDeviceStateChangeNothing)
{
    state.applyInputEvent(eventOfType(InputEventType::None));
    state.applyInputEvent(eventOfType(InputEventType::TextInput));
    state.applyInputEvent(eventOfType(InputEventType::WindowResize));

    EXPECT_FALSE(state.isAnyMouseButtonDown());
    EXPECT_FALSE(state.isKeyDown(SDL_SCANCODE_UNKNOWN));
    EXPECT_FALSE(state.quitRequested());
}
