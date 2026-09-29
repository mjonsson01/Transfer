// File: Tests/DynamoEngine/Input/Test_SDLInputIntake.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// SDL Imports
#include <SDL3/SDL_events.h>
#include <SDL3/SDL_init.h>
#include <SDL3/SDL_keycode.h>
#include <SDL3/SDL_mouse.h>

// Custom Imports
#include "DynamoEngine/Input/SDLInputIntake.hpp"
#include "TestingUtilities/VectorsNear.hpp"

// Standard Library Imports
#include <vector>

// Brings the engine's input types into scope for this test file only (fine in a .cpp, never in a header)
using namespace DynamoEngine;

// --------- FIXTURE --------- //
// translateSDLEvent is a pure function, so most tests just hand-build an SDL_Event and check what comes out.
// No SDL_Init needed for those -- only the pollEvents tests below touch SDL's real event queue.
class SDLInputIntakeTest : public ::testing::Test
{
  protected:
    // --- SDL event builders --- //
    // SDL_Event is a union: set `type` first, then only the member struct that matches it.
    static SDL_Event sdlEventOfType(SDL_EventType type)
    {
        SDL_Event sdl_event{};
        sdl_event.type = type;
        return sdl_event;
    }

    static SDL_Event sdlKeyEvent(SDL_EventType type, SDL_Scancode scancode, bool is_repeat, SDL_Keymod modifiers)
    {
        SDL_Event sdl_event = sdlEventOfType(type);
        sdl_event.key.scancode = scancode;
        sdl_event.key.repeat = is_repeat;
        sdl_event.key.mod = modifiers;
        return sdl_event;
    }

    static SDL_Event sdlMouseButtonEvent(SDL_EventType type, Uint8 button, float x, float y)
    {
        SDL_Event sdl_event = sdlEventOfType(type);
        sdl_event.button.button = button;
        sdl_event.button.x = x;
        sdl_event.button.y = y;
        return sdl_event;
    }

    // An InputEvent pre-filled with recognizable values, to prove a `false` translation leaves it untouched
    static InputEvent sentinelEvent()
    {
        InputEvent event;
        event.type = InputEventType::KeyDown;
        event.key = SDL_SCANCODE_Z;
        return event;
    }

    static void expectUntouched(const InputEvent& event)
    {
        EXPECT_EQ(event.type, InputEventType::KeyDown);
        EXPECT_EQ(event.key, SDL_SCANCODE_Z);
    }

    static constexpr float TOLERANCE = 1e-6f;
};

// --------- APPLICATION / WINDOW --------- //

TEST_F(SDLInputIntakeTest, QuitTranslatesToQuit)
{
    InputEvent event;
    ASSERT_TRUE(translateSDLEvent(sdlEventOfType(SDL_EVENT_QUIT), event));
    EXPECT_EQ(event.type, InputEventType::Quit);
}

TEST_F(SDLInputIntakeTest, FocusLostTranslatesToFocusLost)
{
    InputEvent event;
    ASSERT_TRUE(translateSDLEvent(sdlEventOfType(SDL_EVENT_WINDOW_FOCUS_LOST), event));
    EXPECT_EQ(event.type, InputEventType::FocusLost);
}

TEST_F(SDLInputIntakeTest, WindowResizeCarriesNewSize)
{
    SDL_Event sdl_event = sdlEventOfType(SDL_EVENT_WINDOW_RESIZED);
    sdl_event.window.data1 = 1920;
    sdl_event.window.data2 = 1080;

    InputEvent event;
    ASSERT_TRUE(translateSDLEvent(sdl_event, event));
    EXPECT_EQ(event.type, InputEventType::WindowResize);
    EXPECT_EQ(event.window_width, 1920);
    EXPECT_EQ(event.window_height, 1080);
}

// --------- KEYBOARD --------- //

TEST_F(SDLInputIntakeTest, KeyDownCarriesScancodeRepeatAndModifiers)
{
    InputEvent event;
    ASSERT_TRUE(translateSDLEvent(sdlKeyEvent(SDL_EVENT_KEY_DOWN, SDL_SCANCODE_W, true, SDL_KMOD_LSHIFT), event));
    EXPECT_EQ(event.type, InputEventType::KeyDown);
    EXPECT_EQ(event.key, SDL_SCANCODE_W);
    EXPECT_TRUE(event.is_repeat);
    EXPECT_TRUE(event.key_modifiers & SDL_KMOD_SHIFT);
}

TEST_F(SDLInputIntakeTest, KeyUpTranslatesToKeyUp)
{
    InputEvent event;
    ASSERT_TRUE(translateSDLEvent(sdlKeyEvent(SDL_EVENT_KEY_UP, SDL_SCANCODE_ESCAPE, false, SDL_KMOD_NONE), event));
    EXPECT_EQ(event.type, InputEventType::KeyUp);
    EXPECT_EQ(event.key, SDL_SCANCODE_ESCAPE);
    EXPECT_FALSE(event.is_repeat);
    EXPECT_EQ(event.key_modifiers, SDL_KMOD_NONE);
}

// --------- MOUSE BUTTONS --------- //

TEST_F(SDLInputIntakeTest, EachTrackedMouseButtonMapsToTheRightEngineButton)
{
    struct Case
    {
        Uint8 sdl_button;
        MouseButton expected;
    };
    const Case cases[] = {{SDL_BUTTON_LEFT, MouseButton::Left},
                          {SDL_BUTTON_RIGHT, MouseButton::Right},
                          {SDL_BUTTON_MIDDLE, MouseButton::Middle}};

    for (const Case& test_case : cases)
    {
        SCOPED_TRACE(static_cast<int>(test_case.sdl_button)); // names the failing button if an EXPECT below fails
        InputEvent event;
        ASSERT_TRUE(translateSDLEvent(sdlMouseButtonEvent(SDL_EVENT_MOUSE_BUTTON_DOWN, test_case.sdl_button, 0, 0), event));
        EXPECT_EQ(event.mouse_button, test_case.expected);
    }
}

TEST_F(SDLInputIntakeTest, MouseButtonDownCarriesTypeAndPosition)
{
    InputEvent event;
    ASSERT_TRUE(translateSDLEvent(sdlMouseButtonEvent(SDL_EVENT_MOUSE_BUTTON_DOWN, SDL_BUTTON_LEFT, 12.5f, 40.0f), event));
    EXPECT_EQ(event.type, InputEventType::MouseButtonDown);
    EXPECT_TRUE(VectorsNear(event.mouse_position, {12.5f, 40.0f}, TOLERANCE));
}

TEST_F(SDLInputIntakeTest, MouseButtonUpCarriesTypeAndPosition)
{
    InputEvent event;
    ASSERT_TRUE(translateSDLEvent(sdlMouseButtonEvent(SDL_EVENT_MOUSE_BUTTON_UP, SDL_BUTTON_RIGHT, 300.0f, 5.0f), event));
    EXPECT_EQ(event.type, InputEventType::MouseButtonUp);
    EXPECT_EQ(event.mouse_button, MouseButton::Right);
    EXPECT_TRUE(VectorsNear(event.mouse_position, {300.0f, 5.0f}, TOLERANCE));
}

TEST_F(SDLInputIntakeTest, UntrackedMouseButtonIsRejectedAndLeavesOutputUntouched)
{
    InputEvent event = sentinelEvent();
    EXPECT_FALSE(translateSDLEvent(sdlMouseButtonEvent(SDL_EVENT_MOUSE_BUTTON_DOWN, SDL_BUTTON_X1, 0, 0), event));
    expectUntouched(event);
}

// --------- MOUSE MOTION / WHEEL --------- //

TEST_F(SDLInputIntakeTest, MouseMotionCarriesPositionAndDelta)
{
    SDL_Event sdl_event = sdlEventOfType(SDL_EVENT_MOUSE_MOTION);
    sdl_event.motion.x = 100.0f;
    sdl_event.motion.y = 50.0f;
    sdl_event.motion.xrel = -3.0f;
    sdl_event.motion.yrel = 4.0f;

    InputEvent event;
    ASSERT_TRUE(translateSDLEvent(sdl_event, event));
    EXPECT_EQ(event.type, InputEventType::MouseMove);
    EXPECT_TRUE(VectorsNear(event.mouse_position, {100.0f, 50.0f}, TOLERANCE));
    EXPECT_TRUE(VectorsNear(event.mouse_position_delta, {-3.0f, 4.0f}, TOLERANCE));
}

TEST_F(SDLInputIntakeTest, MouseWheelCarriesScrollAndCursorPosition)
{
    SDL_Event sdl_event = sdlEventOfType(SDL_EVENT_MOUSE_WHEEL);
    sdl_event.wheel.y = 2.0f;
    sdl_event.wheel.direction = SDL_MOUSEWHEEL_NORMAL;
    sdl_event.wheel.mouse_x = 640.0f;
    sdl_event.wheel.mouse_y = 360.0f;

    InputEvent event;
    ASSERT_TRUE(translateSDLEvent(sdl_event, event));
    EXPECT_EQ(event.type, InputEventType::MouseWheel);
    EXPECT_FLOAT_EQ(event.scroll_delta, 2.0f);
    EXPECT_TRUE(VectorsNear(event.mouse_position, {640.0f, 360.0f}, TOLERANCE));
}

TEST_F(SDLInputIntakeTest, FlippedWheelIsNormalizedToTheSameDirection)
{
    SDL_Event sdl_event = sdlEventOfType(SDL_EVENT_MOUSE_WHEEL);
    sdl_event.wheel.y = 2.0f;
    sdl_event.wheel.direction = SDL_MOUSEWHEEL_FLIPPED; // "natural scrolling" on macOS

    InputEvent event;
    ASSERT_TRUE(translateSDLEvent(sdl_event, event));
    EXPECT_FLOAT_EQ(event.scroll_delta, -2.0f);
}

// --------- TEXT --------- //

TEST_F(SDLInputIntakeTest, TextInputCarriesTheTypedText)
{
    SDL_Event sdl_event = sdlEventOfType(SDL_EVENT_TEXT_INPUT);
    sdl_event.text.text = "é!"; // UTF-8: one character can be several bytes

    InputEvent event;
    ASSERT_TRUE(translateSDLEvent(sdl_event, event));
    EXPECT_EQ(event.type, InputEventType::TextInput);
    EXPECT_EQ(event.text_input, "é!");
}

// --------- NON-INPUT EVENTS --------- //

TEST_F(SDLInputIntakeTest, NonInputEventIsRejectedAndLeavesOutputUntouched)
{
    InputEvent event = sentinelEvent();
    EXPECT_FALSE(translateSDLEvent(sdlEventOfType(SDL_EVENT_GAMEPAD_ADDED), event));
    expectUntouched(event);
}

// --------- pollEvents (uses SDL's real event queue) --------- //
// These need SDL's event subsystem running. SetUpTestSuite/TearDownTestSuite run once for the whole suite.
class SDLInputIntakePollTest : public ::testing::Test
{
  protected:
    static void SetUpTestSuite() { ASSERT_TRUE(SDL_Init(SDL_INIT_EVENTS)) << SDL_GetError(); }
    static void TearDownTestSuite() { SDL_Quit(); }

    void SetUp() override { SDL_FlushEvents(SDL_EVENT_FIRST, SDL_EVENT_LAST); } // start every test with an empty queue

    SDLInputIntake intake;
};

TEST_F(SDLInputIntakePollTest, PollTranslatesQueuedEventsInOrderAndSkipsNonInput)
{
    SDL_Event key_down{};
    key_down.type = SDL_EVENT_KEY_DOWN;
    key_down.key.scancode = SDL_SCANCODE_A;
    SDL_Event gamepad{};
    gamepad.type = SDL_EVENT_GAMEPAD_ADDED; // should be skipped
    SDL_Event quit{};
    quit.type = SDL_EVENT_QUIT;

    ASSERT_TRUE(SDL_PushEvent(&key_down));
    ASSERT_TRUE(SDL_PushEvent(&gamepad));
    ASSERT_TRUE(SDL_PushEvent(&quit));

    std::vector<InputEvent> events;
    intake.pollEvents(events);

    ASSERT_EQ(events.size(), 2u);
    EXPECT_EQ(events[0].type, InputEventType::KeyDown);
    EXPECT_EQ(events[0].key, SDL_SCANCODE_A);
    EXPECT_EQ(events[1].type, InputEventType::Quit);
}

TEST_F(SDLInputIntakePollTest, PollAppendsWithoutClearing)
{
    std::vector<InputEvent> events(1); // one pre-existing event

    SDL_Event quit{};
    quit.type = SDL_EVENT_QUIT;
    ASSERT_TRUE(SDL_PushEvent(&quit));

    intake.pollEvents(events);
    ASSERT_EQ(events.size(), 2u);
    EXPECT_EQ(events[1].type, InputEventType::Quit);
}
