// File: Tests/DynamoEngine/Input/Test_InputEvent.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/Constants/GlobalConstants.hpp"
#include "DynamoEngine/Input/InputEvent.hpp"
#include "TestingUtilities/VectorsNear.hpp"

// Standard Library Imports

TEST(InputEvent, DefaultsDoNothing)
{
    DynamoEngine::InputEvent new_input_event;
    EXPECT_EQ(new_input_event.key, SDL_SCANCODE_UNKNOWN);
    EXPECT_EQ(new_input_event.is_repeat, false);
    EXPECT_EQ(new_input_event.key_modifiers, SDL_KMOD_NONE);
    EXPECT_EQ(new_input_event.mouse_button, DynamoEngine::MouseButton::None);
    EXPECT_TRUE(VectorsNear(new_input_event.mouse_position, DynamoEngine::Vector2F(0.0f, 0.0f), DynamoEngine::EPSILON));
    EXPECT_TRUE(
        VectorsNear(new_input_event.mouse_position_delta, DynamoEngine::Vector2F(0.0f, 0.0f), DynamoEngine::EPSILON));
    EXPECT_EQ(new_input_event.scrollDelta, 0);
    EXPECT_EQ(new_input_event.text_input, "");
    EXPECT_EQ(new_input_event.type, DynamoEngine::InputEventType::None);
    EXPECT_EQ(new_input_event.window_height, 0);
    EXPECT_EQ(new_input_event.window_width, 0);
}
