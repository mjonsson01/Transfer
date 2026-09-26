// File: Tests/DynamoEngine/UI/Test_UIRootInput.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/Input/InputEvent.hpp"
#include "DynamoEngine/Input/InputState.hpp"
#include "DynamoEngine/UI/UIRoot.hpp"

// Standard Library Imports
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace DynamoEngine;

// --------- TEST HELPERS --------- //

// An element that records every mouse event it receives, and can be told whether to take presses
class RecordingElement : public UIElement
{
  public:
    RecordingElement(std::string name, std::vector<std::string>& log, bool takes_presses)
        : m_name(std::move(name)), m_log(log), m_takes_presses(takes_presses)
    {
    }

    bool onMousePressed(Vector2F mouse_position) override
    {
        m_log.push_back(m_name + " pressed");
        return m_takes_presses;
    }
    void onMouseDragged(Vector2F mouse_position) override { m_log.push_back(m_name + " dragged"); }
    void onMouseReleased(Vector2F mouse_position, bool released_inside) override
    {
        m_log.push_back(m_name + (released_inside ? " released inside" : " released outside"));
    }
    void onMouseEntered() override { m_log.push_back(m_name + " entered"); }
    void onMouseExited() override { m_log.push_back(m_name + " exited"); }

  private:
    std::string m_name;
    std::vector<std::string>& m_log;
    bool m_takes_presses;
};

// The window is the default 1280 x 720 (UI scale 1), so mouse positions and element rects use the same numbers.
// Elements are 100 x 100: TopLeft covers x 0-100, y 0-100; TopRight covers x 1180-1280, y 0-100.
class UIRootInputTest : public ::testing::Test
{
  protected:
    RecordingElement& addElement(const std::string& name, UIAlign align, bool takes_presses = true)
    {
        auto element = std::make_unique<RecordingElement>(name, log, takes_presses);
        element->setPlacement({.align = align, .size = {100.0f, 100.0f}, .margin = 0.0f});
        RecordingElement& added = *element;
        root.addChild(std::move(element));
        return added;
    }

    // One frame: apply this frame's input events, lay out the UI, then let the UI process the input
    UIInputResult frame(const std::vector<InputEvent>& events)
    {
        input.beginInputFrame();
        for (const InputEvent& event : events)
        {
            input.applyInputEvent(event);
        }
        root.updateElements(0.016f);
        return root.processInput(input);
    }

    static InputEvent moveTo(Vector2F position)
    {
        InputEvent event;
        event.type = InputEventType::MouseMove;
        event.mouse_position = position;
        return event;
    }
    static InputEvent press(Vector2F position)
    {
        InputEvent event;
        event.type = InputEventType::MouseButtonDown;
        event.mouse_button = MouseButton::Left;
        event.mouse_position = position;
        return event;
    }
    static InputEvent release(Vector2F position)
    {
        InputEvent event;
        event.type = InputEventType::MouseButtonUp;
        event.mouse_button = MouseButton::Left;
        event.mouse_position = position;
        return event;
    }
    static InputEvent focusLost()
    {
        InputEvent event;
        event.type = InputEventType::FocusLost;
        return event;
    }

    UIRoot root;
    InputState input;
    std::vector<std::string> log;
};

// --------- HOVER --------- //

TEST_F(UIRootInputTest, HoverEntersAndExitsAsTheMouseMoves)
{
    addElement("button", UIAlign::TopLeft);

    frame({moveTo({50.0f, 50.0f})});
    frame({moveTo({60.0f, 60.0f})}); // still over it: no repeat
    frame({moveTo({500.0f, 500.0f})});

    EXPECT_EQ(log, (std::vector<std::string>{"button entered", "button exited"}));
}

TEST_F(UIRootInputTest, HoverMovesFromOneElementToAnother)
{
    addElement("left", UIAlign::TopLeft);
    addElement("right", UIAlign::TopRight);

    frame({moveTo({50.0f, 50.0f})});
    frame({moveTo({1230.0f, 50.0f})});

    EXPECT_EQ(log, (std::vector<std::string>{"left entered", "left exited", "right entered"}));
}

TEST_F(UIRootInputTest, PointerOverUIIsReported)
{
    addElement("button", UIAlign::TopLeft);

    EXPECT_TRUE(frame({moveTo({50.0f, 50.0f})}).pointer_over_ui);
    EXPECT_FALSE(frame({moveTo({500.0f, 500.0f})}).pointer_over_ui);
}

// --------- PRESS / DRAG / RELEASE --------- //

TEST_F(UIRootInputTest, ClickGoesPressDragReleaseInside)
{
    addElement("button", UIAlign::TopLeft);

    frame({moveTo({50.0f, 50.0f})});
    log.clear();
    frame({press({50.0f, 50.0f})});
    frame({moveTo({55.0f, 55.0f})});
    frame({release({55.0f, 55.0f})});

    EXPECT_EQ(log, (std::vector<std::string>{"button pressed", "button dragged", "button dragged",
                                             "button released inside"}));
}

TEST_F(UIRootInputTest, CapturedElementKeepsTheDragAfterTheMouseLeavesIt)
{
    addElement("slider", UIAlign::TopLeft);

    frame({press({50.0f, 50.0f})});
    log.clear();
    frame({moveTo({900.0f, 50.0f})}); // dragged far off the slider
    frame({release({900.0f, 50.0f})});

    EXPECT_EQ(log, (std::vector<std::string>{"slider exited", "slider dragged", "slider released outside"}));
}

TEST_F(UIRootInputTest, ReleaseFrameStillCountsAsCapturedSoTheGameIgnoresIt)
{
    addElement("slider", UIAlign::TopLeft);

    EXPECT_TRUE(frame({press({50.0f, 50.0f})}).pointer_captured);
    // Released over empty space: without this, the game would treat it as a spawn click
    EXPECT_TRUE(frame({release({900.0f, 500.0f})}).pointer_captured);
    EXPECT_FALSE(frame({moveTo({900.0f, 500.0f})}).pointer_captured); // next frame: free again
}

TEST_F(UIRootInputTest, FastClickPressAndReleaseInOneFrame)
{
    addElement("button", UIAlign::TopLeft);

    UIInputResult result = frame({press({50.0f, 50.0f}), release({50.0f, 50.0f})});

    EXPECT_TRUE(result.pointer_captured);
    EXPECT_EQ(log, (std::vector<std::string>{"button entered", "button pressed", "button released inside"}));
}

TEST_F(UIRootInputTest, PressPassesThroughAnElementThatDeclinesIt)
{
    addElement("button", UIAlign::TopLeft);
    addElement("label_on_top", UIAlign::TopLeft, /*takes_presses=*/false); // drawn above the button

    frame({press({50.0f, 50.0f})});

    EXPECT_EQ(log, (std::vector<std::string>{"label_on_top entered", "label_on_top pressed", "button pressed",
                                             "button dragged"}));
}

TEST_F(UIRootInputTest, PressOnEmptySpaceIsLeftForTheGame)
{
    addElement("button", UIAlign::TopLeft);

    UIInputResult result = frame({press({900.0f, 500.0f})});

    EXPECT_FALSE(result.pointer_captured);
    EXPECT_FALSE(result.pointer_over_ui);
    EXPECT_TRUE(log.empty());
}

TEST_F(UIRootInputTest, LosingWindowFocusReleasesTheCapture)
{
    addElement("slider", UIAlign::TopLeft);

    frame({press({50.0f, 50.0f})});
    log.clear();
    frame({focusLost()}); // e.g. Cmd-Tab mid-drag: no MouseButtonUp ever arrives

    EXPECT_EQ(log, (std::vector<std::string>{"slider released inside"}));
    EXPECT_FALSE(frame({}).pointer_captured);
}

TEST_F(UIRootInputTest, KeyboardIsNeverCapturedYet)
{
    addElement("button", UIAlign::TopLeft);
    EXPECT_FALSE(frame({press({50.0f, 50.0f})}).keyboard_captured);
}
