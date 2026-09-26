// File: Tests/DynamoEngine/UI/Widgets/Test_UIWidgets.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/Rendering/FontAtlas.hpp"
#include "DynamoEngine/UI/UIGeometryBuilder.hpp"
#include "DynamoEngine/UI/UIRoot.hpp"
#include "DynamoEngine/UI/Widgets/UIButton.hpp"
#include "DynamoEngine/UI/Widgets/UILabel.hpp"
#include "DynamoEngine/UI/Widgets/UIRow.hpp"
#include "DynamoEngine/UI/Widgets/UISlider.hpp"

// Standard Library Imports
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace DynamoEngine;

// --------- HELPERS --------- //

// The first vertex's red channel: enough to tell which background color a button or knob drew with
float firstVertexRed(const UIElement& element)
{
    std::vector<UIVertex> vertices;
    FontAtlas empty_atlas;
    UIGeometryBuilder builder(vertices, empty_atlas);
    element.draw(builder);
    return vertices.empty() ? -1.0f : vertices[0].r;
}

// A slider draws its track first (vertices 0-5), then its knob (vertices 6-11): this is the knob's red channel
float knobRed(const UISlider& slider)
{
    std::vector<UIVertex> vertices;
    FontAtlas empty_atlas;
    UIGeometryBuilder builder(vertices, empty_atlas);
    slider.draw(builder);
    return vertices.size() < 12 ? -1.0f : vertices[6].r;
}

// --------- SOUND --------- //

TEST(UISound, RequestsTravelUpToTheRootHandler)
{
    UIRoot root;
    std::vector<UISound> sounds_heard;
    root.setSoundHandler([&sounds_heard](UISound sound) { sounds_heard.push_back(sound); });

    UIElement& panel = root.addChild(std::make_unique<UIElement>());
    UIElement& deep_child = panel.addChild(std::make_unique<UIElement>());
    deep_child.requestSound(UISound::Click);

    EXPECT_EQ(sounds_heard, (std::vector<UISound>{UISound::Click}));
}

TEST(UISound, NoHandlerMeansSilence)
{
    UIRoot root; // no handler set
    UIElement& child = root.addChild(std::make_unique<UIElement>());
    child.requestSound(UISound::Click); // must simply do nothing
    SUCCEED();
}

// --------- LABEL --------- //

TEST(UILabel, TextSourceRefreshesTheTextEveryUpdate)
{
    int frames = 0;
    UILabel label("start");
    label.setTextSource([&frames]() { return "frame " + std::to_string(frames); });

    frames = 7;
    label.update(0.016f);
    EXPECT_EQ(label.text(), "frame 7");
}

// --------- BUTTON --------- //

class UIButtonTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        root.setSoundHandler([this](UISound sound) { sounds_heard.push_back(sound); });
        auto new_button = std::make_unique<UIButton>("Resume");
        new_button->setOnClick([this]() { clicks += 1; });
        button = static_cast<UIButton*>(&root.addChild(std::move(new_button)));
    }

    UIRoot root;
    UIButton* button = nullptr;
    int clicks = 0;
    std::vector<UISound> sounds_heard;
};

TEST_F(UIButtonTest, ReleaseInsideClicksAndPlaysTheClickSound)
{
    EXPECT_TRUE(button->onMousePressed({0.0f, 0.0f})); // buttons always take the press
    button->onMouseReleased({0.0f, 0.0f}, /*released_inside=*/true);

    EXPECT_EQ(clicks, 1);
    EXPECT_EQ(sounds_heard, (std::vector<UISound>{UISound::Click}));
}

TEST_F(UIButtonTest, ReleaseOutsideCancelsTheClick)
{
    button->onMousePressed({0.0f, 0.0f});
    button->onMouseReleased({900.0f, 900.0f}, /*released_inside=*/false);

    EXPECT_EQ(clicks, 0);
    EXPECT_TRUE(sounds_heard.empty());
}

TEST_F(UIButtonTest, DarkensWhenHoveredAndMoreWhenHeld)
{
    const float normal = firstVertexRed(*button);
    button->onMouseEntered();
    const float hovered = firstVertexRed(*button);
    button->onMousePressed({0.0f, 0.0f});
    const float held = firstVertexRed(*button);

    EXPECT_LT(hovered, normal);
    EXPECT_LT(held, hovered);
}

TEST(UIButton, ButtonWithoutAnActionIsSafeToClick)
{
    UIButton button("Nothing");
    button.onMousePressed({0.0f, 0.0f});
    button.onMouseReleased({0.0f, 0.0f}, true);
    SUCCEED();
}

// --------- SLIDER --------- //

// A 220-wide slider at x 0: the knob is 20 wide, so it travels 200 points. Mouse at x = 10 + 200 * p puts
// the knob's center at position p.
class UISliderTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        root.setSoundHandler([this](UISound sound) { sounds_heard.push_back(sound); });
        auto new_slider = std::make_unique<UISlider>("Speed", SliderMapping::linear(0.0, 2.0), 1.0);
        new_slider->setPlacement({.align = UIAlign::TopLeft, .size = {220.0f, 60.0f}, .margin = 0.0f});
        new_slider->setOnValueChanged([this](double value) { values_reported.push_back(value); });
        slider = static_cast<UISlider*>(&root.addChild(std::move(new_slider)));
        root.updateElements(0.016f);
    }

    static float mouseXForPosition(float position) { return 10.0f + 200.0f * position; }

    UIRoot root;
    UISlider* slider = nullptr;
    std::vector<double> values_reported;
    std::vector<UISound> sounds_heard;
};

TEST_F(UISliderTest, StartsAtItsStartingValue)
{
    EXPECT_DOUBLE_EQ(slider->value(), 1.0);
}

TEST_F(UISliderTest, PressingJumpsTheKnobAndReportsTheValue)
{
    EXPECT_TRUE(slider->onMousePressed({mouseXForPosition(0.25f), 10.0f}));
    EXPECT_NEAR(slider->value(), 0.5, 1e-6); // 25% of the way from 0 to 2
    ASSERT_EQ(values_reported.size(), 1u);
    EXPECT_NEAR(values_reported[0], 0.5, 1e-6);
}

TEST_F(UISliderTest, DraggingPastTheEndsClamps)
{
    slider->onMouseDragged({-500.0f, 10.0f});
    EXPECT_DOUBLE_EQ(slider->value(), 0.0);
    slider->onMouseDragged({5000.0f, 10.0f});
    EXPECT_DOUBLE_EQ(slider->value(), 2.0);
}

TEST_F(UISliderTest, TickSoundPlaysWhenCrossingTickMarksOnly)
{
    slider->onMouseDragged({mouseXForPosition(0.5f), 10.0f}); // same spot as the starting value: no tick
    EXPECT_TRUE(sounds_heard.empty());

    slider->onMouseDragged({mouseXForPosition(0.9f), 10.0f}); // moved to a different tick mark
    EXPECT_EQ(sounds_heard, (std::vector<UISound>{UISound::Tick}));
}

TEST_F(UISliderTest, SetValueMovesTheKnobWithoutReporting)
{
    slider->setValue(2.0);
    EXPECT_DOUBLE_EQ(slider->value(), 2.0);
    EXPECT_TRUE(values_reported.empty());
}

TEST_F(UISliderTest, KnobDarkensOnlyWhileTheCursorIsOnTheKnob)
{
    // Starting value 1.0 = position 0.5: the knob covers x 100-120, y 0-30
    const float normal = knobRed(*slider);

    slider->onMouseEntered();
    slider->onMouseHover({200.0f, 10.0f}); // on the slider, but away from the knob
    EXPECT_FLOAT_EQ(knobRed(*slider), normal);

    slider->onMouseHover({110.0f, 10.0f}); // on the knob
    EXPECT_LT(knobRed(*slider), normal);

    slider->onMouseHover({110.0f, 45.0f}); // below the knob, on the label
    EXPECT_FLOAT_EQ(knobRed(*slider), normal);

    slider->onMouseHover({110.0f, 10.0f});
    slider->onMouseExited(); // left the slider entirely
    EXPECT_FLOAT_EQ(knobRed(*slider), normal);
}

TEST_F(UISliderTest, KnobStaysDarkWhileDraggedEvenOffTheSlider)
{
    const float normal = knobRed(*slider);

    slider->onMousePressed({mouseXForPosition(0.2f), 10.0f});
    slider->onMouseExited(); // dragged off the slider
    EXPECT_LT(knobRed(*slider), normal);

    slider->onMouseReleased({900.0f, 900.0f}, /*released_inside=*/false);
    EXPECT_FLOAT_EQ(knobRed(*slider), normal);
}

TEST(SliderMapping, LinearGoesBothWays)
{
    SliderMapping mapping = SliderMapping::linear(-10.0, 30.0);
    EXPECT_DOUBLE_EQ(mapping.to_value(0.0), -10.0);
    EXPECT_DOUBLE_EQ(mapping.to_value(1.0), 30.0);
    EXPECT_DOUBLE_EQ(mapping.to_position(10.0), 0.5);
}

// --------- ROW / COLUMN --------- //

std::unique_ptr<UIElement> boxOfSize(Vector2F size)
{
    auto box = std::make_unique<UIElement>();
    box->setPlacement({.align = UIAlign::TopLeft, .size = size, .margin = 0.0f});
    return box;
}

TEST(UIRow, SizesToItsChildrenAndLinesThemUp)
{
    UIRoot root; // 1280 x 720
    auto row = std::make_unique<UIRow>(/*spacing=*/10.0f);
    row->setPlacement({.align = UIAlign::BottomCenter, .margin = 40.0f});
    UIElement& first = row->addChild(boxOfSize({100.0f, 30.0f}));
    UIElement& second = row->addChild(boxOfSize({50.0f, 60.0f}));
    UIElement& added_row = root.addChild(std::move(row));

    root.updateElements(0.016f);

    // Row: 100 + 10 + 50 = 160 wide, 60 tall, bottom-center 40 up -> x = 640 - 80 = 560, y = 720 - 60 - 40 = 620
    EXPECT_FLOAT_EQ(added_row.rect().x, 560.0f);
    EXPECT_FLOAT_EQ(added_row.rect().y, 620.0f);
    EXPECT_FLOAT_EQ(added_row.rect().w, 160.0f);
    EXPECT_FLOAT_EQ(added_row.rect().h, 60.0f);
    EXPECT_FLOAT_EQ(first.rect().x, 560.0f);
    EXPECT_FLOAT_EQ(second.rect().x, 670.0f); // 560 + 100 + 10
}

TEST(UIColumn, StacksChildrenTopToBottom)
{
    UIRoot root;
    auto column = std::make_unique<UIColumn>(/*spacing=*/5.0f);
    column->setPlacement({.align = UIAlign::TopLeft, .margin = 0.0f});
    UIElement& first = column->addChild(boxOfSize({80.0f, 20.0f}));
    UIElement& second = column->addChild(boxOfSize({120.0f, 20.0f}));
    UIElement& added_column = root.addChild(std::move(column));

    root.updateElements(0.016f);

    EXPECT_FLOAT_EQ(added_column.rect().w, 120.0f); // widest child
    EXPECT_FLOAT_EQ(added_column.rect().h, 45.0f);  // 20 + 5 + 20
    EXPECT_FLOAT_EQ(first.rect().y, 0.0f);
    EXPECT_FLOAT_EQ(second.rect().y, 25.0f);
}
