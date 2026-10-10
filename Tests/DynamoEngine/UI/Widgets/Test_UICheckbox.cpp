// File: Tests/DynamoEngine/UI/Widgets/Test_UICheckbox.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/Rendering/FontAtlas.hpp"
#include "DynamoEngine/UI/UIGeometryBuilder.hpp"
#include "DynamoEngine/UI/UIRoot.hpp"
#include "DynamoEngine/UI/Widgets/UICheckbox.hpp"

// Standard Library Imports
#include <memory>
#include <utility>
#include <vector>

using namespace DynamoEngine;

// --------- TEST HELPERS --------- //

class UICheckboxTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        root.setSoundHandler([this](UISound sound) { sounds_heard.push_back(sound); });
        auto new_checkbox = std::make_unique<UICheckbox>("Shatterable", false);
        new_checkbox->setPlacement({.align = UIAlign::TopLeft, .size = {200.0f, 40.0f}});
        new_checkbox->setOnToggled([this](bool is_checked) { reported.push_back(is_checked); });
        checkbox = static_cast<UICheckbox*>(&root.addChild(std::move(new_checkbox)));
        root.updateElements(0.016f); // lay it out, so m_rect is real
    }

    // Press and release on the checkbox, as a player's click would
    void clickInside()
    {
        checkbox->onMousePressed(ON_CHECKBOX);
        checkbox->onMouseReleased(ON_CHECKBOX, /*released_inside=*/true);
    }

    // How many vertices the checkbox draws (every rectangle is 2 triangles = 6 vertices)
    static size_t vertexCount(const UIElement& element)
    {
        std::vector<UIVertex> vertices;
        FontAtlas empty_atlas; // no glyphs: the label adds no vertices, so only rectangles are counted
        UIGeometryBuilder builder(vertices, empty_atlas);
        element.draw(builder);
        return vertices.size();
    }

    // The first vertex's red channel: the background's color, since the background is drawn first
    static float firstVertexRed(const UIElement& element)
    {
        std::vector<UIVertex> vertices;
        FontAtlas empty_atlas;
        UIGeometryBuilder builder(vertices, empty_atlas);
        element.draw(builder);
        return vertices.empty() ? -1.0f : vertices[0].r;
    }

    static constexpr Vector2F ON_CHECKBOX = {10.0f, 10.0f};
    static constexpr Vector2F FAR_AWAY = {900.0f, 900.0f};

    UIRoot root;
    UICheckbox* checkbox = nullptr;
    std::vector<bool> reported;        // every state the toggled action was called with
    std::vector<UISound> sounds_heard; // every sound that reached the root
};

// --------- TESTS --------- //

TEST_F(UICheckboxTest, StartsAsConstructed)
{
    EXPECT_FALSE(checkbox->isChecked());

    UICheckbox starts_checked("Collidable", true);
    EXPECT_TRUE(starts_checked.isChecked());
}

TEST_F(UICheckboxTest, ReleaseInsideFlipsReportsAndPlaysItsSound)
{
    clickInside();

    EXPECT_TRUE(checkbox->isChecked());
    EXPECT_EQ(reported, (std::vector<bool>{true})); // the action gets the NEW state
    EXPECT_EQ(sounds_heard, (std::vector<UISound>{UISound::Checkbox}));

    clickInside();

    EXPECT_FALSE(checkbox->isChecked());
    EXPECT_EQ(reported, (std::vector<bool>{true, false}));
}

TEST_F(UICheckboxTest, ReleaseOutsideCancels)
{
    checkbox->onMousePressed(ON_CHECKBOX);
    checkbox->onMouseReleased(FAR_AWAY, /*released_inside=*/false);

    EXPECT_FALSE(checkbox->isChecked());
    EXPECT_TRUE(reported.empty());
    EXPECT_TRUE(sounds_heard.empty());
}

TEST_F(UICheckboxTest, SetCheckedDoesNotReport)
{
    checkbox->setChecked(true);

    EXPECT_TRUE(checkbox->isChecked());
    EXPECT_TRUE(reported.empty());
    EXPECT_TRUE(sounds_heard.empty());
}

TEST_F(UICheckboxTest, SourceDecidesWhatIsShown)
{
    bool setting = true;
    checkbox->setCheckedSource([&setting]() { return setting; });

    root.updateElements(0.016f);
    EXPECT_TRUE(checkbox->isChecked());

    setting = false;
    root.updateElements(0.016f);
    EXPECT_FALSE(checkbox->isChecked());
}

// The source is the one source of truth: if the toggled action doesn't change what the source reads,
// the click only lasts until the next update pulls the old value back
TEST_F(UICheckboxTest, ActionThatIgnoresTheSourceIsUndoneNextFrame)
{
    const bool setting = false; // never changed by anything
    checkbox->setCheckedSource([&setting]() { return setting; });

    clickInside();
    EXPECT_TRUE(checkbox->isChecked()); // flipped by the click...

    root.updateElements(0.016f);
    EXPECT_FALSE(checkbox->isChecked()); // ...and pulled back by the source
}

TEST_F(UICheckboxTest, CheckedDrawsOneMoreRect)
{
    const size_t unchecked_vertices = vertexCount(*checkbox);

    checkbox->setChecked(true);

    EXPECT_EQ(vertexCount(*checkbox), unchecked_vertices + 6); // the white fill: one rectangle = 2 triangles
}

TEST_F(UICheckboxTest, DarkensWhenHoveredAndMoreWhenHeld)
{
    const float normal = firstVertexRed(*checkbox);
    checkbox->onMouseEntered();
    const float hovered = firstVertexRed(*checkbox);
    checkbox->onMousePressed(ON_CHECKBOX);
    const float held = firstVertexRed(*checkbox);

    EXPECT_LT(hovered, normal);
    EXPECT_LT(held, hovered);
}

TEST_F(UICheckboxTest, TakesPressesAnywhereOnTheRow)
{
    // The far right of the 200-wide row: on the label's side, nowhere near the box
    EXPECT_EQ(root.topmostElementAt({190.0f, 20.0f}), checkbox);
}
