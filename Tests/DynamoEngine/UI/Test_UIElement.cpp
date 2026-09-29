// File: Tests/DynamoEngine/UI/Test_UIElement.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/UI/UIElement.hpp"

// Standard Library Imports
#include <memory>
#include <utility>

using namespace DynamoEngine;

// --------- CHILDREN --------- //

TEST(UIElement, AddChildTakesOwnershipAndSetsParent)
{
    UIElement parent;
    auto child = std::make_unique<UIElement>();
    UIElement* child_address = child.get();

    UIElement& added = parent.addChild(std::move(child));

    EXPECT_EQ(child, nullptr);        // we handed it over
    EXPECT_EQ(&added, child_address); // the returned reference is the child we added
    EXPECT_EQ(added.parent(), &parent);
    ASSERT_EQ(parent.children().size(), 1u);
    EXPECT_EQ(parent.children()[0].get(), child_address);
    EXPECT_EQ(parent.parent(), nullptr); // the top element has no parent
}

// --------- LAYERS --------- //

TEST(UIElement, RootDefaultsToHUD)
{
    UIElement root;
    EXPECT_EQ(root.layer(), UILayer::HUD);
}

TEST(UIElement, ChildInheritsParentLayerUnlessSet)
{
    UIElement panel;
    panel.setLayer(UILayer::Menu);
    UIElement& button = panel.addChild(std::make_unique<UIElement>());
    UIElement& dropdown_list = panel.addChild(std::make_unique<UIElement>());
    dropdown_list.setLayer(UILayer::Overlay);

    EXPECT_EQ(button.layer(), UILayer::Menu);           // inherited
    EXPECT_EQ(dropdown_list.layer(), UILayer::Overlay); // its own
}

TEST(UIElement, LayerInheritsThroughSeveralLevels)
{
    UIElement root;
    root.setLayer(UILayer::Menu);
    UIElement& middle = root.addChild(std::make_unique<UIElement>());
    UIElement& leaf = middle.addChild(std::make_unique<UIElement>());
    EXPECT_EQ(leaf.layer(), UILayer::Menu);
}

// --------- PLACEMENT --------- //
// Parent: x 100..500, y 50..250 (400 x 200). Element: 40 x 20. Margin: 10.

class UIElementPlacementTest : public ::testing::Test
{
  protected:
    SDL_FRect placed(UIAlign align)
    {
        UIElement element;
        element.setPlacement({.align = align, .size = {40.0f, 20.0f}, .margin = 10.0f});
        element.updateLayout(PARENT_RECT);
        return element.rect();
    }

    static constexpr SDL_FRect PARENT_RECT = {100.0f, 50.0f, 400.0f, 200.0f};
};

TEST_F(UIElementPlacementTest, EveryAlignmentLandsInTheRightSpot)
{
    struct Case
    {
        UIAlign align;
        float expected_x;
        float expected_y;
    };
    // x: left = 100+10, center = 100+200-20 (margin ignored), right = 500-40-10
    // y: top  =  50+10, center =  50+100-10 (margin ignored), bottom = 250-20-10
    const Case cases[] = {
        {UIAlign::TopLeft, 110.0f, 60.0f},      {UIAlign::TopCenter, 280.0f, 60.0f},
        {UIAlign::TopRight, 450.0f, 60.0f},     {UIAlign::CenterLeft, 110.0f, 140.0f},
        {UIAlign::Center, 280.0f, 140.0f},      {UIAlign::CenterRight, 450.0f, 140.0f},
        {UIAlign::BottomLeft, 110.0f, 220.0f},  {UIAlign::BottomCenter, 280.0f, 220.0f},
        {UIAlign::BottomRight, 450.0f, 220.0f},
    };

    for (const Case& test_case : cases)
    {
        SCOPED_TRACE(static_cast<int>(test_case.align));
        SDL_FRect rect = placed(test_case.align);
        EXPECT_FLOAT_EQ(rect.x, test_case.expected_x);
        EXPECT_FLOAT_EQ(rect.y, test_case.expected_y);
        EXPECT_FLOAT_EQ(rect.w, 40.0f); // size is always kept
        EXPECT_FLOAT_EQ(rect.h, 20.0f);
    }
}

TEST_F(UIElementPlacementTest, FillCoversTheParentInsetByTheMargin)
{
    SDL_FRect rect = placed(UIAlign::Fill);
    EXPECT_FLOAT_EQ(rect.x, 110.0f);
    EXPECT_FLOAT_EQ(rect.y, 60.0f);
    EXPECT_FLOAT_EQ(rect.w, 380.0f); // 400 - 2 * 10
    EXPECT_FLOAT_EQ(rect.h, 180.0f); // 200 - 2 * 10
}

TEST(UIElement, ChildrenArePlacedInsideTheirParent)
{
    UIElement panel;
    panel.setPlacement({.align = UIAlign::TopLeft, .size = {200.0f, 100.0f}, .margin = 0.0f});
    UIElement& button = panel.addChild(std::make_unique<UIElement>());
    button.setPlacement({.align = UIAlign::BottomRight, .size = {50.0f, 20.0f}, .margin = 5.0f});

    panel.updateLayout(SDL_FRect{0.0f, 0.0f, 1280.0f, 720.0f}); // "the window"

    EXPECT_FLOAT_EQ(button.rect().x, 145.0f); // 200 - 50 - 5, relative to the PANEL, not the window
    EXPECT_FLOAT_EQ(button.rect().y, 75.0f);  // 100 - 20 - 5
}

// --------- HIT TESTING --------- //

TEST(UIElement, ContainsPointIncludesTopLeftEdgeButNotBottomRightEdge)
{
    UIElement element;
    element.setPlacement({.align = UIAlign::TopLeft, .size = {100.0f, 50.0f}, .margin = 0.0f});
    element.updateLayout(SDL_FRect{0.0f, 0.0f, 1280.0f, 720.0f});

    EXPECT_TRUE(element.containsPoint({50.0f, 25.0f}));   // middle
    EXPECT_TRUE(element.containsPoint({0.0f, 0.0f}));     // top-left corner: inside
    EXPECT_FALSE(element.containsPoint({100.0f, 25.0f})); // right edge: outside
    EXPECT_FALSE(element.containsPoint({50.0f, 50.0f}));  // bottom edge: outside
    EXPECT_FALSE(element.containsPoint({-1.0f, 25.0f}));  // left of it
}

// --------- DEFAULTS --------- //

TEST(UIElement, DefaultsAreVisibleZeroZIndexAndIgnoreClicks)
{
    UIElement element;
    EXPECT_TRUE(element.isVisible());
    EXPECT_EQ(element.zIndex(), 0);
    EXPECT_FALSE(element.onMousePressed({0.0f, 0.0f})); // a plain element never takes a click
}
