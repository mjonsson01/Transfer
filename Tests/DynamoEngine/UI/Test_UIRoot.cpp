// File: Tests/DynamoEngine/UI/Test_UIRoot.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/Rendering/FontAtlas.hpp"
#include "DynamoEngine/UI/UIGeometryBuilder.hpp"
#include "DynamoEngine/UI/UIRoot.hpp"
#include "TestingUtilities/VectorsNear.hpp"

// Standard Library Imports
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace DynamoEngine;

// --------- TEST HELPERS --------- //

// An element with a name that writes that name into a shared log whenever it is drawn or updated,
// so tests can check the ORDER things happened in.
class LoggingElement : public UIElement
{
  public:
    LoggingElement(std::string name, std::vector<std::string>& log) : m_name(std::move(name)), m_log(log) {}

    void draw(UIGeometryBuilder& builder) const override { m_log.push_back("draw " + m_name); }
    void update(float delta_seconds) override { m_log.push_back("update " + m_name); }

  private:
    std::string m_name;
    std::vector<std::string>& m_log;
};

// Adds a full-window LoggingElement to `parent` and returns it, so tests can build trees in one line each
LoggingElement& addLogger(UIElement& parent, const std::string& name, std::vector<std::string>& log)
{
    auto element = std::make_unique<LoggingElement>(name, log);
    element->setPlacement({.align = UIAlign::Fill});
    LoggingElement& added = *element;
    parent.addChild(std::move(element));
    return added;
}

// Draws the whole UI and returns just the order of names drawn
std::vector<std::string> drawOrder(UIRoot& root, std::vector<std::string>& log)
{
    std::vector<UIVertex> vertices;
    FontAtlas empty_atlas;
    UIGeometryBuilder builder(vertices, empty_atlas);

    log.clear();
    root.drawElements(builder);
    return log;
}

// --------- SCALE --------- //

TEST(UIRoot, ScaleIsOneAtTheReferenceHeight)
{
    UIRoot root;
    root.setWindowSize({1280.0f, 720.0f});
    EXPECT_FLOAT_EQ(root.uiScale(), 1.0f);
    EXPECT_TRUE(VectorsNear(root.uiSpaceSize(), {1280.0f, 720.0f}, 1e-4f));
}

TEST(UIRoot, ScaleGrowsWithThewindow_height)
{
    UIRoot root;
    root.setWindowSize({2560.0f, 1440.0f});
    EXPECT_FLOAT_EQ(root.uiScale(), 2.0f);
    EXPECT_TRUE(VectorsNear(root.uiSpaceSize(), {1280.0f, 720.0f}, 1e-4f)); // same layout space as 720p
}

TEST(UIRoot, TallerWindowFitsTheWidth)
{
    // A 16:10-ish window (a MacBook in fullscreen): width is the tight direction, so the scale follows it and the
    // extra height becomes extra UI space at the top and bottom
    UIRoot root;
    root.setWindowSize({1512.0f, 982.0f});
    EXPECT_FLOAT_EQ(root.uiScale(), 1512.0f / 1280.0f);
    EXPECT_FLOAT_EQ(root.uiSpaceSize().x_val, 1280.0f);
    EXPECT_GT(root.uiSpaceSize().y_val, 720.0f);
}

TEST(UIRoot, WiderWindowFitsTheHeight)
{
    // An ultrawide window: height is the tight direction, the extra width becomes extra UI space at the sides
    UIRoot root;
    root.setWindowSize({3440.0f, 1440.0f});
    EXPECT_FLOAT_EQ(root.uiScale(), 2.0f);
    EXPECT_TRUE(VectorsNear(root.uiSpaceSize(), {1720.0f, 720.0f}, 1e-4f));
}

TEST(UIRoot, DesignSizedElementFitsInAnyWindowShape)
{
    // Anything that fits in 1280 x 720 must fit in every window, whatever its shape
    const Vector2F window_sizes[] = {{1280.0f, 720.0f}, {1512.0f, 982.0f}, {800.0f, 800.0f}, {3440.0f, 1440.0f}};
    for (const Vector2F& window_size : window_sizes)
    {
        UIRoot root;
        root.setWindowSize(window_size);
        auto panel = std::make_unique<UIElement>();
        panel->setPlacement({.align = UIAlign::BottomCenter, .size = {1280.0f, 720.0f}, .margin = 0.0f});
        UIElement& added = root.addChild(std::move(panel));

        root.updateElements(0.016f);

        const Vector2F space = root.uiSpaceSize();
        EXPECT_GE(added.rect().x, -1e-3f);
        EXPECT_GE(added.rect().y, -1e-3f);
        EXPECT_LE(added.rect().x + added.rect().w, space.x_val + 1e-3f);
        EXPECT_LE(added.rect().y + added.rect().h, space.y_val + 1e-3f);
    }
}

TEST(UIRoot, PlayerScaleMultipliesTheAutoFit)
{
    UIRoot root;
    root.setWindowSize({2560.0f, 1440.0f});
    root.setPlayerScale(1.5f);
    EXPECT_FLOAT_EQ(root.uiScale(), 3.0f); // 2.0 auto-fit * 1.5 player setting
}

TEST(UIRoot, InvalidPlayerScaleIsIgnored)
{
    UIRoot root;
    root.setPlayerScale(0.0f);
    EXPECT_FLOAT_EQ(root.uiScale(), 1.0f);
}

TEST(UIRoot, MinimizedWindowDoesNotDivideByZero)
{
    UIRoot root;
    root.setWindowSize({0.0f, 0.0f});
    EXPECT_GT(root.uiScale(), 0.0f);
    root.setWindowSize({1280.0f, 0.0f}); // only one side collapsed
    EXPECT_GT(root.uiScale(), 0.0f);
    root.setWindowSize({0.0f, 720.0f});
    EXPECT_GT(root.uiScale(), 0.0f);
}

TEST(UIRoot, ScreenToUISpaceDividesByTheScale)
{
    UIRoot root;
    root.setWindowSize({2560.0f, 1440.0f});
    EXPECT_TRUE(VectorsNear(root.screenToUISpace({1000.0f, 500.0f}), {500.0f, 250.0f}, 1e-4f));
}

TEST(UIRoot, ElementsAreLaidOutInUISpace)
{
    UIRoot root;
    root.setWindowSize({2560.0f, 1440.0f}); // scale 2 -> UI space is 1280 x 720
    auto button = std::make_unique<UIElement>();
    button->setPlacement({.align = UIAlign::BottomRight, .size = {100.0f, 50.0f}, .margin = 0.0f});
    UIElement& added = root.addChild(std::move(button));

    root.updateElements(0.016f);

    EXPECT_FLOAT_EQ(added.rect().x, 1180.0f); // 1280 - 100, in UI space (not 2560 - 100)
    EXPECT_FLOAT_EQ(added.rect().y, 670.0f);  // 720 - 50
}

// --------- DRAW ORDER --------- //

TEST(UIRoot, ParentsDrawBeforeTheirChildren)
{
    std::vector<std::string> log;
    UIRoot root;
    LoggingElement& panel = addLogger(root, "panel", log);
    addLogger(panel, "button", log);

    EXPECT_EQ(drawOrder(root, log), (std::vector<std::string>{"draw panel", "draw button"}));
}

TEST(UIRoot, SiblingsDrawByZIndexAndTiesKeepTheOrderAdded)
{
    std::vector<std::string> log;
    UIRoot root;
    addLogger(root, "a", log).setZIndex(1);
    addLogger(root, "b", log);
    addLogger(root, "c", log);

    // b and c (z 0) keep their added order; a (z 1) goes on top
    EXPECT_EQ(drawOrder(root, log), (std::vector<std::string>{"draw b", "draw c", "draw a"}));
}

TEST(UIRoot, HigherLayerDrawsAboveEverythingInLowerLayers)
{
    std::vector<std::string> log;
    UIRoot root;
    LoggingElement& visor_menu = addLogger(root, "visor_menu", log);
    addLogger(visor_menu, "visor_list", log).setLayer(UILayer::Overlay);
    addLogger(root, "slider", log).setZIndex(99); // even a huge z-index can't beat a higher layer

    EXPECT_EQ(drawOrder(root, log), (std::vector<std::string>{"draw visor_menu", "draw slider", "draw visor_list"}));
}

TEST(UIRoot, HiddenElementHidesItsWholeSubtree)
{
    std::vector<std::string> log;
    UIRoot root;
    LoggingElement& panel = addLogger(root, "panel", log);
    addLogger(panel, "button", log);
    addLogger(root, "fps", log);
    panel.setVisible(false);

    EXPECT_EQ(drawOrder(root, log), (std::vector<std::string>{"draw fps"}));
}

TEST(UIRoot, UpdateReachesEveryVisibleElement)
{
    std::vector<std::string> log;
    UIRoot root;
    LoggingElement& panel = addLogger(root, "panel", log);
    addLogger(panel, "button", log);
    addLogger(root, "hidden", log).setVisible(false);

    root.updateElements(0.016f);

    EXPECT_EQ(log, (std::vector<std::string>{"update panel", "update button"}));
}

// --------- HIT TESTING --------- //

TEST(UIRoot, TopmostElementIsTheOneDrawnLast)
{
    std::vector<std::string> log;
    UIRoot root;
    LoggingElement& bottom = addLogger(root, "bottom", log);
    LoggingElement& top = addLogger(root, "top", log); // both cover the whole window
    root.updateElements(0.016f);

    EXPECT_EQ(root.topmostElementAt({10.0f, 10.0f}), &top);
    (void)bottom;
}

TEST(UIRoot, OverlayWinsTheClickOverLaterHUDElements)
{
    std::vector<std::string> log;
    UIRoot root;
    LoggingElement& visor_menu = addLogger(root, "visor_menu", log);
    LoggingElement& visor_list = addLogger(visor_menu, "visor_list", log);
    visor_list.setLayer(UILayer::Overlay);
    addLogger(root, "slider", log); // added later, same spot, but on HUD
    root.updateElements(0.016f);

    EXPECT_EQ(root.topmostElementAt({10.0f, 10.0f}), &visor_list);
}

TEST(UIRoot, HiddenElementsAreNeverHit)
{
    std::vector<std::string> log;
    UIRoot root;
    LoggingElement& visible = addLogger(root, "visible", log);
    addLogger(root, "hidden_on_top", log).setVisible(false);
    root.updateElements(0.016f);

    EXPECT_EQ(root.topmostElementAt({10.0f, 10.0f}), &visible);
}

TEST(UIRoot, NothingUnderThePointGivesNull)
{
    UIRoot root;
    auto small = std::make_unique<UIElement>();
    small->setPlacement({.align = UIAlign::TopLeft, .size = {10.0f, 10.0f}, .margin = 0.0f});
    root.addChild(std::move(small));
    root.updateElements(0.016f);

    EXPECT_EQ(root.topmostElementAt({500.0f, 500.0f}), nullptr);
}
