// File: Tests/DynamoEngine/Scenes/Test_SceneManager.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/Scenes/SceneManager.hpp"

// Standard Library Imports
#include <memory>
#include <string>
#include <vector>

using namespace DynamoEngine;

// --------- TEST HELPERS --------- //

// Scene names for these tests, written the same way the game will write its own
namespace TestScene
{
enum : SceneId
{
    Menu,
    Game,
    Pause,
};
} // namespace TestScene

// A scene that writes "<name> entered" / "<name> exited" into a shared log
class RecordingScene : public Scene
{
  public:
    RecordingScene(std::string name, std::vector<std::string>& log)
        : Scene(SceneSettings::menu()), m_name(std::move(name)), m_log(log)
    {
    }
    void onEnter() override { m_log.push_back(m_name + " entered"); }
    void onExit() override { m_log.push_back(m_name + " exited"); }

  private:
    std::string m_name;
    std::vector<std::string>& m_log;
};

class SceneManagerTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        scenes.addScene(TestScene::Menu, std::make_unique<RecordingScene>("menu", log));
        scenes.addScene(TestScene::Game, std::make_unique<RecordingScene>("game", log));
        scenes.addScene(TestScene::Pause, std::make_unique<RecordingScene>("pause", log));
    }

    SceneManager scenes;
    std::vector<std::string> log;
};

// --------- SETTINGS --------- //

TEST(SceneSettings, PresetsDescribeMenuAndSimulationScenes)
{
    EXPECT_FALSE(SceneSettings::menu().runs_simulation);
    EXPECT_FALSE(SceneSettings::menu().draws_world);
    EXPECT_TRUE(SceneSettings::simulation().runs_simulation);
    EXPECT_TRUE(SceneSettings::simulation().draws_world);
}

TEST(Scene, KeepsItsSettings)
{
    Scene scene(SceneSettings{.runs_simulation = false, .draws_world = true}); // e.g. a level editor
    EXPECT_FALSE(scene.settings().runs_simulation);
    EXPECT_TRUE(scene.settings().draws_world);
}

// --------- REGISTRATION --------- //

TEST_F(SceneManagerTest, RegisteredScenesAreKnown)
{
    EXPECT_TRUE(scenes.hasScene(TestScene::Game));
    EXPECT_FALSE(scenes.hasScene(42));
}

TEST_F(SceneManagerTest, NoSceneIsActiveUntilTheFirstSwitch) { EXPECT_FALSE(scenes.hasCurrentScene()); }

// --------- SWITCHING --------- //

TEST_F(SceneManagerTest, SwitchWaitsForApplyPendingSwitch)
{
    scenes.requestSwitch(TestScene::Menu);
    EXPECT_FALSE(scenes.hasCurrentScene()); // requested, not yet applied
    EXPECT_TRUE(log.empty());

    scenes.applyPendingSwitch();
    EXPECT_EQ(scenes.currentSceneId(), TestScene::Menu);
    EXPECT_EQ(log, (std::vector<std::string>{"menu entered"}));
}

TEST_F(SceneManagerTest, OldSceneExitsBeforeNewSceneEnters)
{
    scenes.requestSwitch(TestScene::Game);
    scenes.applyPendingSwitch();
    log.clear();

    scenes.requestSwitch(TestScene::Pause);
    scenes.applyPendingSwitch();

    EXPECT_EQ(log, (std::vector<std::string>{"game exited", "pause entered"}));
    EXPECT_EQ(scenes.currentSceneId(), TestScene::Pause);
}

TEST_F(SceneManagerTest, LatestRequestInAFrameWins)
{
    scenes.requestSwitch(TestScene::Game);
    scenes.requestSwitch(TestScene::Pause);
    scenes.applyPendingSwitch();

    EXPECT_EQ(scenes.currentSceneId(), TestScene::Pause);
    EXPECT_EQ(log, (std::vector<std::string>{"pause entered"})); // Game was never entered
}

TEST_F(SceneManagerTest, SwitchingToTheActiveSceneDoesNothing)
{
    scenes.requestSwitch(TestScene::Game);
    scenes.applyPendingSwitch();
    log.clear();

    scenes.requestSwitch(TestScene::Game);
    scenes.applyPendingSwitch();

    EXPECT_TRUE(log.empty()); // no exit/enter churn
}

TEST_F(SceneManagerTest, ApplyWithNoRequestDoesNothing)
{
    scenes.requestSwitch(TestScene::Menu);
    scenes.applyPendingSwitch();
    log.clear();

    scenes.applyPendingSwitch(); // a frame where nothing asked to switch

    EXPECT_TRUE(log.empty());
    EXPECT_EQ(scenes.currentSceneId(), TestScene::Menu);
}

// --------- MISTAKES --------- //
// EXPECT_DEBUG_DEATH: in Debug builds the call must stop the program (the assert fires);
// in Release builds it must simply run without crashing.

TEST_F(SceneManagerTest, SwitchingToAnUnknownSceneIsCaught)
{
    EXPECT_DEBUG_DEATH(scenes.requestSwitch(42), "unknown scene id");
}

TEST_F(SceneManagerTest, RegisteringAnIdTwiceIsCaught)
{
    EXPECT_DEBUG_DEATH(scenes.addScene(TestScene::Game, std::make_unique<Scene>(SceneSettings::menu())),
                       "duplicate scene id");
}
