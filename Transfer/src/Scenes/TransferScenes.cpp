// File: Transfer/src/Scenes/TransferScenes.cpp

#include "Scenes/TransferScenes.hpp"

// Custom Imports
#include "Constants/EngineConstants.hpp"
#include "DynamoEngine/Constants/GlobalConstants.hpp"
#include "DynamoEngine/UI/Widgets/UIButton.hpp"
#include "DynamoEngine/UI/Widgets/UICheckbox.hpp"
#include "DynamoEngine/UI/Widgets/UIDropdown.hpp"
#include "DynamoEngine/UI/Widgets/UILabel.hpp"
#include "DynamoEngine/UI/Widgets/UIRow.hpp"
#include "DynamoEngine/UI/Widgets/UISlider.hpp"
#include "Systems/Physics/PhysicsConstants.hpp"

// Standard Library Imports
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace
{
using DynamoEngine::Scene;
using DynamoEngine::SceneManager;
using DynamoEngine::SceneSettings;
using DynamoEngine::SliderMapping;
using DynamoEngine::UIAlign;
using DynamoEngine::UIButton;
using DynamoEngine::UICheckbox;
using DynamoEngine::UIColumn;
using DynamoEngine::UIDropdown;
using DynamoEngine::UIElement;
using DynamoEngine::UILabel;
using DynamoEngine::UIRow;
using DynamoEngine::UISlider;
using DynamoEngine::UISound;
using DynamoEngine::Vector2D;
using DynamoEngine::Vector2F;

// The bottom of the screen, from the bottom edge up: the spawn panel, then the slider row above it.
// The spawn panel is titled sections side by side: each one a title line over a row of checkboxes.
constexpr Vector2F CHECKBOX_SIZE = {200.0f, 36.0f}; // as tall as the visor dropdown; its label gets 164 of the 200
constexpr float CHECKBOX_SPACING = 12.0f;           // between neighboring checkboxes in a section
constexpr float SECTION_TITLE_HEIGHT = 27.0f;       // one line of the 18 pt UI font
constexpr float READOUT_LINE_HEIGHT = 27.0f;        // the top-left text lines (FPS, planning), same font
constexpr float SECTION_TITLE_GAP = 4.0f;           // between a section's title and its checkboxes
constexpr float SPAWN_SECTION_SPACING = 40.0f;      // between neighboring sections
constexpr float CHECKBOX_ROW_SPACING = 8.0f;        // between the rows of a section with more than one row
constexpr int SPAWN_PANEL_CHECKBOX_ROWS = 1;        // the tallest section's rows (all one row today)
constexpr float SPAWN_PANEL_HEIGHT = SECTION_TITLE_HEIGHT + SECTION_TITLE_GAP +
                                     SPAWN_PANEL_CHECKBOX_ROWS * CHECKBOX_SIZE.y_val +
                                     (SPAWN_PANEL_CHECKBOX_ROWS - 1) * CHECKBOX_ROW_SPACING;
constexpr float SPAWN_PANEL_BOTTOM_MARGIN = 20.0f;
constexpr float SLIDERS_TO_SPAWN_PANEL_GAP = 30.0f; // empty space between the sliders' text and the section titles

constexpr Vector2F SLIDER_SIZE = {300.0f, 60.0f}; // knob (30) + the "Label: value" text underneath
constexpr float SLIDER_SPACING = 125.0f;
// Derived from the spawn panel, so the sliders always sit above it and the two can never overlap
constexpr float SLIDER_ROW_BOTTOM_MARGIN = SPAWN_PANEL_BOTTOM_MARGIN + SPAWN_PANEL_HEIGHT + SLIDERS_TO_SPAWN_PANEL_GAP;

// The top corners: the readouts (top-left), the visor dropdown with the profiler overlay under it (top-right)
constexpr float TOP_CORNER_MARGIN = 10.0f; // from the window's edges
constexpr Vector2F VISOR_DROPDOWN_SIZE = {260.0f, 36.0f};
constexpr float PROFILER_GAP_BELOW_VISOR = 8.0f;
// The profiler's longest line is 46 characters ("%-19s %8.3f %8.3f %8.3f" in Game::profilerReportLines). The UI font,
// Space Mono, is monospaced: every character is 0.612 em wide, 11 points at 18 pt, so 46 of them take ~507 points.
constexpr float PROFILER_TEXT_WIDTH = 520.0f;
// The profiler column must end above the slider row. 1280 x 720 is the tightest the UI space ever gets (UIRoot::fitScale
// only ever adds space), so checking it there covers every window.
constexpr float PROFILER_COLUMN_BOTTOM = TOP_CORNER_MARGIN + VISOR_DROPDOWN_SIZE.y_val + PROFILER_GAP_BELOW_VISOR +
                                         static_cast<float>(MAX_PROFILER_LINES) * READOUT_LINE_HEIGHT;
constexpr float SLIDER_ROW_TOP = DynamoEngine::UI_REFERENCE_HEIGHT - SLIDER_ROW_BOTTOM_MARGIN - SLIDER_SIZE.y_val;
static_assert(PROFILER_COLUMN_BOTTOM < SLIDER_ROW_TOP,
              "The profiler overlay would overlap the sliders: lower MAX_PROFILER_LINES or move the sliders down");

// --- Shared pieces --- //

// Transfer's sound effect for each UI sound. (The AudioSystem already skips them when sound effects are off.)
void queueUISound(UIState& ui_state, UISound sound)
{
    switch (sound)
    {
    case UISound::SpecialClick:
        ui_state.QueueSoundEffect("ButtonClickSpecial");
        break;
    case UISound::Tick:
        ui_state.QueueSoundEffect("SliderTick");
        break;
    case UISound::StandardClick:
        ui_state.QueueSoundEffect("ButtonClickStandard");
        break;
    }
}

// A new scene whose UI sounds go to Transfer's sound effects
std::unique_ptr<Scene> makeScene(SceneSettings settings, UIState& ui_state)
{
    std::unique_ptr<Scene> scene = std::make_unique<Scene>(settings);
    scene->ui().setSoundHandler([&ui_state](UISound sound) { queueUISound(ui_state, sound); });
    return scene;
}

// The mass slider is centered on zero: left of center is negative mass, right of center is positive.
// Moving away from the center, the mass grows exponentially (10^...), so small masses get as much of the track
// as huge ones. Both directions of the curve live here, so they can never disagree.
SliderMapping massMapping()
{
    const double max_mass = MAX_MASS / 10.0;
    const double exponent = std::log10(max_mass); // at either end of the track: 10^exponent = max_mass

    SliderMapping mapping;
    mapping.to_value = [exponent](double position)
    {
        const double centered = (position * 2.0) - 1.0; // -1 (far left) .. 0 (center) .. 1 (far right)
        if (std::abs(centered) < 0.01)
        {
            return 0.0; // snap to zero around the center
        }
        const double sign = (centered < 0.0) ? -1.0 : 1.0;
        return sign * std::pow(10.0, std::abs(centered) * exponent);
    };
    mapping.to_position = [exponent](double value)
    {
        if (value == 0.0)
        {
            return 0.5; // zero sits in the center
        }
        const double sign = (value < 0.0) ? -1.0 : 1.0;
        const double centered = sign * std::log10(std::abs(value)) / exponent;
        return (centered + 1.0) / 2.0;
    };
    return mapping;
}

// One spawn-panel checkbox bound to one SpawnSettings flag: it shows the flag every frame and sets it when clicked.
// `settings` is ui_state's SpawnSettings, which outlives every scene, so the lambdas can keep a reference to it.
void addSpawnCheckbox(UIRow& checkbox_row, const std::string& label, SpawnSettings& settings, SpawnFlag flag)
{
    std::unique_ptr<UICheckbox> checkbox = std::make_unique<UICheckbox>(label, settings.flag(flag));
    checkbox->setPlacement({.size = CHECKBOX_SIZE});
    checkbox->setCheckedSource([&settings, flag]() { return settings.flag(flag); });
    // Through setFlag, so checking one box unchecks the boxes it makes meaningless (they show it next frame)
    checkbox->setOnToggled([&settings, flag](bool is_checked) { settings.setFlag(flag, is_checked); });
    checkbox_row.addChild(std::move(checkbox));
}

// One titled section of the spawn panel: the title on one line, its checkboxes underneath (one row, or a column of
// rows). The title must fit within the checkboxes' width (two checkboxes = 412): it's drawn from its top-left corner.
// A section with fewer rows than SPAWN_PANEL_CHECKBOX_ROWS leaves space underneath (sections line up at the top).
std::unique_ptr<UIColumn> makeSpawnSection(const std::string& title, std::unique_ptr<UIElement> checkboxes)
{
    std::unique_ptr<UILabel> title_label = std::make_unique<UILabel>(title);
    title_label->setPlacement({.size = {0.0f, SECTION_TITLE_HEIGHT}}); // 0 wide: the checkboxes set the width

    std::unique_ptr<UIColumn> section = std::make_unique<UIColumn>(SECTION_TITLE_GAP);
    section->addChild(std::move(title_label));
    section->addChild(std::move(checkboxes));
    return section;
}

// The planning mode line under the FPS counter (empty outside planning mode), with how far the trace has got
std::string planningReadout(const UIState& ui_state, const GameState& game_state)
{
    if (!ui_state.getIsPlanning())
    {
        return "";
    }
    const Trajectories& trajectories = game_state.getTrajectories();
    char text[96];
    switch (trajectories.status)
    {
    case TraceStatus::Idle:
        return "PLANNING: press I to integrate";
    case TraceStatus::Tracing:
        std::snprintf(text, sizeof(text), "PLANNING: Integrating... %.1f / %.0f s", trajectories.predicted_seconds,
                      TRAJECTORY_HORIZON_SECONDS);
        return text;
    case TraceStatus::Finished:
        std::snprintf(text, sizeof(text), "PLANNING: %.0f s predicted", trajectories.predicted_seconds);
        return text;
    }
    return "PLANNING";
}

// The line under that while a planet is being placed: how far its preview reaches (empty otherwise)
std::string previewReadout(const GameState& game_state)
{
    const Trajectories& trajectories = game_state.getTrajectories();
    if (trajectories.spawn_preview_path.empty())
    {
        return "";
    }
    char text[96];
    std::snprintf(text, sizeof(text), "PREVIEW: %.1f / %.0f s", trajectories.preview_predicted_seconds,
                  TRAJECTORY_HORIZON_SECONDS);
    return text;
}

// --- The scenes --- //

std::unique_ptr<Scene> buildStartMenuScene(SceneManager& scenes, UIState& ui_state)
{
    std::unique_ptr<Scene> scene = makeScene(SceneSettings::menu(), ui_state);

    std::unique_ptr<UIButton> play_button = std::make_unique<UIButton>("Play Game");
    play_button->setPlacement({.size = {300.0f, 200.0f}});
    play_button->setOnClick([&scenes]() { scenes.requestSwitch(TransferScene::Game); });

    // The testing sandbox (collision box, optional small reflective world)
    std::unique_ptr<UIButton> test_button = std::make_unique<UIButton>("Test Scene");
    test_button->setPlacement({.size = {300.0f, 60.0f}});
    test_button->setOnClick([&scenes]() { scenes.requestSwitch(TransferScene::TestVisual); });

    std::unique_ptr<UIColumn> buttons = std::make_unique<UIColumn>(20.0f);
    buttons->setPlacement({.align = UIAlign::Center});
    buttons->addChild(std::move(play_button));
    buttons->addChild(std::move(test_button));
    scene->ui().addChild(std::move(buttons));

    return scene;
}

// The sandbox's UI, shared by the game scene and the test scene: readouts (top-left), visor and profiler overlay
// (top-right), the sliders and the spawn panel (bottom)
void addSandboxUI(DynamoEngine::UIRoot& ui, UIState& ui_state, GameState& game_state)
{
    // Top-left, one line each: the FPS counter, the throttle, then the planning mode readouts
    std::unique_ptr<UILabel> fps_label = std::make_unique<UILabel>("FPS: ");
    fps_label->setPlacement({.size = {0.0f, READOUT_LINE_HEIGHT}});
    fps_label->setTextSource([&ui_state]() { return "FPS: " + std::to_string(static_cast<int>(ui_state.getFPS())); });

    // The ship's throttle: it holds when W/S are released, so it needs to be visible
    std::unique_ptr<UILabel> throttle_label = std::make_unique<UILabel>("");
    throttle_label->setPlacement({.size = {0.0f, READOUT_LINE_HEIGHT}});
    throttle_label->setTextSource(
        [&game_state]()
        {
            const UniverseState& universe_state = game_state.getUniverseState();
            if (!universe_state.has_ship)
            {
                return std::string(); // no ship, no throttle
            }
            int percent = static_cast<int>(std::lround(universe_state.ship.throttle() * 100.0));
            return "THROTTLE: " + std::to_string(percent) + "%";
        });

    std::unique_ptr<UILabel> planning_label = std::make_unique<UILabel>("");
    planning_label->setPlacement({.size = {0.0f, READOUT_LINE_HEIGHT}});
    planning_label->setTextSource([&ui_state, &game_state]() { return planningReadout(ui_state, game_state); });

    std::unique_ptr<UILabel> preview_label = std::make_unique<UILabel>("");
    preview_label->setPlacement({.size = {0.0f, READOUT_LINE_HEIGHT}});
    preview_label->setTextSource([&game_state]() { return previewReadout(game_state); });

    std::unique_ptr<UIColumn> readouts = std::make_unique<UIColumn>(0.0f);
    readouts->setPlacement({.align = UIAlign::TopLeft, .margin = TOP_CORNER_MARGIN});
    readouts->addChild(std::move(fps_label));
    readouts->addChild(std::move(throttle_label));
    readouts->addChild(std::move(planning_label));
    readouts->addChild(std::move(preview_label));

    // TEMPORARY (camera limit check, 2026-11-25), shown with F3: where the camera is, the world's edges, and how far
    // keepViewInsideWorld pushed the camera back this frame. Remove (with CameraState::debug_view_push_back) when done.
    std::unique_ptr<UILabel> camera_label = std::make_unique<UILabel>("");
    camera_label->setPlacement({.size = {0.0f, READOUT_LINE_HEIGHT}});
    camera_label->setTextSource(
        [&ui_state, &game_state]()
        {
            if (!ui_state.getRenderDebug())
            {
                return std::string();
            }
            const CameraState& camera = game_state.getCameraState();
            double view_width = camera.window_width / camera.zoom;
            double view_height = camera.window_height / camera.zoom;
            char text[160];
            std::snprintf(text, sizeof(text), "CAMERA: center (%.0f, %.0f)  zoom %.3f  view %.0f x %.0f",
                          view_width / 2.0 - camera.offset.x_val, view_height / 2.0 - camera.offset.y_val, camera.zoom,
                          view_width, view_height);
            return std::string(text);
        });
    readouts->addChild(std::move(camera_label));

    std::unique_ptr<UILabel> world_label = std::make_unique<UILabel>("");
    world_label->setPlacement({.size = {0.0f, READOUT_LINE_HEIGHT}});
    world_label->setTextSource(
        [&ui_state, &game_state]()
        {
            if (!ui_state.getRenderDebug())
            {
                return std::string();
            }
            const UniverseState& universe = game_state.getUniverseState();
            const Vector2D& moved = game_state.getCameraState().debug_view_push_back;
            char text[160];
            std::snprintf(text, sizeof(text), "WORLD: x [%.0f, %.0f] y [%.0f, %.0f]  pushed back by (%.1f, %.1f)",
                          universe.world_min_corner.x_val, universe.world_max_corner.x_val,
                          universe.world_min_corner.y_val, universe.world_max_corner.y_val, moved.x_val, moved.y_val);
            return std::string(text);
        });
    readouts->addChild(std::move(world_label));

    ui.addChild(std::move(readouts));

    // Visor Dropdown, shows what bodies are colored by
    // camera_state.visor_view is the source of truth, the dropdown reads it every frame
    // and makes requests to change it
    CameraState& camera_state = game_state.getCameraStateMutable();
    std::unique_ptr<UIDropdown> visor_dropdown = std::make_unique<UIDropdown>(
        "Visor", std::vector<std::string>{"Realistic", "Mass", "Charge", "Temperature"}); // Same order as VisorView
    visor_dropdown->setPlacement({.align = UIAlign::TopRight, .size = VISOR_DROPDOWN_SIZE});
    visor_dropdown->setOnOptionChosen([&camera_state](int option_index)
                                      { camera_state.visor_view = static_cast<VisorView>(option_index); });
    visor_dropdown->setSelectionSource([&camera_state]() { return static_cast<int>(camera_state.visor_view); });

    // Top-right column: the visor dropdown, then the profiler overlay (F4 frame / F5 physics) under it, clear of the
    // sliders (see the static_assert at the top). The dropdown's open list is drawn over the profiler text.
    std::unique_ptr<UIColumn> top_right = std::make_unique<UIColumn>(0.0f);
    top_right->setPlacement({.align = UIAlign::TopRight, .margin = TOP_CORNER_MARGIN});
    top_right->addChild(std::move(visor_dropdown));
    std::unique_ptr<UIElement> profiler_gap = std::make_unique<UIElement>(); // empty: just space
    profiler_gap->setPlacement({.size = {0.0f, PROFILER_GAP_BELOW_VISOR}});
    top_right->addChild(std::move(profiler_gap));
    // One label per line of the latest report. Wide enough for the longest line, so the column is that wide and the
    // text (drawn from each label's left edge) stays on screen.
    for (size_t line_index = 0; line_index < MAX_PROFILER_LINES; ++line_index)
    {
        std::unique_ptr<UILabel> profiler_label = std::make_unique<UILabel>("");
        profiler_label->setPlacement({.size = {PROFILER_TEXT_WIDTH, READOUT_LINE_HEIGHT}});
        profiler_label->setTextSource(
            [&ui_state, line_index]()
            {
                ProfilerOverlay overlay = ui_state.getProfilerOverlay();
                if (overlay == ProfilerOverlay::Off)
                {
                    return std::string();
                }
                const std::vector<std::string>& lines = (overlay == ProfilerOverlay::Frame)
                                                            ? ui_state.frameProfilerLines()
                                                            : ui_state.physicsProfilerLines();
                if (line_index >= lines.size())
                {
                    return std::string();
                }
                return lines[line_index];
            });
        top_right->addChild(std::move(profiler_label));
    }
    ui.addChild(std::move(top_right));

    // The three sliders, left to right: simulation speed, radius, mass
    std::unique_ptr<UISlider> speed_slider = std::make_unique<UISlider>(
        "Simulation Speed", SliderMapping::linear(MIN_TIME_SCALE_FACTOR, MAX_TIME_SCALE_FACTOR),
        REGULAR_TIME_SCALE_FACTOR);
    speed_slider->setPlacement({.size = SLIDER_SIZE});
    speed_slider->setOnValueChanged([&ui_state](double value) { ui_state.setSimulationSpeed(value); });

    std::unique_ptr<UISlider> radius_slider =
        std::make_unique<UISlider>("Radius", SliderMapping::linear(0.0, MAX_RADIUS - 100.0), 0.0);
    radius_slider->setPlacement({.size = SLIDER_SIZE});
    radius_slider->setOnValueChanged([&ui_state](double value) { ui_state.spawnSettings().radius = value; });

    std::unique_ptr<UISlider> mass_slider = std::make_unique<UISlider>("Mass", massMapping(), 0.0);
    mass_slider->setPlacement({.size = SLIDER_SIZE});
    mass_slider->setOnValueChanged([&ui_state](double value) { ui_state.spawnSettings().mass = value; });

    // Spawn panel, under the sliders: settings for the NEXT planet or cluster.
    // ui_state's SpawnSettings is the source of truth: each checkbox reads it every frame and writes it when clicked.
    SpawnSettings& spawn_settings = ui_state.spawnSettings();
    // Global: settings for planets AND clusters
    std::unique_ptr<UIRow> global_checkboxes = std::make_unique<UIRow>(CHECKBOX_SPACING);
    addSpawnCheckbox(*global_checkboxes, "Accretable", spawn_settings, SpawnFlag::Accretable);
    addSpawnCheckbox(*global_checkboxes, "Bounce", spawn_settings, SpawnFlag::Bounce);
    addSpawnCheckbox(*global_checkboxes, "Ghost", spawn_settings, SpawnFlag::Ghost);

    // Macro body: settings only planets have
    std::unique_ptr<UIRow> macro_checkboxes = std::make_unique<UIRow>(CHECKBOX_SPACING);
    addSpawnCheckbox(*macro_checkboxes, "Static", spawn_settings, SpawnFlag::Static);
    addSpawnCheckbox(*macro_checkboxes, "Shatterable", spawn_settings, SpawnFlag::Shatterable);

    std::unique_ptr<UIRow> spawn_panel = std::make_unique<UIRow>(SPAWN_SECTION_SPACING);
    spawn_panel->setPlacement({.align = UIAlign::BottomCenter, .margin = SPAWN_PANEL_BOTTOM_MARGIN});
    spawn_panel->addChild(makeSpawnSection("Global Instantiation Properties", std::move(global_checkboxes)));
    spawn_panel->addChild(makeSpawnSection("Macro Body Instantiation Properties", std::move(macro_checkboxes)));

    ui.addChild(std::move(spawn_panel));

    std::unique_ptr<UIRow> slider_row = std::make_unique<UIRow>(SLIDER_SPACING);
    slider_row->setPlacement({.align = UIAlign::BottomCenter, .margin = SLIDER_ROW_BOTTOM_MARGIN});
    slider_row->addChild(std::move(speed_slider));
    slider_row->addChild(std::move(radius_slider));
    slider_row->addChild(std::move(mass_slider));
    ui.addChild(std::move(slider_row));
}

std::unique_ptr<Scene> buildGameScene(UIState& ui_state, GameState& game_state)
{
    std::unique_ptr<Scene> scene = makeScene(SceneSettings::simulation(), ui_state);
    addSandboxUI(scene->ui(), ui_state, game_state);
    return scene;
}

std::unique_ptr<Scene> buildPauseScene(SceneManager& scenes, UIState& ui_state)
{
    std::unique_ptr<Scene> scene = makeScene(SceneSettings::menu(), ui_state);

    std::unique_ptr<UIButton> resume_button = std::make_unique<UIButton>("Resume");
    resume_button->setPlacement({.align = UIAlign::Center, .size = {200.0f, 100.0f}});
    // Back to the scene that was paused (the game or the test scene)
    resume_button->setOnClick([&scenes, &ui_state]() { scenes.requestSwitch(ui_state.getResumeScene()); });
    scene->ui().addChild(std::move(resume_button));

    return scene;
}

// A test bed: the simulation runs and takes game input, but nothing is drawn yet
// The test scene: the game's sandbox plus testing tools. Entering it puts a collision box around the view (and turns on
// REFLECTIVE_WORLD_WALLS if set); leaving it takes them away again, so the game scene never has them.
class TestVisualScene : public Scene
{
  public:
    TestVisualScene(UIState& ui_state, GameState& game_state)
        : Scene(SceneSettings::simulation()), m_game_state(game_state)
    {
        ui().setSoundHandler([&ui_state](UISound sound) { queueUISound(ui_state, sound); });
    }

    void onEnter() override
    {
        UniverseState& universe_state = m_game_state.getUniverseStateMutable();
        const CameraState& camera_state = m_game_state.getCameraState();
        // The middle of what's on screen right now, in world coordinates
        const Vector2D view_center = {camera_state.window_width / 2.0 / camera_state.zoom - camera_state.offset.x_val,
                                      camera_state.window_height / 2.0 / camera_state.zoom - camera_state.offset.y_val};

        const Vector2D box_half_size = {TEST_COLLISION_BOX_WIDTH / 2.0, TEST_COLLISION_BOX_HEIGHT / 2.0};
        universe_state.has_collision_box = true;
        universe_state.collision_box_min_corner = view_center - box_half_size;
        universe_state.collision_box_max_corner = view_center + box_half_size;

        // Remember whether the world's edges were walls, to put it back on exit
        m_saved_reflective_walls = universe_state.reflective_walls;
        universe_state.reflective_walls = REFLECTIVE_WORLD_WALLS;
    }

    void onExit() override
    {
        UniverseState& universe_state = m_game_state.getUniverseStateMutable();
        universe_state.has_collision_box = false;
        universe_state.reflective_walls = m_saved_reflective_walls;
    }

  private:
    GameState& m_game_state;
    bool m_saved_reflective_walls = false;
};

std::unique_ptr<Scene> buildTestVisualScene(UIState& ui_state, GameState& game_state)
{
    std::unique_ptr<TestVisualScene> scene = std::make_unique<TestVisualScene>(ui_state, game_state);
    addSandboxUI(scene->ui(), ui_state, game_state);

    // Say which scene this is, top-center
    std::unique_ptr<UILabel> test_label = std::make_unique<UILabel>("TEST SCENE: collision box");
    test_label->setPlacement({.align = UIAlign::TopCenter, .size = {0.0f, READOUT_LINE_HEIGHT}, .margin = 10.0f});
    scene->ui().addChild(std::move(test_label));

    return scene;
}
} // namespace

void addTransferScenes(SceneManager& scenes, UIState& ui_state, GameState& game_state)
{
    scenes.addScene(TransferScene::StartMenu, buildStartMenuScene(scenes, ui_state));
    scenes.addScene(TransferScene::Game, buildGameScene(ui_state, game_state));
    scenes.addScene(TransferScene::Pause, buildPauseScene(scenes, ui_state));
    scenes.addScene(TransferScene::TestVisual, buildTestVisualScene(ui_state, game_state));
}
