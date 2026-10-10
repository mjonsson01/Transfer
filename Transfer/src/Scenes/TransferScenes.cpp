// File: Transfer/src/Scenes/TransferScenes.cpp

#include "Scenes/TransferScenes.hpp"

// Custom Imports
#include "Core/DEPRECATED_InputState.hpp"
#include "DynamoEngine/UI/Widgets/UIButton.hpp"
#include "DynamoEngine/UI/Widgets/UICheckbox.hpp"
#include "DynamoEngine/UI/Widgets/UIDropdown.hpp"
#include "DynamoEngine/UI/Widgets/UILabel.hpp"
#include "DynamoEngine/UI/Widgets/UIRow.hpp"
#include "DynamoEngine/UI/Widgets/UISlider.hpp"
#include "Utilities/Constants/EngineConstants.hpp"

// Standard Library Imports
#include <cmath>
#include <memory>
#include <string>
#include <utility>

namespace
{
using DynamoEngine::Scene;
using DynamoEngine::SceneManager;
using DynamoEngine::SceneSettings;
using DynamoEngine::SliderMapping;
using DynamoEngine::UIAlign;
using DynamoEngine::UIButton;
using DynamoEngine::UICheckbox;
using DynamoEngine::UIDropdown;
using DynamoEngine::UILabel;
using DynamoEngine::UIRow;
using DynamoEngine::UISlider;
using DynamoEngine::UISound;
using DynamoEngine::Vector2F;

// The slider row sits where the old sliders were: tracks 4/5 of the way down a 720p window, spread across it
constexpr Vector2F SLIDER_SIZE = {300.0f, 60.0f}; // knob (30) + the "Label: value" text underneath
constexpr float SLIDER_SPACING = 125.0f;
constexpr float SLIDER_ROW_BOTTOM_MARGIN = 93.0f;

// --- Shared pieces --- //

// Transfer's sound effect for each UI sound. (The AudioSystem already skips them when sound effects are off.)
void queueUISound(UIState& ui_state, UISound sound)
{
    switch (sound)
    {
    case UISound::Click:
        ui_state.QueueSoundEffect("ButtonClick");
        break;
    case UISound::Tick:
        ui_state.QueueSoundEffect("SliderTick");
        break;
    case UISound::Checkbox:
        ui_state.QueueSoundEffect("CheckboxToggle");
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

// --- The scenes --- //

std::unique_ptr<Scene> buildStartMenuScene(SceneManager& scenes, UIState& ui_state)
{
    std::unique_ptr<Scene> scene = makeScene(SceneSettings::menu(), ui_state);

    std::unique_ptr<UIButton> play_button = std::make_unique<UIButton>("Play Game");
    play_button->setPlacement({.align = UIAlign::Center, .size = {300.0f, 200.0f}});
    play_button->setOnClick([&scenes]() { scenes.requestSwitch(TransferScene::Game); });
    scene->ui().addChild(std::move(play_button));

    return scene;
}

std::unique_ptr<Scene> buildGameScene(UIState& ui_state, GameState& game_state)
{
    std::unique_ptr<Scene> scene = makeScene(SceneSettings::simulation(), ui_state);
    DEPRECATED_InputState& choices = ui_state.getMutableDEPRECATED_InputState(); // where the slider values go

    // FPS counter, top-left
    std::unique_ptr<UILabel> fps_label = std::make_unique<UILabel>("FPS: ");
    fps_label->setPlacement({.align = UIAlign::TopLeft, .margin = 10.0f});
    fps_label->setTextSource([&ui_state]() { return "FPS: " + std::to_string(static_cast<int>(ui_state.getFPS())); });
    scene->ui().addChild(std::move(fps_label));

    // Visor Dropdown, shows what bodies are colored by
    // camera_state.visor_view is the source of truth, the dropdown reads it every frame
    // and makes requests to change it
    CameraState& camera_state = game_state.getCameraStateMutable();
    std::unique_ptr<UIDropdown> visor_dropdown = std::make_unique<UIDropdown>(
        "Visor", std::vector<std::string>{"Realistic", "Mass", "Charge", "Temperature"}); // Same order as VisorView
    visor_dropdown->setPlacement({.align = UIAlign::TopRight, .size = {260.0f, 36.0f}, .margin = 10.0f});
    visor_dropdown->setOnOptionChosen([&camera_state](int option_index)
                                      { camera_state.visor_view = static_cast<VisorView>(option_index); });
    visor_dropdown->setSelectionSource([&camera_state]() { return static_cast<int>(camera_state.visor_view); });
    scene->ui().addChild(std::move(visor_dropdown));

    // The three sliders, left to right: simulation speed, radius, mass
    std::unique_ptr<UISlider> speed_slider = std::make_unique<UISlider>(
        "Simulation Speed", SliderMapping::linear(MIN_TIME_SCALE_FACTOR, MAX_TIME_SCALE_FACTOR),
        REGULAR_TIME_SCALE_FACTOR);
    speed_slider->setPlacement({.size = SLIDER_SIZE});
    speed_slider->setOnValueChanged([&choices](double value) { choices.selectedSimSpeedScale = value; });

    std::unique_ptr<UISlider> radius_slider =
        std::make_unique<UISlider>("Radius", SliderMapping::linear(0.0, MAX_RADIUS - 100.0), 0.0);
    radius_slider->setPlacement({.size = SLIDER_SIZE});
    radius_slider->setOnValueChanged([&choices](double value) { choices.selectedRadius = value; });

    std::unique_ptr<UISlider> mass_slider = std::make_unique<UISlider>("Mass", massMapping(), 0.0);
    mass_slider->setPlacement({.size = SLIDER_SIZE});
    mass_slider->setOnValueChanged([&choices](double value) { choices.selectedMass = value; });

    std::unique_ptr<UIRow> slider_row = std::make_unique<UIRow>(SLIDER_SPACING);
    slider_row->setPlacement({.align = UIAlign::BottomCenter, .margin = SLIDER_ROW_BOTTOM_MARGIN});
    slider_row->addChild(std::move(speed_slider));
    slider_row->addChild(std::move(radius_slider));
    slider_row->addChild(std::move(mass_slider));
    scene->ui().addChild(std::move(slider_row));

    return scene;
}

std::unique_ptr<Scene> buildPauseScene(SceneManager& scenes, UIState& ui_state)
{
    std::unique_ptr<Scene> scene = makeScene(SceneSettings::menu(), ui_state);

    std::unique_ptr<UIButton> resume_button = std::make_unique<UIButton>("Resume");
    resume_button->setPlacement({.align = UIAlign::Center, .size = {200.0f, 100.0f}});
    resume_button->setOnClick([&scenes]() { scenes.requestSwitch(TransferScene::Game); });
    scene->ui().addChild(std::move(resume_button));

    return scene;
}

// A test bed: the simulation runs and takes game input, but nothing is drawn yet
std::unique_ptr<Scene> buildTestVisualScene(UIState& ui_state)
{
    std::unique_ptr<Scene> scene = makeScene(SceneSettings{.runs_simulation = true, .draws_world = false}, ui_state);

    // A lone checkbox to look at and click (its sound, hover and held colors). It isn't wired to anything yet:
    // it just logs its new state.
    std::unique_ptr<UICheckbox> test_checkbox = std::make_unique<UICheckbox>("Test Checkbox", false);
    test_checkbox->setPlacement({.align = UIAlign::Center, .size = {260.0f, 40.0f}});
    test_checkbox->setOnToggled([](bool is_checked)
                                { SDL_Log("Test checkbox is now %s", is_checked ? "checked" : "unchecked"); });
    scene->ui().addChild(std::move(test_checkbox));

    return scene;
}
} // namespace

void addTransferScenes(SceneManager& scenes, UIState& ui_state, GameState& game_state)
{
    scenes.addScene(TransferScene::StartMenu, buildStartMenuScene(scenes, ui_state));
    scenes.addScene(TransferScene::Game, buildGameScene(ui_state, game_state));
    scenes.addScene(TransferScene::Pause, buildPauseScene(scenes, ui_state));
    scenes.addScene(TransferScene::TestVisual, buildTestVisualScene(ui_state));
}
