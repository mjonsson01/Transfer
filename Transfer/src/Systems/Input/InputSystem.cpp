// File: Transfer/src/Systems/Input/InputSystem.cpp

// Custom Includes
#include "Systems/Input/InputSystem.hpp"

// SDL Includes
#include <SDL3/SDL_init.h>
#include <SDL3/SDL_scancode.h>

// Blank namespace helper functions
namespace
{
void requestShutdown(GameState& game_state)
{
    game_state.SetPlaying(false);
    game_state.setIsShuttingDownAudioSystem(true);
}
void zoomAroundCursor(CameraState& camera_state, float scroll, DynamoEngine::Vector2D mouse_position)
{
    if (!firstWithinEpsilonOfSecond(scroll, 0.0))
    {
        DynamoEngine::Vector2D world_under_cursor = ScreenToWorldCoordinates(mouse_position, camera_state);
        DynamoEngine::Vector2D star_world_under_cursor =
            mouse_position / camera_state.zoom - camera_state.twinkling_star_offset;

        camera_state.zoom *= std::pow(1.1, scroll);
        camera_state.zoom = std::clamp(camera_state.zoom, MIN_ZOOM, MAX_ZOOM);

        camera_state.offset = mouse_position / camera_state.zoom - world_under_cursor;
        camera_state.twinkling_star_offset = mouse_position / camera_state.zoom - star_world_under_cursor;
    }
}
void panCamera(CameraState& camera_state, DynamoEngine::Vector2D screen_delta)
{
    camera_state.offset += screen_delta / camera_state.zoom;
    camera_state.twinkling_star_offset += (screen_delta / camera_state.zoom) * STAR_PARALLAX_FACTOR;
}
void clampCameraToStarField(CameraState& camera_state)
{
    // Prevent panning (and the star field's own independent pan) past the edge of the generated star field.
    double star_field_half_width = camera_state.max_display_width / (2.0 * MIN_ZOOM);
    double star_field_half_height = camera_state.max_display_height / (2.0 * MIN_ZOOM);
    DynamoEngine::Vector2D star_field_center = {SCREEN_WIDTH / 2.0, SCREEN_HEIGHT / 2.0};

    double view_half_width = (camera_state.window_width / 2.0) / camera_state.zoom;
    double view_half_height = (camera_state.window_height / 2.0) / camera_state.zoom;

    double slack_x = std::max(0.0, star_field_half_width - view_half_width);
    double slack_y = std::max(0.0, star_field_half_height - view_half_height);

    auto clamp_offset_to_star_field = [&](DynamoEngine::Vector2D& offset_to_clamp)
    {
        DynamoEngine::Vector2D view_center_world = {view_half_width - offset_to_clamp.x_val,
                                                    view_half_height - offset_to_clamp.y_val};

        view_center_world.x_val =
            std::clamp(view_center_world.x_val, star_field_center.x_val - slack_x, star_field_center.x_val + slack_x);
        view_center_world.y_val =
            std::clamp(view_center_world.y_val, star_field_center.y_val - slack_y, star_field_center.y_val + slack_y);

        offset_to_clamp.x_val = view_half_width - view_center_world.x_val;
        offset_to_clamp.y_val = view_half_height - view_center_world.y_val;
    };

    clamp_offset_to_star_field(camera_state.offset);
    clamp_offset_to_star_field(camera_state.twinkling_star_offset);
}

} // namespace

InputSystem::InputSystem()
{
    // Initialize input system variables if needed
    SDL_InitSubSystem(SDL_INIT_EVENTS);
}

InputSystem::~InputSystem() {}

// --- SYSTEM-LEVEL METHOD --- //
void InputSystem::processSystemInputFrame(GameState& game_state, UIState& ui_state, DynamoEngine::SceneManager& scenes,
                                          PlayerInput& player_input, float frame_seconds)
{
    m_input.beginInputFrame();
    m_frame_events.clear(); // Clear the last frame's events
    m_intake.pollEvents(m_frame_events);

    for (const DynamoEngine::InputEvent& event : m_frame_events)
    {
        trackDragAnchor(event);
        m_input.applyInputEvent(event);

        if (event.type == DynamoEngine::InputEventType::WindowResize)
        {
            CameraState& camera_state = game_state.getCameraStateMutable();
            camera_state.window_width = static_cast<float>(event.window_width);
            camera_state.window_height = static_cast<float>(event.window_height);
        }
    }
    if (m_input.quitRequested())
    {
        requestShutdown(game_state);
        return;
    }
    // The current scene's UI gets the mouse first; the game only acts on what the UI didn't take
    DynamoEngine::Scene& scene = scenes.currentScene();
    const DynamoEngine::UIInputResult ui_result = updateSceneUI(scene.ui(), game_state.getCameraState(), frame_seconds);

    // One-frame requests last exactly one frame: start every frame with none
    player_input.commands = PlayerCommands();

    const bool runs_simulation = scene.settings().runs_simulation;
    if (runs_simulation)
    {
        updateCamera(game_state);
        translateGameInputs(game_state, ui_state, scenes, player_input, ui_result.pointer_captured);
    }
    else
    {
        translateMenuInputs(scenes);
    }

    // A see-through preview follows the cursor while a spawn button is held down on the world (not on the UI):
    // left = the planet, right = the cluster's particles. If both are held, left wins (one preview at a time).
    const bool on_world = runs_simulation && !ui_result.pointer_captured;
    const bool aiming_planet = on_world && m_input.isMouseButtonDown(DynamoEngine::MouseButton::Left);
    const bool aiming_cluster =
        on_world && m_input.isMouseButtonDown(DynamoEngine::MouseButton::Right) && !aiming_planet;
    player_input.is_aiming_spawn = aiming_planet || aiming_cluster;
    if (player_input.is_aiming_spawn)
    {
        player_input.spawn_aim = aimSpawn(aiming_cluster, m_input.isShiftDown(), game_state.getCameraState());
    }
}

SpawnProperties InputSystem::aimSpawn(bool is_cluster, bool with_velocity, const CameraState& camera_state) const
{
    // Screen positions as doubles first: the world math is done in doubles
    const DynamoEngine::Vector2D cursor = m_input.mousePosition();
    const DynamoEngine::Vector2D drag_start = m_mouse_drag_anchor;

    SpawnProperties aim;
    aim.is_cluster = is_cluster;
    aim.with_velocity = with_velocity;
    if (with_velocity)
    {
        aim.world_position = ScreenToWorldCoordinates(drag_start, camera_state);
        aim.world_velocity = (cursor - drag_start) / camera_state.zoom; // the drag, in world units per second
    }
    else
    {
        aim.world_position = ScreenToWorldCoordinates(cursor, camera_state);
    }
    return aim;
}

DynamoEngine::UIInputResult InputSystem::updateSceneUI(DynamoEngine::UIRoot& ui, const CameraState& camera_state,
                                                       float frame_seconds)
{
    ui.setWindowSize({camera_state.window_width, camera_state.window_height});
    ui.updateElements(frame_seconds); // lays everything out for this window size, then e.g. refreshes the FPS text
    return ui.processInput(m_input);
}

void InputSystem::trackDragAnchor(const DynamoEngine::InputEvent& event)
{
    using DynamoEngine::InputEventType;
    using DynamoEngine::MouseButton;
    bool is_first_button_of_drag = event.type == InputEventType::MouseButtonDown && !m_input.isAnyMouseButtonDown();
    bool is_shift_press = event.type == InputEventType::KeyDown && !event.is_repeat &&
                          (event.key == SDL_SCANCODE_LSHIFT || event.key == SDL_SCANCODE_RSHIFT);
    bool is_creation_drag_active =
        m_input.isMouseButtonDown(MouseButton::Left) || m_input.isMouseButtonDown(MouseButton::Right);

    if (is_first_button_of_drag)
    {
        m_mouse_drag_anchor = event.mouse_position;
    }
    else if (is_shift_press && is_creation_drag_active)
    {
        m_mouse_drag_anchor = m_input.mousePosition(); // pressing shift mid-drag re-anchors the velocity arrow
    }
}
// --- CAMERA UPDATE --- //
void InputSystem::updateCamera(GameState& game_state)
{
    using DynamoEngine::MouseButton;
    CameraState& camera_state = game_state.getCameraStateMutable();
    const DynamoEngine::Vector2D mouse_position = m_input.mousePosition(); // Camera math uses doubles
    float scroll = m_input.mouseScrollDeltaThisFrame();

    zoomAroundCursor(camera_state, scroll, mouse_position);

    // Middle-mouse pan: the engine already sums this frame's motion, so no "previous position" bookkeeping
    if (m_input.isMouseButtonDown(MouseButton::Middle))
    {
        panCamera(camera_state, m_input.mousePositionDeltaThisFrame());
    }

    clampCameraToStarField(camera_state);
}

void InputSystem::translateMenuInputs(DynamoEngine::SceneManager& scenes)
{
    if (m_input.wasKeyPressed(SDL_SCANCODE_ESCAPE))
    {
        scenes.requestSwitch(TransferScene::Game);
    }
}

void InputSystem::translateGameInputs(GameState& game_state, UIState& ui_state, DynamoEngine::SceneManager& scenes,
                                      PlayerInput& player_input, bool pointer_captured)
{
    using DynamoEngine::MouseButton;

    PlayerCommands& commands = player_input.commands;

    // Delete / Backspace: undo the newest spawn; with Shift: clear everything. Once per tap.
    if (m_input.wasKeyPressed(SDL_SCANCODE_BACKSPACE) || m_input.wasKeyPressed(SDL_SCANCODE_DELETE))
    {
        if (m_input.isShiftDown())
        {
            commands.clear_all_bodies = true;
        }
        else
        {
            commands.undo_newest_spawn = true;
        }
        return;
    }
    if (m_input.wasKeyPressed(SDL_SCANCODE_ESCAPE))
    {
        scenes.requestSwitch(TransferScene::Pause);
        return;
    }

    // Starship controls: each pair counts only when exactly one of its two keys is held
    bool w_down = m_input.isKeyDown(SDL_SCANCODE_W);
    bool s_down = m_input.isKeyDown(SDL_SCANCODE_S);
    bool a_down = m_input.isKeyDown(SDL_SCANCODE_A);
    bool d_down = m_input.isKeyDown(SDL_SCANCODE_D);
    player_input.ship.thrust = (w_down != s_down) ? (w_down ? 1 : -1) : 0;
    player_input.ship.rotation = (a_down != d_down) ? (a_down ? 1 : -1) : 0;

    // Spawning happens on release, on the world (not on the UI): left = planet, right = particle cluster (left wins
    // if both). Shift gives it an initial velocity.
    bool left_released = m_input.wasMouseButtonReleased(MouseButton::Left);
    bool right_released = m_input.wasMouseButtonReleased(MouseButton::Right);
    if (!pointer_captured && (left_released || right_released))
    {
        commands.spawn = true;
        commands.spawn_aim = aimSpawn(!left_released, m_input.isShiftDown(), game_state.getCameraState());
    }

    // F3: show or hide the debug overlay (hitboxes)
    if (m_input.wasKeyPressed(SDL_SCANCODE_F3))
    {
        ui_state.setRenderDebug(!ui_state.getRenderDebug());
    }
    // F4: show or hide the profiler (how long each system takes: average and 1% low)
    if (m_input.wasKeyPressed(SDL_SCANCODE_F4))
    {
        ui_state.setShowProfiler(!ui_state.getShowProfiler());
    }
    // P: planning mode on/off (time stops; leaving it resumes at the speed slider's value)
    if (m_input.wasKeyPressed(SDL_SCANCODE_P))
    {
        ui_state.setIsPlanning(!ui_state.getIsPlanning());
    }
    // I, only while planning: trace the predicted trajectories, or clear them (PhysicsSystem::updateTrajectories)
    if (ui_state.getIsPlanning() && m_input.wasKeyPressed(SDL_SCANCODE_I))
    {
        commands.toggle_trajectories = true;
    }
    // Update the visor if necessary
    updateVisor(game_state.getCameraStateMutable());
}

void InputSystem::updateVisor(CameraState& camera_state)
{
    if (m_input.wasKeyPressed(SDL_SCANCODE_TAB))
    {
        camera_state.visor_view = nextVisorView(camera_state.visor_view);
    }
}
// --------- CLEANUP HELPER METHOD --------- //

void InputSystem::cleanUp()
{
    // Any necessary cleanup code for the input system
}

// --------- ADDITIONAL METHODS --------- //
