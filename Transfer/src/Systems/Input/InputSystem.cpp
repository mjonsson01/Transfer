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
        updateCamera(game_state, frame_seconds);
        translateGameInputs(game_state, ui_state, scenes, player_input, ui_result.pointer_captured);
    }
    else
    {
        translateMenuInputs(ui_state, scenes);
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
void InputSystem::updateCamera(GameState& game_state, float frame_seconds)
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

    // Arrow keys pan too (trackpads have no middle button). A drag moves the world WITH the mouse, so looking right
    // means moving the world left: the screen delta is the opposite of the arrow's direction.
    DynamoEngine::Vector2D arrow_direction = {0.0, 0.0};
    if (m_input.isKeyDown(SDL_SCANCODE_RIGHT))
        arrow_direction.x_val += 1.0;
    if (m_input.isKeyDown(SDL_SCANCODE_LEFT))
        arrow_direction.x_val -= 1.0;
    if (m_input.isKeyDown(SDL_SCANCODE_DOWN))
        arrow_direction.y_val += 1.0;
    if (m_input.isKeyDown(SDL_SCANCODE_UP))
        arrow_direction.y_val -= 1.0;
    if (arrow_direction.x_val != 0.0 || arrow_direction.y_val != 0.0)
    {
        // The speed is in UI points per second: it grows with the window like the UI does, so panning across the
        // screen takes the same time on a small window and a big one
        const double window_scale =
            DynamoEngine::UIRoot::fitScale({camera_state.window_width, camera_state.window_height});
        const double pan_speed = CAMERA_KEY_PAN_SPEED * window_scale;
        panCamera(camera_state, arrow_direction * (-pan_speed * frame_seconds));
    }

    const UniverseState& universe_state = game_state.getUniverseState();
    camera_state.debug_view_push_back =
        keepViewInsideWorld(camera_state, universe_state.world_min_corner, universe_state.world_max_corner);
}

void InputSystem::translateMenuInputs(const UIState& ui_state, DynamoEngine::SceneManager& scenes)
{
    if (m_input.wasKeyPressed(SDL_SCANCODE_ESCAPE))
    {
        // From the pause menu: back to the scene that was paused. From the start menu: into the game.
        bool is_paused = scenes.currentSceneId() == TransferScene::Pause;
        scenes.requestSwitch(is_paused ? ui_state.getResumeScene() : TransferScene::Game);
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
        ui_state.setResumeScene(scenes.currentSceneId()); // the pause menu returns here
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
    // F4: the frame profiler (each system); F5: the physics profiler (each part of the tick). The other key switches
    // over, the same key again hides it.
    if (m_input.wasKeyPressed(SDL_SCANCODE_F4))
    {
        bool showing = ui_state.getProfilerOverlay() == ProfilerOverlay::Frame;
        ui_state.setProfilerOverlay(showing ? ProfilerOverlay::Off : ProfilerOverlay::Frame);
    }
    if (m_input.wasKeyPressed(SDL_SCANCODE_F5))
    {
        bool showing = ui_state.getProfilerOverlay() == ProfilerOverlay::Physics;
        ui_state.setProfilerOverlay(showing ? ProfilerOverlay::Off : ProfilerOverlay::Physics);
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
