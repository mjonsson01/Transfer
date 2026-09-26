// File: Transfer/src/Systems/InputSystem.cpp

// Custom Includes
#include "Systems/InputSystem.hpp"

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
bool isSimulationScene(SceneIdentifier scene_id)
{
    return (scene_id == SceneIdentifier::GAME_SCENE || scene_id == SceneIdentifier::TEST_VISUAL_SCENE);
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
void InputSystem::processSystemInputFrame(GameState& game_state, UIState& ui_state)
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
    // Pass Engine State off into Transfer's meaning for current scene
    ui_state.getMutableDEPRECATED_InputState().resetTransientFlags();

    if (isSimulationScene(ui_state.getCurrentSceneID()))
    {
        updateCamera(game_state);
        translateGameInputs(ui_state);
    }
    else
    {
        translateMenuInputs(ui_state);
    }
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

void InputSystem::copySharedPointerState(DEPRECATED_InputState& legacy_state)
{
    using DynamoEngine::MouseButton;

    legacy_state.mouseCurrPosition = m_input.mousePosition();
    legacy_state.isDragging = m_input.isAnyMouseButtonDown();
    legacy_state.mouseDragStartPosition = m_mouse_drag_anchor;
    legacy_state.isClickingLeftMouseButton = m_input.isMouseButtonDown(MouseButton::Left);
    legacy_state.isClickingRightMouseButton = m_input.isMouseButtonDown(MouseButton::Right);
    legacy_state.isClickingMiddleMouseButton = m_input.isMouseButtonDown(MouseButton::Middle);
    legacy_state.isPressingShift = m_input.isShiftDown(); // either Shift now works, not just left
    legacy_state.leftMouseButtonJustPressed = m_input.wasMouseButtonPressed(MouseButton::Left);
    legacy_state.leftMouseButtonJustReleased = m_input.wasMouseButtonReleased(MouseButton::Left);
}

void InputSystem::translateMenuInputs(UIState& ui_state)
{
    DEPRECATED_InputState& legacy_state = ui_state.getMutableDEPRECATED_InputState();
    copySharedPointerState(legacy_state);

    if (m_input.wasKeyPressed(SDL_SCANCODE_ESCAPE))
    {
        ui_state.setCurrentScene(SceneIdentifier::GAME_SCENE);
        legacy_state.resetFlagsForSceneChange();
    }
}

void InputSystem::translateGameInputs(UIState& ui_state)
{
    using DynamoEngine::MouseButton;

    DEPRECATED_InputState& legacy_state = ui_state.getMutableDEPRECATED_InputState();

    // Clear the screen once per tap (was: every frame while the key was held)
    if (m_input.wasKeyPressed(SDL_SCANCODE_BACKSPACE) || m_input.wasKeyPressed(SDL_SCANCODE_DELETE))
    {
        legacy_state.clearAllBodies();
        legacy_state.resetTransientFlags();
        return;
    }
    if (m_input.wasKeyPressed(SDL_SCANCODE_ESCAPE))
    {
        ui_state.setCurrentScene(SceneIdentifier::PAUSE_SCENE);
        legacy_state.resetFlagsForSceneChange();
        return;
    }

    copySharedPointerState(legacy_state);

    // Starship controls
    bool w_down = m_input.isKeyDown(SDL_SCANCODE_W);
    bool s_down = m_input.isKeyDown(SDL_SCANCODE_S);
    bool a_down = m_input.isKeyDown(SDL_SCANCODE_A);
    bool d_down = m_input.isKeyDown(SDL_SCANCODE_D);
    legacy_state.isRequestingThrust = w_down != s_down; // exactly one of the pair (same as the old xor)
    legacy_state.positiveThrust = w_down;
    legacy_state.negativeThrust = s_down;
    legacy_state.isRequestingRotation = d_down != a_down;
    legacy_state.positiveRotation = a_down;
    legacy_state.negativeRotation = d_down;

    // Spawning happens on release: left = macro body, right = particle cluster. Shift adds initial velocity.
    if (m_input.wasMouseButtonReleased(MouseButton::Left))
    {
        legacy_state.isCreatingCollidable = true;
        legacy_state.isCreatingWithInitialVelocity = m_input.isShiftDown();
        legacy_state.isCreatingMacro = true;
    }
    if (m_input.wasMouseButtonReleased(MouseButton::Right))
    {
        legacy_state.isCreatingCollidable = true;
        legacy_state.isCreatingWithInitialVelocity = m_input.isShiftDown();
        legacy_state.isCreatingParticleCluster = true;
    }
}
// --------- CLEANUP HELPER METHOD --------- //

void InputSystem::cleanUp()
{
    // Any necessary cleanup code for the input system
}

// --------- ADDITIONAL METHODS --------- //
