// File: Transfer/src/Systems/InputSystem.cpp

#include "Systems/InputSystem.hpp"

InputSystem::InputSystem()
{
    // Initialize input system variables if needed
    SDL_InitSubSystem(SDL_INIT_EVENTS);
}

InputSystem::~InputSystem() {}

// --------- SYSTEM-LEVEL METHOD --------- //

void InputSystem::ProcessSystemInputFrame(GameState& gameState, UIState& UIState)
{

    transferInputs.resetJustPressed();
    UIState.getMutableInputState().resetTransientFlags(); // clean the input state before polling for new events.
    SDL_Event event;
    // int eventCount = 0;
    while (SDL_PollEvent(&event))
    {
        // eventCount++;
        if (event.type == SDL_EVENT_QUIT)
        {
            gameState.SetPlaying(false);
            gameState.setIsShuttingDownAudioSystem(true);
        }
        else if (event.type == SDL_EVENT_WINDOW_RESIZED)
        {
            CameraState& camera_state_local = gameState.getCameraStateMutable();
            camera_state_local.windowWidth = (float)event.window.data1;
            camera_state_local.windowHeight = (float)event.window.data2;
        }
        else
        {
            // Will ensure the event is not yet consumed by the UI.
            UIState.getMutableInputState().UIInputConsumed = false;
            // First check if in start menu. If so, route input to start menu behaviors
            SceneIdentifier current_scene = UIState.getCurrentSceneID();

            if (current_scene == SceneIdentifier::GAME_SCENE)
            {
                routeSDL_EventInputInGame(&event); // writes to internal member transferInputs;
            }
            else
            {
                routeSDL_EventInputInMenu(&event);
            }
        }
    }

    SceneIdentifier current_scene = UIState.getCurrentSceneID();
    if (current_scene == SceneIdentifier::GAME_SCENE)
    {
        CameraState& camera_state = gameState.getCameraStateMutable();

        if (!firstWithinEpsilonOfSecond(transferInputs.pendingScrollData, 0.0f))
        {
            Vector2D worldUnderCursor = ScreenToWorldCoordinates(transferInputs.mouseCurrPosition, camera_state);
            Vector2D starWorldUnderCursor =
                transferInputs.mouseCurrPosition / camera_state.zoom - camera_state.twinklingStarOffset;

            camera_state.zoom *= std::pow(1.1, transferInputs.pendingScrollData);
            camera_state.zoom = std::clamp(camera_state.zoom, MIN_ZOOM, MAX_ZOOM);

            camera_state.offset = transferInputs.mouseCurrPosition / camera_state.zoom - worldUnderCursor;
            camera_state.twinklingStarOffset =
                transferInputs.mouseCurrPosition / camera_state.zoom - starWorldUnderCursor;
            transferInputs.pendingScrollData = 0.0f;
        }

        if (transferInputs.middleMouseJustPressed)
        {
            transferInputs.previousMiddleDragPosition = transferInputs.mouseCurrPosition;
        }
        else if (transferInputs.middleMousePressed)
        {
            Vector2D dragDelta = transferInputs.mouseCurrPosition - transferInputs.previousMiddleDragPosition;
            camera_state.offset += dragDelta / camera_state.zoom;
            camera_state.twinklingStarOffset += (dragDelta / camera_state.zoom) * STAR_PARALLAX_FACTOR;
            transferInputs.previousMiddleDragPosition = transferInputs.mouseCurrPosition;
        }
        // Prevent panning (and the star field's own independent pan) past the edge of
        // the generated star field.
        double starFieldHalfWidth = camera_state.maxDisplayWidth / (2.0 * MIN_ZOOM);
        double starFieldHalfHeight = camera_state.maxDisplayHeight / (2.0 * MIN_ZOOM);
        Vector2D starFieldCenter = {SCREEN_WIDTH / 2.0, SCREEN_HEIGHT / 2.0};

        double viewHalfWidth = (camera_state.windowWidth / 2.0) / camera_state.zoom;
        double viewHalfHeight = (camera_state.windowHeight / 2.0) / camera_state.zoom;

        double slackX = std::max(0.0, starFieldHalfWidth - viewHalfWidth);
        double slackY = std::max(0.0, starFieldHalfHeight - viewHalfHeight);

        auto clampOffsetToStarField = [&](Vector2D& offsetToClamp)
        {
            Vector2D viewCenterWorld = {viewHalfWidth - offsetToClamp.xVal, viewHalfHeight - offsetToClamp.yVal};

            viewCenterWorld.xVal =
                std::clamp(viewCenterWorld.xVal, starFieldCenter.xVal - slackX, starFieldCenter.xVal + slackX);
            viewCenterWorld.yVal =
                std::clamp(viewCenterWorld.yVal, starFieldCenter.yVal - slackY, starFieldCenter.yVal + slackY);

            offsetToClamp.xVal = viewHalfWidth - viewCenterWorld.xVal;
            offsetToClamp.yVal = viewHalfHeight - viewCenterWorld.yVal;
        };

        clampOffsetToStarField(camera_state.offset);
        clampOffsetToStarField(camera_state.twinklingStarOffset);
        translateAndPassTransferInputsOff(UIState);
    }
    else
    {
        translateAndPassMenuInputsOff(UIState);
    }
}

// Only need to worry about clicking
void InputSystem::routeSDL_EventInputInMenu(SDL_Event* e)
{
    switch (e->type)
    {
    case SDL_EVENT_MOUSE_BUTTON_DOWN:
        transferInputs.mouseCurrPosition = {e->button.x, e->button.y};

        switch (e->button.button)
        {
        case SDL_BUTTON_LEFT:
            transferInputs.leftMouseJustPressed = true;
            transferInputs.leftMousePressed = true;
            break;
        default:
            break;
        }
        break;
    case SDL_EVENT_MOUSE_BUTTON_UP:
        transferInputs.mouseCurrPosition = {e->button.x, e->button.y};

        switch (e->button.button)
        {
        case SDL_BUTTON_LEFT:
            transferInputs.leftMousePressed = false;
            transferInputs.leftMouseJustReleased = true;
            break;
        default:
            break;
        }
        break;
    case SDL_EVENT_KEY_DOWN:
        switch (e->key.scancode)
        {
        case SDL_SCANCODE_ESCAPE:
            if (e->key.repeat == 0)
            {
                transferInputs.escJustPressed = true;
            }
            transferInputs.escPressed = true;
            break;
        default:
            break;
        }
        break;
    case SDL_EVENT_KEY_UP:
        switch (e->key.scancode)
        {
        case SDL_SCANCODE_ESCAPE:
            transferInputs.escPressed = false;
            break;
        default:
            break;
        }
        break;
    default:
        break;
    }
}
void InputSystem::routeSDL_EventInputInGame(SDL_Event* e)
{
    switch (e->type)
    {
    case SDL_EVENT_MOUSE_BUTTON_DOWN:
        transferInputs.mouseCurrPosition = {e->button.x, e->button.y};

        // if no mouse buttons are currently being dragged, set the start of the drag location. otherwise just
        // continue tracking buttons (the start position will remain fixed)
        if (!transferInputs.leftMousePressed && !transferInputs.rightMousePressed && !transferInputs.middleMousePressed)
        {
            // This is a dragging event
            transferInputs.mouseDragStartPosition = transferInputs.mouseCurrPosition;
            transferInputs.isDragging = true;
        }
        switch (e->button.button)
        {
        case SDL_BUTTON_LEFT:
            transferInputs.leftMousePressed = true;
            transferInputs.leftMouseJustPressed = true;
            break;
        case SDL_BUTTON_RIGHT:
            transferInputs.rightMousePressed = true;
            break;
        case SDL_BUTTON_MIDDLE:
            transferInputs.middleMousePressed = true;
            transferInputs.middleMouseJustPressed = true;
            break;
        }
        break;
    case SDL_EVENT_MOUSE_BUTTON_UP:
        transferInputs.mouseCurrPosition = {e->button.x, e->button.y};
        transferInputs.isDragging = false;
        switch (e->button.button)
        {
        case SDL_BUTTON_LEFT:
            transferInputs.leftMousePressed = false;
            transferInputs.leftMouseJustReleased = true;
            break;
        case SDL_BUTTON_RIGHT:
            transferInputs.rightMousePressed = false;
            break;
        case SDL_BUTTON_MIDDLE:
            transferInputs.middleMousePressed = false;
            break;
        }
        break;

    case SDL_EVENT_MOUSE_MOTION:
        transferInputs.mouseCurrPosition = {e->motion.x, e->motion.y};
        break;

    case SDL_EVENT_KEY_DOWN:
        switch (e->key.scancode)
        {
        case SDL_SCANCODE_W:
            transferInputs.wPressed = true;
            break;
        case SDL_SCANCODE_A:
            transferInputs.aPressed = true;
            break;
        case SDL_SCANCODE_S:
            transferInputs.sPressed = true;
            break;
        case SDL_SCANCODE_D:
            transferInputs.dPressed = true;
            break;
        //  Left Alt Preferred
        case SDL_SCANCODE_LALT:
            transferInputs.altPressed = true;
            break;
        // Left Shift Preferred
        case SDL_SCANCODE_LSHIFT:
            // update drag start position if repressing shift
            if (e->key.repeat == 1)
                break;
            if (transferInputs.leftMousePressed || transferInputs.rightMousePressed)
            {
                if (transferInputs.isDragging)
                {
                    transferInputs.mouseDragStartPosition = transferInputs.mouseCurrPosition;
                }
            }
            transferInputs.shiftPressed = true;
            break;
        case SDL_SCANCODE_SPACE:
            transferInputs.spacePressed = true;
            break;
        case SDL_SCANCODE_ESCAPE:
            if (e->key.repeat == 0)
            {
                transferInputs.escJustPressed = true;
            }
            transferInputs.escPressed = true;
            // pause menu flag set?
            break;
        case SDL_SCANCODE_BACKSPACE:
        case SDL_SCANCODE_DELETE:
            transferInputs.clearParticlesPressed = true;
            break;
        default:
            break;
        }
        break;
    case SDL_EVENT_KEY_UP:
        switch (e->key.scancode)
        {
        case SDL_SCANCODE_W:
            transferInputs.wPressed = false;
            break;
        case SDL_SCANCODE_A:
            transferInputs.aPressed = false;
            break;
        case SDL_SCANCODE_S:
            transferInputs.sPressed = false;
            break;
        case SDL_SCANCODE_D:
            transferInputs.dPressed = false;
            break;
        // Left Alt Preferred
        case SDL_SCANCODE_LALT:
            transferInputs.altPressed = false;
            break;
        // Left Shift Preferred
        case SDL_SCANCODE_LSHIFT:
            transferInputs.shiftPressed = false;
            break;
        case SDL_SCANCODE_SPACE:
            transferInputs.spacePressed = false;
            break;
        case SDL_SCANCODE_ESCAPE:
            transferInputs.escPressed = false;
            break;
        case SDL_SCANCODE_BACKSPACE:
        case SDL_SCANCODE_DELETE:
            transferInputs.clearParticlesPressed = false;
            break;
        default:
            break;
        }
        break;
    case SDL_EVENT_MOUSE_WHEEL:
    {
        float scroll_y = e->wheel.y;
        if (e->wheel.direction == SDL_MOUSEWHEEL_FLIPPED)
        {
            scroll_y *= -1.0f;
        }
        transferInputs.pendingScrollData += scroll_y;
        break;
    }
    default:
        // no fall through behavior necessary
        break;
    }
}

void InputSystem::translateAndPassMenuInputsOff(UIState& UIState)
{
    InputState& updated_input_state = UIState.getMutableInputState();
    updated_input_state.mouseCurrPosition = transferInputs.mouseCurrPosition;
    updated_input_state.isDragging = transferInputs.isDragging;
    updated_input_state.mouseDragStartPosition = transferInputs.mouseDragStartPosition;
    updated_input_state.isClickingLeftMouseButton = transferInputs.leftMousePressed;
    updated_input_state.isClickingRightMouseButton = transferInputs.rightMousePressed;
    updated_input_state.isClickingMiddleMouseButton = transferInputs.middleMousePressed;
    updated_input_state.isPressingShift = transferInputs.shiftPressed;

    updated_input_state.leftMouseButtonJustPressed = transferInputs.leftMouseJustPressed;
    updated_input_state.leftMouseButtonJustReleased = transferInputs.leftMouseJustReleased;
    if (transferInputs.escJustPressed)
    {
        UIState.setCurrentScene(SceneIdentifier::GAME_SCENE);
        transferInputs.resetAllInputsForSceneChange();
        updated_input_state.resetFlagsForSceneChange();
        return;
    }
}
void InputSystem::translateAndPassTransferInputsOff(UIState& UIState)
{
    // Check for clear all particle orders
    InputState& updated_input_state = UIState.getMutableInputState();
    if (transferInputs.clearParticlesPressed)
    {
        updated_input_state.clearAllBodies();
        transferInputs.resetAllKeyPressedVars();
        transferInputs.resetAllMousePressedVars();
        transferInputs.resetJustPressed();
        updated_input_state.resetTransientFlags();
        return;
    }
    if (transferInputs.escJustPressed)
    {
        UIState.setCurrentScene(SceneIdentifier::PAUSE_SCENE);
        transferInputs.resetAllInputsForSceneChange();
        updated_input_state.resetFlagsForSceneChange();
        return;
    }
    // Always pass off these

    updated_input_state.mouseCurrPosition = transferInputs.mouseCurrPosition;
    updated_input_state.isDragging = transferInputs.isDragging;
    updated_input_state.mouseDragStartPosition = transferInputs.mouseDragStartPosition;
    updated_input_state.isClickingLeftMouseButton = transferInputs.leftMousePressed;
    updated_input_state.isClickingRightMouseButton = transferInputs.rightMousePressed;
    updated_input_state.isClickingMiddleMouseButton = transferInputs.middleMousePressed;
    updated_input_state.isPressingShift = transferInputs.shiftPressed;

    updated_input_state.leftMouseButtonJustPressed = transferInputs.leftMouseJustPressed;
    updated_input_state.leftMouseButtonJustReleased = transferInputs.leftMouseJustReleased;
    if (transferInputs.leftMouseJustReleased)
    {
        updated_input_state.isCreatingCollidable = true;
        if (transferInputs.shiftPressed)
        {
            updated_input_state.isCreatingWithInitialVelocity = true;
        }
        updated_input_state.isCreatingMacro = true;
    }
}

// --------- CLEANUP HELPER METHOD --------- //

void InputSystem::CleanUp()
{
    // Any necessary cleanup code for the input system
}

// --------- ADDITIONAL METHODS --------- //
