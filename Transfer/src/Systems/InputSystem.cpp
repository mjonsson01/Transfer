// File: Transfer/src/Systems/InputSystem.cpp

#include "Systems/InputSystem.h"

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
    while (SDL_PollEvent(&event))
    {
        if (event.type == SDL_EVENT_QUIT)
        {
            gameState.SetPlaying(false);
            gameState.setIsShuttingDownAudioSystem(true);
        }
        else
        {
            // Will ensure the event is not yet consumed by the UI.
            UIState.getMutableInputState().UIInputConsumed = false;
            // First check if in start menu. If so, route input to start menu behaviors
            if (UIState.getStartMenuActive())
            {
                // Handle start menu specific input
                break;
            }
            // If not in start menu, double check we are not in a pause menu. If so, route to pause menu behaviors
            if (UIState.getPauseMenuActive())
            {
                // Handle pause menu specific input
                break;
            }
            // Check if we are in a level editor. Not implemented yet so will never enter for now.
            if (UIState.getLevelEditorScene())
            {
                // Handle level editor specific input
                break;
            }
            // Check if we are in a level. If so, route commands accordingly.
            if (UIState.getLevelScene())
            {
                // Handle level scene specific input
                routeSDL_EventInputInGame(&event); // writes to internal member transferInputs;

                // UI Only interactible through mouse clicks, but pass through the transfer game inputs so that the
                // other systems have knowledge of the desired areas.

                // pass off data to the UIState InputState sub data. UISystem may or may not consume this input. If it
                // does not, the UISystem will flip flags as necessary to make sure the input is consumed
                translateAndPassTransferInputsOff(UIState);

                break;
            }
        }
    }
    // std::cout<< "Transfer Inputs -> " << transferInputs << std::endl;
    // Now that our event is routed to an InputRoutedEvent (and we haven't quit), we can pass off our details to the UI
    // for checking first, then to the rest of the game if
}

void InputSystem::routeSDL_EventInputInGame(SDL_Event* e)
{
    switch (e->type)
    {
    case SDL_EVENT_MOUSE_BUTTON_DOWN:
        // if no mouse buttons are currently being dragged, set the start of the drag location. otherwise just continue
        // tracking buttons (the start position will remain fixed)
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
            break;
        }
        break;
    case SDL_EVENT_MOUSE_BUTTON_UP:
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
            transferInputs.shiftPressed = true;
            break;
        case SDL_SCANCODE_SPACE:
            transferInputs.spacePressed = true;
            break;
        case SDL_SCANCODE_ESCAPE:
            transferInputs.escPressed = false;
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
    default:
        // no fall through behavior necessary
        break;
    }
}

void InputSystem::translateAndPassTransferInputsOff(UIState& UIState)
{
    // Check for clear all particle orders
    if (transferInputs.clearParticlesPressed)
    {
        UIState.getMutableInputState().clearAllBodies();
        return;
    }
    // Always pass off these
    InputState& updated_input_state = UIState.getMutableInputState();
    updated_input_state.mouseCurrPosition = transferInputs.mouseCurrPosition;
    updated_input_state.isDragging = transferInputs.isDragging;
    updated_input_state.mouseDragStartPosition = transferInputs.mouseDragStartPosition;
    updated_input_state.isClickingLeftMouseButton = transferInputs.leftMousePressed;
    updated_input_state.isClickingRightMouseButton = transferInputs.rightMousePressed;
    updated_input_state.isClickingMiddleMouseButton = transferInputs.middleMousePressed;
    updated_input_state.isPressingShift = transferInputs.shiftPressed;

    updated_input_state.leftMouseButtonJustPressed = transferInputs.leftMouseJustPressed;

    if (transferInputs.leftMouseJustReleased)
    {
        if (transferInputs.shiftPressed)
        {
            updated_input_state.isCreatingWithInitialVelocity = true;
        }
        updated_input_state.isCreatingMacro = true;
        updated_input_state.isCreatingCollidable = true;
    }
}

// --------- CLEANUP HELPER METHOD --------- //

void InputSystem::CleanUp()
{
    // Any necessary cleanup code for the input system
}

// --------- ADDITIONAL METHODS --------- //
