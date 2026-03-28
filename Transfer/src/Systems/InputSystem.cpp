// File: Transfer/src/Systems/InputSystem.cpp

#include "Systems/InputSystem.h"

InputSystem::InputSystem()
{
    // Initialize input system variables if needed
    SDL_InitSubSystem(SDL_INIT_EVENTS);
}

InputSystem::~InputSystem()
{

}



// --------- SYSTEM-LEVEL METHOD --------- //

// void InputSystem::ProcessSystemInputFrame(GameState& gameState, UIState& UIState)
// {
//     // Process input events and update the Game State, UI State, and internal input gameState accordingly
//     SDL_Event event;
//     while (SDL_PollEvent(&event)) {
//         auto& updated_input_state = UIState.getMutableInputState();
//         auto& current_input_state = UIState.getInputState(); 

//         switch (event.type)
//         {
//             // Handle different event types here
//             case SDL_EVENT_QUIT:
//                 gameState.SetPlaying(false);
//                 gameState.setIsShuttingDownAudioSystem(true);
//                 break;
//             case SDL_EVENT_MOUSE_MOTION:
//             {
//                 updated_input_state.mouseCurrPosition = {event.motion.x, event.motion.y};
//                 // If dragging the mass knob, update the drag end position and knob location
//                 // if (updated_input_state.isDraggingMassKnob) {
//                 //     updated_input_state.mouseDragEndPosition = {event.motion.x, event.motion.y};
//                 //     UIState.setMassKnobRectPosition(updated_input_state.mouseDragEndPosition);
//                 // }
//                 break;
//             }
//             case SDL_EVENT_MOUSE_BUTTON_DOWN:
//             {

//                 // if held down and dragging inside the space of the knob rectangle location, capture the motion
//                 // and pull the knob rect along with it.
                
//                 // if (UIState.inRect(UIState.getMassKnobRect()))
//                 // {
//                 //     std::cout << "Dragging mass knob!" << std::endl;
//                 //     updated_input_state.isDraggingMassKnob = true;
//                 //     updated_input_state.mouseDragStartPosition = {event.motion.x, event.motion.y};
//                 //     UIState.setMassKnobRectPosition(updated_input_state.mouseDragStartPosition);
//                 //     break; // need to update position
//                 // }
//                 if (event.button.button == SDL_BUTTON_LEFT)
//                 {
//                     updated_input_state.isHoldingLeftMouseButton = true;
//                     updated_input_state.dirty = true;
//                     updated_input_state.isCreatingMacro = true;
//                     updated_input_state.selectedMass = MAX_MASS/10.0;
//                     updated_input_state.selectedRadius = 50.0;
//                     break;
//                 }
//                 if (event.button.button == SDL_BUTTON_RIGHT)
//                 {
//                     updated_input_state.isHoldingRightMouseButton = true;
//                     break;
//                 }
                
//             }
//             case SDL_EVENT_MOUSE_BUTTON_UP:
//             {
//                 auto& updated_input_state = UIState.getMutableInputState();
//                 updated_input_state.isHoldingRightMouseButton = false;
//                 updated_input_state.isHoldingLeftMouseButton = false;
//                 // updated_input_state.isDraggingMassKnob = false;
//                 updated_input_state.spawnAccumulator = 0.0;
//                 break;
//             }
//             case SDL_EVENT_KEY_DOWN:
//                 if (event.key.scancode == SDL_SCANCODE_SPACE)
//                 {
//                     auto& updated_input_state = UIState.getMutableInputState();
//                     updated_input_state.dirty = true;
//                     updated_input_state.isCreatingMacro = true;
//                     updated_input_state.isCreatingStatic = true;
//                     updated_input_state.selectedMass = MAX_MASS;
//                     updated_input_state.selectedRadius = 50.0;
//                     break;
//                 }
//                 else if (event.key.scancode == SDL_SCANCODE_TAB)
//                 {
//                     gameState.invertToggleSlow();
//                 }
//                 else if (event.key.scancode == SDL_SCANCODE_P)
//                 {
//                     UIState.getMutableInputState().togglePhysicsPause();
//                 }
//                 else if (event.key.scancode == SDL_SCANCODE_BACKSPACE || event.key.scancode == SDL_SCANCODE_DELETE)
//                 {
//                     UIState.getMutableInputState().clearAllBodies();
//                 }
//                 else if (event.key.scancode == SDL_SCANCODE_0)
//                 {
//                     gameState.invertPlayMusic();
//                     // std::cout<<"0 hit"<<std::endl;
//                 }
//                 else if (event.key.scancode == SDL_SCANCODE_MINUS)
//                 {
//                     UIState.invertUIElementsVisibility();
//                 }
//                 else if (event.key.scancode == SDL_SCANCODE_EQUALS)
//                 {
//                     gameState.invertToggleFast();
//                 }
//                 else if (event.key.scancode == SDL_SCANCODE_T)
//                 {
//                     auto& updated_input_state = UIState.getMutableInputState();
//                     updated_input_state.dirty = true;
//                     updated_input_state.isCreatingParticleCluster = true ;
//                     updated_input_state.selectedMass = MAX_MASS;
//                     updated_input_state.selectedRadius = 100.0;
//                     break;
//                 }


//         }
        
//     }

//     auto& updated_input_state = UIState.getMutableInputState();
//     if (updated_input_state.isHoldingRightMouseButton)
//     {
//         updated_input_state.spawnAccumulator += FRAME_DELAY_MS;

//         while (updated_input_state.spawnAccumulator >= SPAWN_DELAY_MS)
//         {
//             updated_input_state.spawnAccumulator -= SPAWN_DELAY_MS;

//             updated_input_state.dirty = true;
//             updated_input_state.isCreatingParticle = true;
//             updated_input_state.selectedMass = MAX_MASS/100000.0;
//             updated_input_state.selectedRadius = 1.0;
//         }
//     }
//     else
//     {
//         updated_input_state.spawnAccumulator = 0.0;
//     }

//     // Additional input handling logic to be implemented later
// }

void InputSystem::ProcessSystemInputFrame(GameState& gameState, UIState& UIState)
{
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
            if  (UIState.getStartMenuActive())
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
                routeSDL_EventInputInGame(&event); //writes to internal member transferInputs;

                // UI Only interactible through mouse clicks, but pass through the transfer game inputs so that the other systems have knowledge of the desired areas.

                // pass off data to the UIState InputState sub data. UISystem may or may not consume this input. If it does not, the UISystem will flip flags as necessary to make sure the input is consumed
                translateAndPassTransferInputsOff(UIState);

                break;
            }
        }
    }
    // std::cout<< "Transfer Inputs -> " << transferInputs << std::endl;
    // Now that our event is routed to an InputRoutedEvent (and we haven't quit), we can pass off our details to the UI for checking first, then to the rest of the game
    // if 
}

void InputSystem::routeSDL_EventInputInGame(SDL_Event* e)
{
    switch (e->type)
    {
        case SDL_EVENT_MOUSE_BUTTON_DOWN:
            // if no mouse buttons are currently being dragged, set the start of the drag location. otherwise just continue tracking buttons (the start position will remain fixed)
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
    updated_input_state.isHoldingLeftMouseButton = transferInputs.leftMousePressed;
    updated_input_state.isHoldingRightMouseButton = transferInputs.rightMousePressed;
    updated_input_state.isHoldingMiddleMouseButton = transferInputs.rightMousePressed;
    
    if (transferInputs.leftMousePressed)
    {
        updated_input_state.isCreatingMacro = true;
    //  updated_input_state.selectedMass = MAX_MASS/10.0;
        updated_input_state.selectedRadius = 50.0;
    }
    
}


// --------- CLEANUP HELPER METHOD --------- //

void InputSystem::CleanUp()
{
    // Any necessary cleanup code for the input system
}

// --------- ADDITIONAL METHODS --------- //

