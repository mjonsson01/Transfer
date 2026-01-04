// File: Transfer/src/Systems/InputSystem.cpp

#include "Systems/InputSystem.h"

InputSystem::InputSystem()
{
    // Initialize input system variables if needed
}

InputSystem::~InputSystem()
{
    // Cleanup if necessary
}



// --------- SYSTEM-LEVEL METHOD --------- //

void InputSystem::ProcessSystemInputFrame(GameState& gameState, UIState& UIState)
{
    // Process input events and update the Game State, UI State, and internal input gameState accordingly
    SDL_Event event;
    while (SDL_PollEvent(&event)) {

        switch (event.type) {
            // Handle different event types here
            case SDL_EVENT_QUIT:
                gameState.SetPlaying(false);
                break;
            case SDL_EVENT_MOUSE_MOTION:
            {
                auto& updated_input_state = UIState.getMutableInputState();
                updated_input_state.mouseCurrPosition = {event.motion.x, event.motion.y};
                break;
            }
            case SDL_EVENT_MOUSE_BUTTON_DOWN:
            {
                if (event.button.button == SDL_BUTTON_LEFT) {
                    auto& updated_input_state = UIState.getMutableInputState();
                    updated_input_state.isHoldingLeftMouseButton = true;
                    updated_input_state.dirty = true;
                    updated_input_state.isCreatingMacro = true;
                    // updated_input_state.mouseCurrPosition = {event.button.x, event.button.y};
                    updated_input_state.selectedMass = MAX_MASS/10.0;
                    updated_input_state.selectedRadius = 50.0;
                    break;
                }
                if (event.button.button == SDL_BUTTON_RIGHT) {
                    auto& updated_input_state = UIState.getMutableInputState();
                    updated_input_state.isHoldingRightMouseButton = true;
                    updated_input_state.mouseCurrPosition = {event.button.x, event.button.y};
                    break;
                }
                
            }
            case SDL_EVENT_MOUSE_BUTTON_UP:
            {
                auto& updated_input_state = UIState.getMutableInputState();
                updated_input_state.isHoldingRightMouseButton = false;
                updated_input_state.isHoldingLeftMouseButton = false;
                updated_input_state.spawnAccumulator = 0.0;
                break;
            }
            case SDL_EVENT_KEY_DOWN:
                if (event.key.scancode == SDL_SCANCODE_SPACE){
                    auto& updated_input_state = UIState.getMutableInputState();
                    updated_input_state.dirty = true;
                    updated_input_state.isCreatingMacro = true;
                    updated_input_state.isCreatingStatic = true;
                    updated_input_state.selectedMass = MAX_MASS;
                    updated_input_state.selectedRadius = 50.0;
                    break;
                }
                else if (event.key.scancode == SDL_SCANCODE_TAB)
                {
                    gameState.invertToggleSlow();
                }
                else if (event.key.scancode == SDL_SCANCODE_P)
                {
                    UIState.getMutableInputState().togglePhysicsPause();
                }
                else if (event.key.scancode == SDL_SCANCODE_BACKSPACE || event.key.scancode == SDL_SCANCODE_DELETE)
                {
                    UIState.getMutableInputState().clearAllBodies();
                }
                // else if (event.key.scancode == SDL_SCANCODE_T)
                // {
                //     auto& updated_input_state = UIState.getMutableInputState();
                //     updated_input_state.dirty = true;
                //     updated_input_state.isCreatingParticleCluster = true ;
                //     updated_input_state.selectedMass = MAX_MASS/100000.0;
                //     updated_input_state.selectedRadius = 100.0;
                //     break;
                // }


        }
        
    }

    auto& updated_input_state = UIState.getMutableInputState();
    if (updated_input_state.isHoldingRightMouseButton)
    {
        updated_input_state.spawnAccumulator += FRAME_DELAY_MS;

        while (updated_input_state.spawnAccumulator >= SPAWN_DELAY_MS)
        {
            updated_input_state.spawnAccumulator -= SPAWN_DELAY_MS;

            updated_input_state.dirty = true;
            updated_input_state.isCreatingParticle = true;
            updated_input_state.selectedMass = MAX_MASS/100000.0;
            updated_input_state.selectedRadius = 1.0;
        }
    }
    else
    {
        updated_input_state.spawnAccumulator = 0.0;
    }

    // Additional input handling logic to be implemented later
}


// --------- CLEANUP HELPER METHOD --------- //

void InputSystem::CleanUp()
{
    // Any necessary cleanup code for the input system
}

// --------- ADDITIONAL METHODS --------- //
// void InputSystem::handleMassSliderInput(SDL_Event& event, UIState& UIState)
// {
//     // Implementation for handling mass slider input events
// }



