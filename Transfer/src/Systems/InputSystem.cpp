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

void InputSystem::CleanUp()
{
    // Any necessary cleanup code for the input system
}
void InputSystem::ProcessSystemInputFrame(GameState& state, UIState& UIState)
{
    // Process input events and update the game state accordingly
    // (Implementation details would go here)

    SDL_Event event;
    while (SDL_PollEvent(&event)) {

        switch (event.type) {
            // Handle different event types here
            case SDL_EVENT_QUIT:
                state.SetPlaying(false);
                break;
            case SDL_EVENT_MOUSE_MOTION:
            {
                auto& newState = UIState.getMutableInputState();
                newState.mouseCurrPosition = {event.motion.x, event.motion.y};
                break;
            }
            case SDL_EVENT_MOUSE_BUTTON_DOWN:
            {
                if (event.button.button == SDL_BUTTON_LEFT) {
                    auto& newState = UIState.getMutableInputState();
                    newState.isHoldingLeftMouseButton = true;
                    newState.dirty = true;
                    newState.isCreatingPlanet = true;
                    newState.isCreatingDust = false;
                    newState.isCreatingCluster = false;
                    newState.isCreatingStatic = false;
                    newState.mouseCurrPosition = {event.button.x, event.button.y};
                    newState.selectedMass = MAX_MASS/10.0;
                    newState.selectedRadius = 50.0;
                    break;
                }
                if (event.button.button == SDL_BUTTON_RIGHT) {
                    // Right mouse button logic (if any) goes here
                    auto& newState = UIState.getMutableInputState();
                    newState.isHoldingRightMouseButton = true;
                    newState.mouseCurrPosition = {event.button.x, event.button.y};
                    break;
                }
                
            }
            case SDL_EVENT_MOUSE_BUTTON_UP:
            {
                // std::cout<<"mouse button up"<<std::endl;
                auto& newState = UIState.getMutableInputState();
                newState.isHoldingRightMouseButton = false;
                newState.isHoldingLeftMouseButton = false;
                newState.spawnAccumulator = 0.0;
                break;
            }
            case SDL_EVENT_KEY_DOWN:
                if (event.key.scancode == SDL_SCANCODE_SPACE){
                    auto& newState = UIState.getMutableInputState();
                    newState.dirty = true;
                    newState.isCreatingPlanet = true;
                    newState.isCreatingStatic = true;
                    // newState.mouseCurrPosition = {event.button.x, event.button.y};
                    // newState.mouseCurrPosition = { SCREEN_WIDTH / 2, 2 * SCREEN_HEIGHT / 3 };
                    newState.selectedMass = MAX_MASS/2.0;
                    newState.selectedRadius = 50.0;
                    break;
                }
                else if (event.key.scancode == SDL_SCANCODE_TAB)
                {
                    state.invertToggleSlow();
                }
                else if (event.key.scancode == SDL_SCANCODE_P)
                {
                    UIState.getMutableInputState().isPaused = !UIState.getMutableInputState().isPaused;
                }
                else if (event.key.scancode == SDL_SCANCODE_BACKSPACE || event.key.scancode == SDL_SCANCODE_DELETE)
                {
                    UIState.getMutableInputState().clearAll = true;
                }


        }
        
    }

    auto& inputState = UIState.getMutableInputState();
    double dt = FRAME_DELAY_MS;

    if (inputState.isHoldingRightMouseButton)
    {
        inputState.spawnAccumulator += dt;

        while (inputState.spawnAccumulator >= SPAWN_DELAY_MS)
        {
            inputState.spawnAccumulator -= SPAWN_DELAY_MS;

            inputState.dirty = true;
            inputState.isCreatingDust = true;
            inputState.isCreatingPlanet = false;
            inputState.selectedMass = MAX_MASS/1000.0;
            inputState.selectedRadius = 1.0;
        }
    }
    else
    {
        inputState.spawnAccumulator = 0.0;
    }

        // Additional input handling logic can be added here
}

// void InputSystem::createNewBody(SDL_Event& event, GameState& state)
// {
//     GravitationalBody newBody;

//     float radius = 15.0;
//     newBody.setRadius(radius);

//     double mass = MAX_MASS; // in kilograms
//     newBody.setMass(mass);

//     // Vector2D bodyVelocity = {100,0};
//     // newBody.setNetVelocity(bodyVelocity);

//     newBody.setCollisionEnabled(true);
    
//     Vector2D bodyPosition = { event.button.x, event.button.y }; // position measured in pixel xy cooriifnates measured from top-left corner.
//     newBody.setPosition(bodyPosition);

//     Vector2D prevPosition = bodyPosition; // previous pixel position
//     newBody.setPrevPosition(prevPosition);

//     state.addBody(newBody);

//     std::cout<<"New Body Created: "<<std::endl;
//     std::cout<<newBody<<std::endl;

// };

void InputSystem::handleMassSliderInput(SDL_Event& event, UIState& UIState)
{
    // Implementation for handling mass slider input events
    // (This function can be expanded based on specific requirements)
}


// void InputSystem::createNewGravitationalCluster(SDL_Event& event, GameState& state)
// {
//     // Implementation for creating a new cluster based on input events
//     // (This function can be expanded based on specific requirements)
    
//     GravitationalCluster newCluster;
//     // float radius = 15.0;
//     // newBody.setRadius(radius);

//     // double mass = MAX_MASS; // in kilograms
//     // newBody.setMass(mass);

//     // // Vector2D bodyVelocity = {100,0};
//     // // newBody.setNetVelocity(bodyVelocity);

//     // newBody.setCollisionEnabled(true);
    
//     // Vector2D bodyPosition = { event.button.x, event.button.y }; // position measured in pixel xy cooriifnates measured from top-left corner.
//     // newBody.setPosition(bodyPosition);

//     // Vector2D prevPosition = bodyPosition; // previous pixel position
//     // newBody.setPrevPosition(prevPosition);

//     // state.addBody(newBody);

//     // std::cout<<"New Body Created: "<<std::endl;
//     // std::cout<<newBody<<std::endl;

// }