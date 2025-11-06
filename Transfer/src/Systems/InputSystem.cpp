#include "Systems/InputSystem.h"
#include <iostream>

InputSystem::InputSystem()
{
    // Initialize input system variables if needed
}

InputSystem::~InputSystem()
{
    // Cleanup if necessary
}
void InputSystem::ProcessSystemInputFrame(GameState& state)
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
            case SDL_EVENT_MOUSE_BUTTON_DOWN:
                createNewBody(event, state);

            case SDL_EVENT_KEY_DOWN:
                if (event.key.scancode == SDL_SCANCODE_SPACE){
                    GravitationalBody newBody;
                    float radius = 50.0;
                    newBody.setRadius(radius);

                    double mass = MAX_MASS; // in kilograms
                    newBody.setMass(mass);

                    Vector2D bodyPosition = { SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 }; // position measured in pixel xy coordinates measured from top-left corner.
                    newBody.setPosition(bodyPosition);

                    Vector2D prevPosition = bodyPosition; // previous pixel position
                    newBody.setPrevPosition(prevPosition);

                    state.addBody(newBody);
                }


        }
        
    }
        
        // Additional input handling logic can be added here
}

void InputSystem::createNewBody(SDL_Event& event, GameState& state)
{
    GravitationalBody newBody;

    float radius = 50.0;
    newBody.setRadius(radius);

    double mass = 1e12; // in kilograms
    newBody.setMass(mass);
    
    Vector2D bodyPosition = { event.button.x, event.button.y }; // position measured in pixel xy coordinates measured from top-left corner.
    newBody.setPosition(bodyPosition);

    Vector2D prevPosition = bodyPosition; // previous pixel position
    newBody.setPrevPosition(prevPosition);

    state.addBody(newBody);
};
