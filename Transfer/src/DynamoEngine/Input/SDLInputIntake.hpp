// File: Transfer/src/DynamoEngine/Input/SDLInputIntake.hpp

#pragma once

// SDL Imports
#include <SDL3/SDL_events.h> // IWYU pragma: export

// Custom Imports
#include "DynamoEngine/Input/InputEvent.hpp" // IWYU pragma: export

// Standard Library Imports
#include <vector> // IWYU pragma: export

namespace DynamoEngine
{

// Translates one SDL event into the engine's InputEvent vocabulary.
// Returns false for SDL events that input doesn't care about (out_event is then left untouched).
bool translateSDLEvent(const SDL_Event& sdl_event, InputEvent& output_event);

class SDLInputIntake
{
  public:
    SDLInputIntake() = default;
    ~SDLInputIntake() = default;

  public:
    // Drains SDL's queue, appending one InputEvent per relevant SDL event. Does NOT clear out_events.
    void pollEvents(std::vector<InputEvent>& output_events);
};

} // namespace DynamoEngine