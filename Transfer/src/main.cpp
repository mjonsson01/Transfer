// File: Transfer/src/main.cpp
// SDL3 imports
#include <SDL3/SDL.h>
#include <SDL3/SDL_events.h>

// Custom Imports
#include "Core/Game.h"

int main(int argcc, char* argv[])
{
    Game gameInstance;
    gameInstance.StartGame();
}