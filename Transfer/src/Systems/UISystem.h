#pragma once

#include "SDL3/SDL.h"

#include "Core/GameState.h"
#include "Core/UIState.h"
// Handles Logic of UI Components
class UISystem
{
    public:
        UISystem();
        ~UISystem();

        void ProcessUIFrame(GameState& state, UIState& UIState);

        void RenderUIElements(SDL_Renderer* renderer, UIState& UIState);
    private:
        // Add UI system members and methods here
};