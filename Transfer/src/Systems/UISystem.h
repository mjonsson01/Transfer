// File: Transfer/src/Systems/UISystem.h

#pragma once

// SDL3 Imports
#include "SDL3/SDL.h"
#include "SDL3_ttf/SDL_ttf.h"

// Custom Imports
#include "Core/GameState.h"
#include "Core/UIState.h"
#include "Entities/UIElement.h"

// Handles Logic of UI Components
class UISystem
{
    public:
        UISystem();
        ~UISystem();

        void ProcessUIFrame(GameState& state, UIState& UIState);

        void RenderUIElements(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont);
        
        void InitializeUIElements(UIState& UIState);
        
        void DeleteUIElements(UIState& UIState);
    private:
        // Add UI system members and methods here
};