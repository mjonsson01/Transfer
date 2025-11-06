#pragma once

#include "Core/GameState.h"
#include "Core/UIState.h"

// Handles Logic of UI Components
class UISystem
{
    public:
        UISystem();
        ~UISystem();

        void ProcessUIFrame(GameState& state, UIState& UIState);
    private:
        // Add UI system members and methods here
};