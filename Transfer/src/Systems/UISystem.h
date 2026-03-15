// File: Transfer/src/Systems/UISystem.h

#pragma once

// SDL3 Imports
#include "SDL3/SDL.h"
#include "SDL3_ttf/SDL_ttf.h"

// Custom Imports
#include "Core/GameState.h"
#include "Core/UIState.h"
#include "Entities/UIElements/UIElement.h"
#include "Entities/UIElements/Overlay/FPSCounter.h"
#include "Entities/UIElements/Sliders/MassSlider.h"


// Standard Library Imports
#include <vector>

// Owns Logic of UI Components and stores the UIElements while dispatching information to the GameState and UIState
class UISystem
{
    public:
        UISystem();
        ~UISystem();
        void CleanUp();
        std::vector<UIElement*>& getUIElements() { return allUIElements;}
    private:
        std::vector<UIElement*> allUIElements;

};
