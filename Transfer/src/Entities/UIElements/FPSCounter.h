// File: Transfer/src/Entities/FPSCounter.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Entities/UIElement.h"
#include "Core/UIState.h"
#include "Utilities/Colors.h"

// Standard Library Imports
#include <string>

class FPSCounter : public UIElement
{
    public:
        FPSCounter();
        ~FPSCounter() = default;
        
        virtual void renderElement(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont) override; // FPS Counter Specialty Render Method
    
    private:
};
