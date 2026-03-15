// File: // File: Transfer/src/Entities/UIElements/Sliders/Slider.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Entities/UIElements/UIElement.h"
#include "Core/UIState.h"

// Standard Library Imports

enum Orientation 
{
    Horizontal = 0,
    Vertical = 1
};

class Slider : public UIElement
{
    public:
        Slider();
        ~Slider() = default;
        void renderMe(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont) override;
        void updateKnobPosition();
        double convertKnobPositionToValue(SDL_FRect knobRect);
        SDL_FRect convertValueToKnobPosition(UIState& UIState) const;
        double getSliderValue() { return sliderValue;}

    public:
    private:
        Orientation orientation;
        const SDL_FRect trackRect;
        SDL_FRect knobRect;
        double sliderValue;
};