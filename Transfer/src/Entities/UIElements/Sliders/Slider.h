// File: // File: Transfer/src/Entities/UIElements/Sliders/Slider.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Entities/UIElements/UIElement.h"
#include "Entities/UIElements/UIElementTypeEnum.h"
#include "Core/UIState.h"
#include "Entities/MathStructures.h"

// Standard Library Imports
#include <string>

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
        virtual std::string getDisplayText() const { return std::to_string(sliderValue);}
        void updateMe(Vector2D positionOfEvent, double& returnedElementValue) override;
        double getSliderValue() { return sliderValue;}
    protected:
        Orientation orientation;
        SDL_FRect trackRect;
        SDL_FRect knobRect;
        double sliderValue;
        double minValue;
        double maxValue;
};