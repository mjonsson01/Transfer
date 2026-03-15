// File: Transfer/Entities/src/UIElements/Sliders/Slider.cpp

#include "Entities/UIElements/Sliders/Slider.h"


Slider::Slider()
{

}


void Slider::renderMe(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UI_Font)
{

}

void Slider::updateKnobPosition()
{

}

double Slider::convertKnobPositionToValue(SDL_FRect knobRect)
{
    return 0.0;
}

SDL_FRect Slider::convertValueToKnobPosition(UIState& UIState) const
{
    return SDL_FRect{0,0,0,0};
}