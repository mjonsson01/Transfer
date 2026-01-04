// File: Transfer/src/Entities/UIElements/MassSlider.cpp

#include "MassSlider.h"

MassSlider::MassSlider()
{
    // Set default position and size for the Mass Slider
    setPosition(50.0f, 3*SCREEN_HEIGHT/4); // top left corner of track rect position, will need to scale based on resolution
    setSize(200.0f, 20.0f);    // Example size, will need to update based on resolution
}

MassSlider::~MassSlider()
{
    // Cleanup if necessary
}

// --------- RENDER METHOD --------- //

void MassSlider::renderElement(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont)
{
    // Get the positions for the trackRect and knobRect
    getTrackAndKnobPositions();

    SDL_SetRenderDrawColor(renderer, ColorLibrary::Gray.r, ColorLibrary::Gray.g, ColorLibrary::Gray.b, ColorLibrary::Gray.a);
    SDL_RenderFillRect(renderer, &trackRect);
    SDL_SetRenderDrawColor(renderer, ColorLibrary::White.r, ColorLibrary::White.g, ColorLibrary::White.b, ColorLibrary::White.a);
    SDL_RenderFillRect(renderer, &knobRect);
}

// --------- UTILITY METHOD FOR POSITIONS --------- //
void MassSlider::getTrackAndKnobPositions()
{
    // Get the track rect positions
    trackRect.x = getX();
    trackRect.y = getY();
    trackRect.w = getWidth();
    trackRect.h = getHeight();

    // Get the slider knob positions
    knobRect.x = getX();
    knobRect.y = getY()-0.5f*getHeight();
    knobRect.w = 0.1f * getWidth();
    knobRect.h = 2.0f*getHeight();
}