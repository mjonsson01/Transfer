// File: Transfer/src/Entities/UIElements/MassSlider.cpp

#include "MassSlider.h"

MassSlider::MassSlider()
{
    // Set default position and size for the Mass Slider
    setPosition(50.0f, 1040.0f); // top left corner of track rect position //probably need to scale based on resolution?
    setSize(200.0f, 20.0f);    // Example size // Will also need to scale on resolution?
}

MassSlider::~MassSlider()
{
    // Cleanup if necessary
}

void MassSlider::renderElement(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont)
{
    // Render the slider track
    trackRect.x = getX();
    trackRect.y = getY();
    trackRect.w = getWidth();
    trackRect.h = getHeight();

    // Render the slider knob
    knobRect.x = getX();
    knobRect.y = getY()-0.5f*getHeight();
    knobRect.w = 0.1f * getWidth();
    knobRect.h = 2.0f*getHeight();

    SDL_SetRenderDrawColor(renderer, ColorLibrary::Gray.r, ColorLibrary::Gray.g, ColorLibrary::Gray.b, ColorLibrary::Gray.a);
    SDL_RenderFillRect(renderer, &trackRect);
    // SDL_Log("Rendered Mass Slider Track at (%.2f, %.2f) with size (%.2f x %.2f)", trackRect.x, trackRect.y, trackRect.w, trackRect.h);  
    SDL_SetRenderDrawColor(renderer, ColorLibrary::White.r, ColorLibrary::White.g, ColorLibrary::White.b, ColorLibrary::White.a);
    SDL_RenderFillRect(renderer, &knobRect);
    // SDL_Log("Rendered Mass Slider Knob at (%.2f, %.2f) with size (%.2f x %.2f)", knobRect.x, knobRect.y, knobRect.w, knobRect.h);
    // Calculate knob position based on selected mass value
    // float massRange = maxMassValue - minMassValue;
    // float normalizedMass = (selected

}