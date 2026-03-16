// File: Transfer/src/Entities/UIElements/Sliders/MassSlider.cpp

#include "Entities/UIElements/Sliders/MassSlider.h"


MassSlider::MassSlider()
    : Slider(Orientation::Horizontal, SDL_FRect{100,300,200,6})
{
        minValue = 0.0;
        maxValue = MAX_MASS;
        knobRect = SDL_FRect{100, 300-7, 14, 20}; // default knob position and size
        setVisibility(true);
        setPosition(trackRect.x, trackRect.y);
        sliderValue = 0.0;
        // updateKnobPosition(MAX_MASS/2); WORKS :DDDDD
}




// MassSlider::MassSlider()
// {
//     setPosition(50.0f, 3*SCREEN_HEIGHT/4); //upper left corner of element. Other dimensions determined at the system level, or as global constants?
// }

// MassSlider::MassSlider()
// {
//     // Set default position and size for the Mass Slider
//     setPosition(50.0f, 3*SCREEN_HEIGHT/4); // top left corner of track rect position, will need to scale based on resolution
//     setSize(200.0f, 20.0f);    // Example size, will need to update based on resolution
//     knobRect = {50.0f, 3*SCREEN_HEIGHT/4-10.0f, 20.0f, 40.0f}; 
// }

// MassSlider::~MassSlider()
// {
//     // Cleanup if necessary
// }

// // --------- RENDER METHOD --------- //

// void MassSlider::renderElement(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont)
// {
//     // Get the positions for the trackRect and knobRect
//     // getTrackAndKnobPositions();

//     SDL_SetRenderDrawColor(renderer, ColorLibrary::Gray.r, ColorLibrary::Gray.g, ColorLibrary::Gray.b, ColorLibrary::Gray.a);
//     SDL_RenderFillRect(renderer, &trackRect);
//     SDL_SetRenderDrawColor(renderer, ColorLibrary::White.r, ColorLibrary::White.g, ColorLibrary::White.b, ColorLibrary::White.a);
//     SDL_RenderFillRect(renderer, &knobRect);
// }

// // --------- UTILITY METHOD FOR POSITIONS --------- //
// // void MassSlider::getTrackAndKnobPositions()
// // {
// //     // Get the track rect positions

// //     // // Get the slider knob positions
// // }