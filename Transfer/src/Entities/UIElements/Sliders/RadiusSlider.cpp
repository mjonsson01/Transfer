// File: Transfer/src/Entities/UIElements/Sliders/RadiusSlider.cpp

#include "Entities/UIElements/Sliders/RadiusSlider.hpp"

RadiusSlider::RadiusSlider() : Slider()
{
    orientation = Orientation::Horizontal;

    // Track size and position
    // trackRect = SDL_FRect{1200, 650, 300, 12};
    trackRect = SDL_FRect{9 * SCREEN_WIDTH / 24, 2 * SCREEN_HEIGHT / 3, 300, 12};
    knobRect = SDL_FRect{0, 0, 20, 30}; // will set x,y below
    // Dead zone
    float deadzonePaddingX = 30.0f;
    float deadzonePaddingY = 60.0f;
    hotZoneRect = SDL_FRect{trackRect.x - deadzonePaddingX / 2.0f, trackRect.y - deadzonePaddingY / 2.0f,
                            trackRect.w + deadzonePaddingX, trackRect.h + deadzonePaddingY};

    // Slider range
    maxValue = MAX_RADIUS - 50;
    minValue = 0.0;
    sliderValue = 0.0; // start centered

    // Position knob based on slider value
    knobRect.x = trackRect.x + (sliderValue - minValue) / (maxValue - minValue) * trackRect.w - knobRect.w / 2.0f;
    knobRect.y = trackRect.y - (knobRect.h - trackRect.h) / 2.0f;

    setVisibility(true);
    setPosition(trackRect.x, trackRect.y);
    UIElementID = UIElementIdentifier::RADIUS_SLIDER_INDEX;
}
