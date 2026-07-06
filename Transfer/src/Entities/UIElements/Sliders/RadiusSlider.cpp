// File: Transfer/src/Entities/UIElements/Sliders/RadiusSlider.cpp

#include "Entities/UIElements/Sliders/RadiusSlider.hpp"

RadiusSlider::RadiusSlider() : Slider()
{
    orientation = Orientation::Horizontal;
    knobRect = SDL_FRect{0, 0, 20, 30}; // will set x,y below
    updateLayout(SCREEN_WIDTH, SCREEN_HEIGHT);
    // Slider range
    maxValue = MAX_RADIUS - 100; // For throttling
    minValue = 0.0;
    sliderValue = 0.0; // start centered
    setVisibility(true);
    UIElementID = UIElementIdentifier::RADIUS_SLIDER_INDEX;
}

void RadiusSlider::updateLayout(float windowWidth, float windowHeight)
{

    // Track size and position
    trackRect = SDL_FRect{9 * windowWidth / 24, 4 * windowHeight / 5, 300, 12};

    // Dead zone
    float deadzonePaddingX = 30.0f;
    float deadzonePaddingY = 60.0f;
    hotZoneRect = SDL_FRect{trackRect.x - deadzonePaddingX / 2.0f, trackRect.y - deadzonePaddingY / 2.0f,
                            trackRect.w + deadzonePaddingX, trackRect.h + deadzonePaddingY};

    // Position knob based on slider value
    float track_length_x = trackRect.w - knobRect.w;
    knobRect.x = trackRect.x + (float)((sliderValue - minValue) / (maxValue - minValue)) * track_length_x;
    // knobRect.x = trackRect.x + (sliderValue - minValue) / (maxValue - minValue) * trackRect.w - knobRect.w / 2.0f;
    knobRect.y = trackRect.y - (knobRect.h - trackRect.h) / 2.0f;

    setPosition(trackRect.x, trackRect.y);
}