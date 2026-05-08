// File: Transfer/src/Entities/UIElements/Sliders/MassSlider.cpp

#include "Entities/UIElements/Sliders/MassSlider.h"

MassSlider::MassSlider() : Slider()
{
    orientation = Orientation::Horizontal;

    // Track size and position
    // trackRect = SDL_FRect{1200, 650, 300, 12};
    trackRect = SDL_FRect{17 * SCREEN_WIDTH / 24, 2 * SCREEN_HEIGHT / 3, 300, 12};
    knobRect = SDL_FRect{0, 0, 20, 30}; // will set x,y below

    float hotzonePaddingX = 30.0f;
    float hotzonePaddingY = 60.0f;
    hotZoneRect = SDL_FRect{trackRect.x - hotzonePaddingX / 2.0f, trackRect.y - hotzonePaddingY / 2.0f,
                            trackRect.w + hotzonePaddingX, trackRect.h + hotzonePaddingY};

    // Slider range
    maxValue = MAX_MASS;
    // minValue = 0.0;
    minValue = -MAX_MASS; // now supports negative values
    sliderValue = 0.0;    // start centered

    // Position knob based on slider value
    knobRect.x = trackRect.x + (sliderValue - minValue) / (maxValue - minValue) * trackRect.w - knobRect.w / 2.0f;
    knobRect.y = trackRect.y - (knobRect.h - trackRect.h) / 2.0f;

    setVisibility(true);
    setPosition(trackRect.x, trackRect.y);
    UIElementTypeIdentifier = UIElementType::MASS_SLIDER_INDEX;
}

void MassSlider::slideMe(Vector2D positionOfEvent, double& returnedElementValue)
{
    float track_start_x = trackRect.x;
    float track_length_x = trackRect.w - knobRect.w;

    // 1. Get normalized 0.0 to 1.0 position on the track
    float mouse_x = positionOfEvent.xVal - (knobRect.w / 2.0f);
    float t = std::clamp((mouse_x - track_start_x) / track_length_x, 0.0f, 1.0f);

    // 2. Map t (0..1) to a centered range (-1.0 to 1.0)
    // -1 is far left (min mass), 0 is center (zero mass), 1 is far right (max mass)
    double centered_t = (t * 2.0) - 1.0;
    double sign = (centered_t < 0) ? -1.0 : 1.0;

    // 3. Apply the Exponential Curve
    // We use a power function to give detail to small numbers.
    // 10^15 is massive, so we raise 10 to the power of (abs(centered_t) * 15)
    if (std::abs(centered_t) < 0.01)
    {
        sliderValue = 0.0; // "Snap" to zero in the middle
    }
    else
    {
        // This gives you a range of +/- 1 to +/- 1e15
        sliderValue = sign * std::pow(10.0, std::abs(centered_t) * 15.0);
    }

    // 4. Update the knob position (Standard linear for visual consistency)
    knobRect.x = track_start_x + (t * track_length_x);

    returnedElementValue = sliderValue;
}