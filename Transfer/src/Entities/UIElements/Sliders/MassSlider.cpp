// File: Transfer/src/Entities/UIElements/Sliders/MassSlider.cpp

#include "Entities/UIElements/Sliders/MassSlider.hpp"

MassSlider::MassSlider() : Slider()
{
    orientation = Orientation::Horizontal;
    knobRect = SDL_FRect{0, 0, 20, 30}; // will set x,y below
    updateLayout(SCREEN_WIDTH, SCREEN_HEIGHT);

    // Slider range
    maxValue = MAX_MASS / 10;
    // minValue = 0.0;
    minValue = -MAX_MASS / 10; // now supports negative values
    sliderValue = 0.0;         // start centered
    setVisibility(true);
    UIElementID = UIElementIdentifier::MASS_SLIDER_INDEX;
}

void MassSlider::slideMe(Vector2D positionOfEvent, double& returnedElementValue, UIState& UIState)
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
    playTickSoundIfMoved(UIState);
    return;
}

void MassSlider::updateLayout(float windowWidth, float windowHeight)
{
    trackRect = SDL_FRect{17 * windowWidth / 24, 4 * windowHeight / 5, 300, 12};
    knobRect = SDL_FRect{0, 0, 20, 30};

    float hotzonePaddingX = 30.0f;
    float hotzonePaddingY = 60.0f;
    hotZoneRect = SDL_FRect{trackRect.x - hotzonePaddingX / 2.0f, trackRect.y - hotzonePaddingY / 2.0f,
                            trackRect.w + hotzonePaddingX, trackRect.h + hotzonePaddingY};

    double centered_t;
    if (sliderValue == 0.0)
    {
        centered_t = 0.0;
    }
    else
    {
        double sign = (sliderValue < 0) ? -1.0 : 1.0;
        centered_t = sign * (std::log10(std::abs(sliderValue)) / 15.0);
    }
    double t = (centered_t + 1.0) / 2.0;

    float track_length_x = trackRect.w - knobRect.w;
    knobRect.x = trackRect.x + (float)t * track_length_x;
    knobRect.y = trackRect.y - (knobRect.h - trackRect.h) / 2.0f;

    setPosition(trackRect.x, trackRect.y);
}

void MassSlider::playTickSoundIfMoved(UIState& UIState)
{
    double centered_t;
    if (sliderValue == 0.0)
    {
        centered_t = 0.0;
    }
    else
    {
        double sign = (sliderValue < 0) ? -1.0 : 1.0;
        centered_t = sign * (std::log10(std::abs(sliderValue)) / 15.0);
    }

    int currentTick = static_cast<int>(std::round(centered_t * NUM_SLIDER_TICKS));

    if (currentTick != lastTickIndex)
    {
        UIState.QueueSoundEffect("SliderTick");
        lastTickIndex = currentTick;
    }
}