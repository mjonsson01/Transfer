// File: Transfer/Entities/src/UIElements/Sliders/Slider.cpp

#include "Entities/UIElements/Sliders/Slider.hpp"

Slider::Slider()
{
    orientation = Orientation::Horizontal;
    trackRect = SDL_FRect{0, 0, 0, 0};
    knobRect = {0, 0, 0, 0};
    hotZoneRect = {0, 0, 0, 0};
    sliderValue = 0.0;
    minValue = 0.0;
    max_value = 0.0;
}

void Slider::slideMe(DynamoEngine::Vector2D positionOfEvent, double& returnedElementValue, UIState& ui_state)
{

    // Track start positions
    float track_start_x = trackRect.x;
    float track_start_y = trackRect.y;

    // Usable track lengths (accounting for knob size)
    float track_length_x = trackRect.w - knobRect.w;
    float track_length_y = trackRect.h - knobRect.h;

    if (orientation == Orientation::Horizontal)
    {
        float new_x = positionOfEvent.x_val - (knobRect.w / 2.0f);

        // Clamp the new centered position
        if (new_x < track_start_x)
            new_x = track_start_x;
        if (new_x > track_start_x + track_length_x)
            new_x = track_start_x + track_length_x;

        // Map knob position to slider value (handles negative minValue)
        sliderValue = minValue + ((new_x - track_start_x) / track_length_x) * (max_value - minValue);

        // Update knob position to reflect sliderValue
        knobRect.x = track_start_x + ((sliderValue - minValue) / (max_value - minValue)) * track_length_x;
    }
    else // Vertical
    {
        float new_y = positionOfEvent.y_val - (knobRect.h / 2.0f);
        if (new_y < track_start_y)
            new_y = track_start_y;
        if (new_y > track_start_y + track_length_y)
            new_y = track_start_y + track_length_y;

        // Vertical sliders usually invert direction (top = max, bottom = min)
        sliderValue = max_value - ((new_y - track_start_y) / track_length_y) * (max_value - minValue);

        // Update knob position to match sliderValue
        knobRect.y = track_start_y + ((max_value - sliderValue) / (max_value - minValue)) * track_length_y;
    }

    // Return updated value
    returnedElementValue = sliderValue;
    playTickSoundIfMoved(ui_state);
    return;
}

void Slider::buildGeometry(DynamoEngine::UIGeometryBuilder& builder)
{
    builder.addRect(trackRect, ColorLibrary::Gray);
    builder.addRect(knobRect, ColorLibrary::White);
    builder.addText(getDisplayText(), {getX(), getY() + knobRect.h});
}

void Slider::playTickSoundIfMoved(UIState& ui_state)
{
    if (max_value == minValue)
        return; // avoid division by zero on an uninitialized/degenerate slider

    int current_tick =
        static_cast<int>(std::round((sliderValue - minValue) / (max_value - minValue) * NUM_SLIDER_TICKS));

    if (current_tick != lastTickIndex)
    {
        ui_state.QueueSoundEffect("SliderTick");
        lastTickIndex = current_tick;
    }
}