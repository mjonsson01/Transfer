// File: Transfer/Entities/src/UIElements/Sliders/Slider.cpp

#include "Entities/UIElements/Sliders/Slider.h"

Slider::Slider()
{
    orientation = Orientation::Horizontal;
    trackRect = SDL_FRect{0, 0, 0, 0};
    knobRect = {0, 0, 0, 0};
    hotZoneRect = {0, 0, 0, 0};
    sliderValue = 0.0;
    minValue = 0.0;
    maxValue = 0.0;
}

void Slider::slideMe(Vector2D positionOfEvent, double& returnedElementValue)
{

    // Track start positions
    float track_start_x = trackRect.x;
    float track_start_y = trackRect.y;

    // Usable track lengths (accounting for knob size)
    float track_length_x = trackRect.w - knobRect.w;
    float track_length_y = trackRect.h - knobRect.h;

    if (orientation == Orientation::Horizontal)
    {
        float new_x = positionOfEvent.xVal - (knobRect.w / 2.0f);

        // Clamp the new centered position
        if (new_x < track_start_x)
            new_x = track_start_x;
        if (new_x > track_start_x + track_length_x)
            new_x = track_start_x + track_length_x;

        // Map knob position to slider value (handles negative minValue)
        sliderValue = minValue + ((new_x - track_start_x) / track_length_x) * (maxValue - minValue);

        // Update knob position to reflect sliderValue
        knobRect.x = track_start_x + ((sliderValue - minValue) / (maxValue - minValue)) * track_length_x;
    }
    else // Vertical
    {
        float new_y = positionOfEvent.yVal - (knobRect.h / 2.0f);
        if (new_y < track_start_y)
            new_y = track_start_y;
        if (new_y > track_start_y + track_length_y)
            new_y = track_start_y + track_length_y;

        // Vertical sliders usually invert direction (top = max, bottom = min)
        sliderValue = maxValue - ((new_y - track_start_y) / track_length_y) * (maxValue - minValue);

        // Update knob position to match sliderValue
        knobRect.y = track_start_y + ((maxValue - sliderValue) / (maxValue - minValue)) * track_length_y;
    }

    // Return updated value
    returnedElementValue = sliderValue;
    return;
}

void Slider::renderMe(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont)
{
    SDL_SetRenderDrawColor(renderer, ColorLibrary::Gray.r, ColorLibrary::Gray.g, ColorLibrary::Gray.b,
                           ColorLibrary::Gray.a);
    SDL_RenderFillRect(renderer, &trackRect);
    SDL_SetRenderDrawColor(renderer, ColorLibrary::White.r, ColorLibrary::White.g, ColorLibrary::White.b,
                           ColorLibrary::White.a);
    SDL_RenderFillRect(renderer, &knobRect);
    if (UIState.getRenderDebug())
    {
        SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 60); // lighter alpha
        SDL_RenderFillRect(renderer, &hotZoneRect);
    }

    std::string slider_text = getDisplayText();
    SDL_Surface* text_surface =
        TTF_RenderText_Blended(UIFont, slider_text.c_str(), slider_text.length(), ColorLibrary::White);
    if (!text_surface)
    {
        SDL_Log("Text surface creation failed: %s", SDL_GetError());
        return;
    }
    SDL_Texture* text_texture = SDL_CreateTextureFromSurface(renderer, text_surface);
    if (!text_texture)
    {
        SDL_Log("Text texture creation failed: %s", SDL_GetError());
        return;
    }
    float width = static_cast<float>(text_surface->w);
    float height = static_cast<float>(text_surface->h);
    SDL_FRect dst_rect = {getX(), getY() + knobRect.h, width, height};
    SDL_RenderTexture(renderer, text_texture, nullptr, &dst_rect);
    SDL_DestroySurface(text_surface);
    SDL_DestroyTexture(text_texture);
}
