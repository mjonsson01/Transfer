// File: Transfer/Entities/src/UIElements/Sliders/Slider.cpp

#include "Entities/UIElements/Sliders/Slider.h"


Slider::Slider(Orientation orientation, SDL_FRect track)
    : orientation(orientation), trackRect(track)
{
    knobRect = {0.0, 0.0, 0.0, 0.0}; // default;
    sliderValue = 0.0;
}


void Slider::renderMe(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont)
{    
    SDL_SetRenderDrawColor(renderer, ColorLibrary::Gray.r, ColorLibrary::Gray.g, ColorLibrary::Gray.b, ColorLibrary::Gray.a);
    SDL_RenderFillRect(renderer, &trackRect);
    SDL_SetRenderDrawColor(renderer, ColorLibrary::White.r, ColorLibrary::White.g, ColorLibrary::White.b, ColorLibrary::White.a);
    SDL_RenderFillRect(renderer, &knobRect);
    
    std::string slider_text = getDisplayText();
    SDL_Surface* text_surface = TTF_RenderText_Blended(UIFont, slider_text.c_str(), slider_text.length(), ColorLibrary::White);
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
    SDL_FRect dst_rect = {getX(), getY()+knobRect.h, width, height};
    SDL_RenderTexture(renderer, text_texture, nullptr, &dst_rect);
    SDL_DestroySurface(text_surface);
    SDL_DestroyTexture(text_texture);    
}

void Slider::updateKnobPosition(double desiredValue)
{
    if (desiredValue < minValue)
    {
        sliderValue = minValue;
    }
    else if (desiredValue > maxValue)
    {
        sliderValue = maxValue;
    }
    else
        sliderValue = desiredValue;

    Vector2D new_knob_position = convertValueToKnobPosition(sliderValue);
    
    if (orientation == Orientation::Horizontal)
    {
        knobRect.x = new_knob_position.xVal;
    }
    else if (orientation == Orientation::Vertical)
    {
        knobRect.y = new_knob_position.yVal;
    }
}

// double Slider::convertKnobPositionToValue(SDL_FRect knobPosition)
// {
//     return 0.0;
// }

Vector2D Slider::convertValueToKnobPosition(double valueToConvert) const
{
    double x = getX() + (valueToConvert/maxValue)*trackRect.w;
    double y = getY() + (valueToConvert/maxValue)*trackRect.h;
    return Vector2D(x,y);
}