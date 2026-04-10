// File: Transfer/Entities/src/UIElements/Sliders/Slider.cpp

#include "Entities/UIElements/Sliders/Slider.h"

Slider::Slider() {
    orientation = Orientation::Horizontal;
    trackRect = SDL_FRect{0, 0, 0, 0};
    knobRect = {0, 0, 0, 0};
    hotZoneRect = {0, 0, 0, 0};
    sliderValue = 0.0;
    minValue = 0.0;
    maxValue = 0.0;
}

void Slider::updateMe(Vector2D positionOfEvent, double& returnedElementValue) {
    // Track start positions
    float trackStartX = trackRect.x;
    float trackStartY = trackRect.y;

    // Usable track lengths (accounting for knob size)
    float trackLengthX = trackRect.w - knobRect.w;
    float trackLengthY = trackRect.h - knobRect.h;

    // Clamp the event position to the track bounds
    float clampedX = positionOfEvent.xVal;
    float clampedY = positionOfEvent.yVal;

    if (orientation == Orientation::Horizontal) {
        if (clampedX < trackStartX)
            clampedX = trackStartX;
        if (clampedX > trackStartX + trackLengthX)
            clampedX = trackStartX + trackLengthX;

        // Map knob position to slider value (handles negative minValue)
        sliderValue = minValue + ((clampedX - trackStartX) / trackLengthX) *
                                     (maxValue - minValue);

        // Update knob position to reflect sliderValue
        knobRect.x =
            trackStartX +
            ((sliderValue - minValue) / (maxValue - minValue)) * trackLengthX;
    } else // Vertical
    {
        if (clampedY < trackStartY)
            clampedY = trackStartY;
        if (clampedY > trackStartY + trackLengthY)
            clampedY = trackStartY + trackLengthY;

        // Vertical sliders usually invert direction (top = max, bottom = min)
        sliderValue = maxValue - ((clampedY - trackStartY) / trackLengthY) *
                                     (maxValue - minValue);

        // Update knob position to match sliderValue
        knobRect.y =
            trackStartY +
            ((maxValue - sliderValue) / (maxValue - minValue)) * trackLengthY;
    }

    // Return updated value
    returnedElementValue = sliderValue;
}

void Slider::renderMe(SDL_Renderer* renderer, UIState& UIState,
                      TTF_Font* UIFont) {
    SDL_SetRenderDrawColor(renderer, ColorLibrary::Gray.r, ColorLibrary::Gray.g,
                           ColorLibrary::Gray.b, ColorLibrary::Gray.a);
    SDL_RenderFillRect(renderer, &trackRect);
    SDL_SetRenderDrawColor(renderer, ColorLibrary::White.r,
                           ColorLibrary::White.g, ColorLibrary::White.b,
                           ColorLibrary::White.a);
    SDL_RenderFillRect(renderer, &knobRect);
    if (UIState.getRenderDebug()) {
        SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 60); // lighter alpha
        SDL_RenderFillRect(renderer, &hotZoneRect);
    }

    std::string slider_text = getDisplayText();
    SDL_Surface* text_surface = TTF_RenderText_Blended(
        UIFont, slider_text.c_str(), slider_text.length(), ColorLibrary::White);
    if (!text_surface) {
        SDL_Log("Text surface creation failed: %s", SDL_GetError());
        return;
    }
    SDL_Texture* text_texture =
        SDL_CreateTextureFromSurface(renderer, text_surface);
    if (!text_texture) {
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
