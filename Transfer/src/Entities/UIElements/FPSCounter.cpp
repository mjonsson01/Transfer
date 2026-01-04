// File: Transfer/src/Entities/FPSCounter.cpp

#include "FPSCounter.h"

FPSCounter::FPSCounter()
{
    // Set default position and size for the FPS counter
    setPosition(10, 10); // Top-left corner
    // setSize(100, 30);   // Width and height
}

// --------- RENDER METHOD --------- //

void FPSCounter::renderElement(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont)
{
    if (!UIState.getShowFPSCounter()) return;
    
    else
    {
        float fps = UIState.getFPS();
        std::string fps_text = "FPS: " + std::to_string(static_cast<int>(fps));
        SDL_Surface* text_surface = TTF_RenderText_Blended(UIFont, fps_text.c_str(), fps_text.length(), ColorLibrary::White);
        if (!text_surface) {
            SDL_Log("Text surface creation failed: %s", SDL_GetError());
            return;
        }

        SDL_Texture* text_texture = SDL_CreateTextureFromSurface(renderer, text_surface);
        if (!text_texture) {
            SDL_Log("Text texture creation failed: %s", SDL_GetError());
            return;
        }
        setSize(getScaleFactor()*static_cast<float>(text_surface->w), getScaleFactor()*static_cast<float>(text_surface->h));
        SDL_FRect dst_rect = {getX(), getY(), getWidth(), getHeight()};
        SDL_RenderTexture(renderer, text_texture, nullptr, &dst_rect);
        SDL_DestroySurface(text_surface);
        SDL_DestroyTexture(text_texture);
    }
}
