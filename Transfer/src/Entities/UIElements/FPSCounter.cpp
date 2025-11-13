// File: Transfer/src/Entities/FPSCounter.cpp

#include "FPSCounter.h"


FPSCounter::FPSCounter()
{
    // Set default position and size for the FPS counter
    setPosition(10, 10); // Top-left corner
    // setSize(100, 30);   // Width and height
}
void FPSCounter::renderElement(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont)
{
    if (!UIState.getShowFPSCounter())
    {
        return;
    }
    else
    {
        float fps = UIState.getFPS();
        std::string fpsText = "FPS: " + std::to_string(static_cast<int>(fps));
        SDL_Surface* textSurface = TTF_RenderText_Blended(UIFont, fpsText.c_str(), fpsText.length(), ColorLibrary::White);
        if (!textSurface) {
            SDL_Log("Text surface creation failed: %s", SDL_GetError());
            return;
        }

        SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
        if (!textTexture) {
            SDL_Log("Text texture creation failed: %s", SDL_GetError());
            return;
        }
        setSize(getScaleFactor()*static_cast<float>(textSurface->w), getScaleFactor()*static_cast<float>(textSurface->h));
        SDL_FRect dstRect = {getX(), getY(), getWidth(), getHeight()};
        SDL_RenderTexture(renderer, textTexture, nullptr, &dstRect);
        SDL_DestroySurface(textSurface);
        SDL_DestroyTexture(textTexture);
    }
}
