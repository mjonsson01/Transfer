// File: Transfer/src/Entities/UIElements/Overlay/FPSCounter.cpp

#include "Entities/UIElements/Overlay/FPSCounter.h"

FPSCounter::FPSCounter()
{
    setPosition(10.0f, 10.0f);
    setVisibility(true);
    UIElementIdentifierIdentifier = UIElementIdentifier::FPS_COUNTER_INDEX;
}

void FPSCounter::renderMe(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont)
{
    float fps = UIState.getFPS();
    std::string fps_text = "FPS: " + std::to_string(static_cast<int>(fps));
    SDL_Surface* text_surface =
        TTF_RenderText_Blended(UIFont, fps_text.c_str(), fps_text.length(), ColorLibrary::White);
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
    SDL_FRect dst_rect = {getX(), getY(), width, height};
    hotZoneRect = dst_rect;

    if (UIState.getRenderDebug())
    {
        SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 60); // lighter alpha
        SDL_RenderFillRect(renderer, &hotZoneRect);
    }
    SDL_RenderTexture(renderer, text_texture, nullptr, &dst_rect);
    SDL_DestroySurface(text_surface);
    SDL_DestroyTexture(text_texture);
}
