// File: Transfer/src/Entities/UIElements/Buttons/Button.cpp


#include "Entities/UIElements/Buttons/Button.h"


Button::Button()
{
    boundingRect = SDL_FRect{0.0f, 0.0f, 0.0f, 0.0};
    buttonSelected = false;
    setPosition(boundingRect.x, boundingRect.y);
    hotZoneRect = boundingRect;
}


void Button::clickMe(Vector2D positionOfEvent)
{
    std::string temp = altText;
    altText = displayText;
    displayText = temp;
    return;
}

void Button::updateMe(Vector2D positionOfEvent, double& returnedElementValue)
{
    clickMe(positionOfEvent);
}
void Button::renderMe(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont)
{
    SDL_SetRenderDrawColor(renderer, ColorLibrary::Gray.r, ColorLibrary::Gray.g, ColorLibrary::Gray.b,
                           ColorLibrary::Gray.a);
    SDL_RenderFillRect(renderer, &boundingRect);
    if (UIState.getRenderDebug())
    {
        SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 60); // lighter alpha
        SDL_RenderFillRect(renderer, &hotZoneRect);
    }
    std::string button_text = getDisplayText();
    SDL_Surface* text_surface =
        TTF_RenderText_Blended(UIFont, button_text.c_str(), button_text.length(), ColorLibrary::White);
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
    float text_width = static_cast<float>(text_surface->w);
    float text_height = static_cast<float>(text_surface->h);
    SDL_FRect dst_rect;
    dst_rect.w = text_width;
    dst_rect.h = text_height;
    dst_rect.x = boundingRect.x + (boundingRect.w - text_width) / 2.0f;
    dst_rect.y = boundingRect.y + (boundingRect.h - text_height) / 2.0f;
    SDL_RenderTexture(renderer, text_texture, nullptr, &dst_rect);
    SDL_DestroySurface(text_surface);
    SDL_DestroyTexture(text_texture);
    return;
}