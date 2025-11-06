#pragma once
#include <SDL3/SDL.h>

struct Star {
    SDL_Texture* texture;   // pre-created texture
    SDL_FRect dstRect;      // position and size
    float baseAlpha;        // 0.0 - 1.0
    float twinkleSpeed;     // multiplier for sine oscillation
    float currentAlpha;     // updated each frame
};