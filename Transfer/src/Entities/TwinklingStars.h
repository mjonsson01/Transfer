// File: Transfer/src/Entities/TwinklingStars.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>

struct TwinklingStar {
    SDL_Texture* texture;   // pre-created texture
    SDL_FRect dstRect;      // position and size
    float baseAlpha;        // 0.0 - 1.0 handles opacity baseline
    float twinkleSpeed;     // multiplier for sine oscillation
    float currentAlpha;     // updated each frame
};