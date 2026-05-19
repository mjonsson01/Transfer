// File: Transfer/src/Entities/Sound/SoundEffect.h

#pragma once

#include <SDL3/SDL.h>

struct SoundEffect
{
    SDL_AudioSpec spec{};
    Uint8* buffer = nullptr;
    Uint32 length = 0;
};