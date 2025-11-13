// File: Transfer/src/Utilities/Colors.h

#pragma once

// SDL3 Imports
#include "SDL3/SDL_pixels.h"


struct ColorLibrary {
    inline static constexpr SDL_Color Red         {255, 0, 0, 255};
    inline static constexpr SDL_Color Green       {0, 255, 0, 255};
    inline static constexpr SDL_Color Blue        {0, 0, 255, 255};
    inline static constexpr SDL_Color White       {255, 255, 255, 255};
    inline static constexpr SDL_Color Black       {0, 0, 0, 255};
    inline static constexpr SDL_Color Transparent {0, 0, 0, 0};
    inline static constexpr SDL_Color Yellow      {255, 255, 0, 255};
    inline static constexpr SDL_Color Cyan        {0, 255, 255, 255};
    inline static constexpr SDL_Color Magenta     {255, 0, 255, 255};
    inline static constexpr SDL_Color Gray        {128, 128, 128, 255};
};