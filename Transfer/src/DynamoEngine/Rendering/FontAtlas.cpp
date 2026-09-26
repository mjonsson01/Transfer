// File: Transfer/src/DynamoEngine/Rendering/FontAtlas.cpp

#include "DynamoEngine/Rendering/FontAtlas.hpp"

// SDL Imports
#include <SDL3/SDL_pixels.h>
#include <SDL3/SDL_rect.h>
#include <SDL3/SDL_stdinc.h>

// Standard Library Imports
#include <algorithm>

namespace DynamoEngine
{
SDL_Surface* FontAtlas::buildAtlas(TTF_Font* font)
{
    if (font == nullptr)
    {
        return nullptr;
    }

    constexpr int ATLAS_WIDTH = 512;
    constexpr int ATLAS_HEIGHT = 512;
    constexpr int PADDING = 2; // gap between glyphs so linear filtering doesn't bleed neighbors in

    SDL_Surface* atlas = SDL_CreateSurface(ATLAS_WIDTH, ATLAS_HEIGHT, SDL_PIXELFORMAT_RGBA32);
    if (atlas == nullptr)
    {
        return nullptr;
    }
    SDL_FillSurfaceRect(atlas, nullptr, SDL_MapSurfaceRGBA(atlas, 255, 255, 255, 0));

    m_font_height = static_cast<float>(TTF_GetFontHeight(font));

    int cursor_x = 0;
    int cursor_y = 0;
    int row_height = 0;

    for (char character = FIRST_GLYPH; character <= LAST_GLYPH; ++character)
    {
        SDL_Surface* glyph_surface =
            TTF_RenderGlyph_Blended(font, static_cast<Uint32>(character), SDL_Color{255, 255, 255, 255});
        if (glyph_surface == nullptr)
        {
            continue;
        }

        // Wrap to the next row when this glyph wouldn't fit
        if (cursor_x + glyph_surface->w + PADDING > ATLAS_WIDTH)
        {
            cursor_x = 0;
            cursor_y += row_height + PADDING;
            row_height = 0;
        }

        SDL_Rect destination = {cursor_x, cursor_y, glyph_surface->w, glyph_surface->h};
        SDL_BlitSurface(glyph_surface, nullptr, atlas, &destination);

        int advance = 0;
        TTF_GetGlyphMetrics(font, static_cast<Uint32>(character), nullptr, nullptr, nullptr, nullptr, &advance);

        GlyphMetrics& metrics = m_glyphs[static_cast<std::size_t>(character - FIRST_GLYPH)];
        metrics.u1 = static_cast<float>(cursor_x) / ATLAS_WIDTH;
        metrics.v1 = static_cast<float>(cursor_y) / ATLAS_HEIGHT;
        metrics.u2 = static_cast<float>(cursor_x + glyph_surface->w) / ATLAS_WIDTH;
        metrics.v2 = static_cast<float>(cursor_y + glyph_surface->h) / ATLAS_HEIGHT;
        metrics.width = static_cast<float>(glyph_surface->w);
        metrics.height = static_cast<float>(glyph_surface->h);
        metrics.advance_x = static_cast<float>(advance);

        cursor_x += glyph_surface->w + PADDING;
        row_height = std::max(row_height, glyph_surface->h);

        SDL_DestroySurface(glyph_surface);
    }

    return atlas;
}

GlyphMetrics FontAtlas::glyph(char character) const
{
    if (character < FIRST_GLYPH || character > LAST_GLYPH)
    {
        return GlyphMetrics{};
    }
    return m_glyphs[static_cast<std::size_t>(character - FIRST_GLYPH)];
}

float FontAtlas::measureTextWidth(std::string_view text) const
{
    float width = 0.0f;
    for (char character : text)
    {
        width += glyph(character).advance_x;
    }
    return width;
}
} // namespace DynamoEngine