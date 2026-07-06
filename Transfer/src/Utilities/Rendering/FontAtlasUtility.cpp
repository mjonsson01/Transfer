// File: Transfer/src/Utilities/Rendering/FontAtlasUtility.cpp

#include "Utilities/Rendering/FontAtlasUtility.hpp"

SDL_Surface* FontAtlasUtility::BuildAtlas(TTF_Font* font)
{
    if (!font)
        return nullptr;

    const int atlasWidth = 512;
    const int atlasHeight = 512;
    const int padding = 2;

    SDL_Surface* atlas = SDL_CreateSurface(atlasWidth, atlasHeight, SDL_PIXELFORMAT_RGBA32);
    if (!atlas)
        return nullptr;

    SDL_FillSurfaceRect(atlas, nullptr, SDL_MapSurfaceRGBA(atlas, 255, 255, 255, 0));

    fontHeight = static_cast<float>(TTF_GetFontHeight(font));

    int cursorX = 0;
    int cursorY = 0;
    int rowHeight = 0;

    for (char c = 32; c < 127; ++c)
    {
        SDL_Surface* glyphSurface =
            TTF_RenderGlyph_Blended(font, static_cast<Uint32>(c), SDL_Color{255, 255, 255, 255});
        if (!glyphSurface)
            continue;

        if (cursorX + glyphSurface->w + padding > atlasWidth)
        {
            cursorX = 0;
            cursorY += rowHeight + padding;
            rowHeight = 0;
        }

        SDL_Rect dstRect = {cursorX, cursorY, glyphSurface->w, glyphSurface->h};
        SDL_BlitSurface(glyphSurface, nullptr, atlas, &dstRect);

        int minx, maxx, miny, maxy, advance;
        TTF_GetGlyphMetrics(font, static_cast<Uint32>(c), &minx, &maxx, &miny, &maxy, &advance);

        GlyphMetrics metrics;
        metrics.u1 = static_cast<float>(cursorX) / atlasWidth;
        metrics.v1 = static_cast<float>(cursorY) / atlasHeight;
        metrics.u2 = static_cast<float>(cursorX + glyphSurface->w) / atlasWidth;
        metrics.v2 = static_cast<float>(cursorY + glyphSurface->h) / atlasHeight;
        metrics.width = static_cast<float>(glyphSurface->w);
        metrics.height = static_cast<float>(glyphSurface->h);
        metrics.offsetX = 0.0f;
        metrics.offsetY = 0.0f;
        metrics.advanceX = static_cast<float>(advance);

        glyphs[c] = metrics;

        cursorX += glyphSurface->w + padding;
        rowHeight = std::max(rowHeight, glyphSurface->h);

        SDL_DestroySurface(glyphSurface);
    }

    return atlas;
}
