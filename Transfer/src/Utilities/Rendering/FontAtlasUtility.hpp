// File: Transfer/src/Utilities/Rendering/FontAtlasUtility.hpp

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Standard Library Imports
#include <algorithm>
#include <string>
#include <unordered_map>

struct GlyphMetrics
{
    float u1 = 0.0f, v1 = 0.0f, u2 = 0.0f, v2 = 0.0f; // atlas UV rect
    float width = 0.0f, height = 0.0f;                // glyph quad size in pixels
    float offsetX = 0.0f, offsetY = 0.0f;             // offset from cursor to quad origin
    float advanceX = 0.0f;                            // how far to move the cursor after this glyph
};

class FontAtlasUtility
{
  public:
    // Bakes glyphs 32-126 from `font` into a single RGBA32 surface and fills in
    // glyphs/fontHeight. Caller owns the returned surface and must
    // SDL_DestroySurface() it once it has been uploaded to the GPU.
    SDL_Surface* BuildAtlas(TTF_Font* font);

    GlyphMetrics GetGlyph(char c) const
    {
        auto it = glyphs.find(c);
        return it != glyphs.end() ? it->second : GlyphMetrics{};
    }

    float CalculateTextWidth(const std::string& text) const
    {
        float width = 0.0f;
        for (char c : text)
            width += GetGlyph(c).advanceX;
        return width;
    }

    float GetFontHeight() const { return fontHeight; }

  private:
    std::unordered_map<char, GlyphMetrics> glyphs;
    float fontHeight = 0.0f;
};
