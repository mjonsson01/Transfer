// File: Transfer/src/DynamoEngine/Rendering/FontAtlas.hpp

#pragma once

// SDL Imports
#include <SDL3/SDL_surface.h>
#include <SDL3_ttf/SDL_ttf.h>

// Standard Library Imports
#include <array>
#include <cstddef>
#include <string_view>

namespace DynamoEngine
{

// Where one character lives in the font atlas texture and how to place it on the screen
struct GlyphMetrics
{
    float u1 = 0.0f, v1 = 0.0f, u2 = 0.0f, v2 = 0.0f; // Atlas UV quad
    float width = 0.0f, height = 0.0f;                // Specific glyph quad size in screen points
    float offset_x = 0.0f, offset_y = 0.0f;           // offset from the text cursor to the quad's top-left
    float advance_x = 0.0f;                           // how far to move the cursor after this glyph
};

// Printable ASCII (32-126) baked into one texture, plus the metrics needed to lay out text with it.
class FontAtlas
{
  public:
    // Bakes the glyphs from `font` into a new RGBA32 surface and fills in the metrics.
    // The caller owns the returned surface and must SDL_DestroySurface() it after uploading it to the GPU.
    SDL_Surface* buildAtlas(TTF_Font* font);

    // Metrics for one character; characters outside printable ASCII get empty (all-zero) metrics
    GlyphMetrics glyph(char character) const;

    // Total advance width of a string, in points
    float measureTextWidth(std::string_view text) const;

    float fontHeight() const { return m_font_height; }

  private:
    static constexpr char FIRST_GLYPH = 32; // ' '
    static constexpr char LAST_GLYPH = 126; // '~'
    static constexpr std::size_t GLYPH_COUNT = LAST_GLYPH - FIRST_GLYPH + 1;

    std::array<GlyphMetrics, GLYPH_COUNT> m_glyphs{}; // indexed by character - FIRST_GLYPH
    float m_font_height = 0.0f;
};

} // namespace DynamoEngine