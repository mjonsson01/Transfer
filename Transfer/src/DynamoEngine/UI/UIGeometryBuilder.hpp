// File: Transfer/src/DynamoEngine/UI/UIGeometryBuilder.hpp

#pragma once

// SDL Imports
#include <SDL3/SDL_pixels.h>
#include <SDL3/SDL_rect.h>

// Custom Imports
#include "DynamoEngine/Math/Vector2.hpp"
#include "DynamoEngine/Rendering/FontAtlas.hpp"
#include "DynamoEngine/Rendering/UIVertex.hpp"

// Standard Library Imports
#include <string_view>
#include <vector>

namespace DynamoEngine
{
class UIGeometryBuilder
{
  public:
    UIGeometryBuilder(std::vector<UIVertex>& output_vertices,
                      const FontAtlas& font_atlas); // appends to the output_vertices
    // High level function to add an axis-aligned rectangle
    void addRect(const SDL_FRect& rect, SDL_Color color);

    // Single-line text with top left corner at top_left
    void addText(std::string_view text, Vector2F top_left, SDL_Color color = SDL_Color{255, 255, 255, 255});

    // Single-line text centered inside area
    void addTextCentered(std::string_view text, const SDL_FRect& area, SDL_Color color = SDL_Color{255, 255, 255, 255});

    // Layout helpers to size elements around specific text
    float measureTextWidth(std::string_view text) const { return m_font_atlas.measureTextWidth(text); }
    float fontHeight() const { return m_font_atlas.fontHeight(); }

  private:
    // Pushes 2 triangles to cover rect textured with the uv_rect from the atlas
    void addQuad(const SDL_FRect& rect, const SDL_FRect& uv_rect, SDL_Color color, UIVertexMode);

    std::vector<UIVertex>& m_vertices;
    const FontAtlas& m_font_atlas;
};
} // namespace DynamoEngine