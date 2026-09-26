// File: Transfer/src/DynamoEngine/UI/UIGeometryBuilder.cpp

#include "DynamoEngine/UI/UIGeometryBuilder.hpp"

// Standard Library Imports
#include <cstdint>

namespace DynamoEngine
{
UIGeometryBuilder::UIGeometryBuilder(std::vector<UIVertex>& output_vertices, const FontAtlas& font_atlas)
    : m_vertices(output_vertices), m_font_atlas(font_atlas)
{
}

void UIGeometryBuilder::addRect(const SDL_FRect& rect, SDL_Color color)
{
    addQuad(rect, SDL_FRect{0.0f, 0.0f, 0.0f, 0.0f}, color, UIVertexMode::Solid);
}

void UIGeometryBuilder::addText(std::string_view text, Vector2F top_left, SDL_Color color)
{
    float cursor_x = top_left.x_val;
    for (char character : text)
    {
        GlyphMetrics metrics = m_font_atlas.glyph(character);

        SDL_FRect glyph_rect = {cursor_x + metrics.offset_x, top_left.y_val + metrics.offset_y, metrics.width,
                                metrics.height};
        SDL_FRect uv_rect = {metrics.u1, metrics.v1, metrics.u2 - metrics.u1, metrics.v2 - metrics.v1};
        addQuad(glyph_rect, uv_rect, color, UIVertexMode::Textured);

        cursor_x += metrics.advance_x;
    }
}

void UIGeometryBuilder::addTextCentered(std::string_view text, const SDL_FRect& area, SDL_Color color)
{
    Vector2F top_left = {area.x + (area.w - measureTextWidth(text)) / 2.0f, area.y + (area.h - fontHeight()) / 2.0f};
    addText(text, top_left, color);
}

void UIGeometryBuilder::addQuad(const SDL_FRect& rect, const SDL_FRect& uv_rect, SDL_Color color, UIVertexMode mode)
{
    const float x1 = rect.x;
    const float y1 = rect.y;
    const float x2 = rect.x + rect.w;
    const float y2 = rect.y + rect.h;

    const float u1 = uv_rect.x;
    const float v1 = uv_rect.y;
    const float u2 = uv_rect.x + uv_rect.w;
    const float v2 = uv_rect.y + uv_rect.h;

    const float r = color.r / 255.0f;
    const float g = color.g / 255.0f;
    const float b = color.b / 255.0f;
    const float a = color.a / 255.0f;
    const auto mode_value = static_cast<uint32_t>(mode);

    // Two triangles: (top-left, top-right, bottom-left) and (top-right, bottom-right, bottom-left)
    m_vertices.push_back({x1, y1, u1, v1, r, g, b, a, 0, mode_value});
    m_vertices.push_back({x2, y1, u2, v1, r, g, b, a, 0, mode_value});
    m_vertices.push_back({x1, y2, u1, v2, r, g, b, a, 0, mode_value});
    m_vertices.push_back({x2, y1, u2, v1, r, g, b, a, 0, mode_value});
    m_vertices.push_back({x2, y2, u2, v2, r, g, b, a, 0, mode_value});
    m_vertices.push_back({x1, y2, u1, v2, r, g, b, a, 0, mode_value});
}
} // namespace DynamoEngine