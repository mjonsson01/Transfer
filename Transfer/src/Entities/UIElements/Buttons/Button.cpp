// File: Transfer/src/Entities/UIElements/Buttons/Button.cpp

#include "Entities/UIElements/Buttons/Button.hpp"

Button::Button()
{
    boundingRect = SDL_FRect{0.0f, 0.0f, 0.0f, 0.0};
    buttonSelected = false;
    setPosition(boundingRect.x, boundingRect.y);
    hotZoneRect = boundingRect;
}

void Button::clickMe(Vector2D positionOfEvent, UIState& UIState)
{
    std::string temp = altText;
    altText = displayText;
    displayText = temp;
    return;
}

void Button::buildGeometry(std::vector<UIElementVertex>& vertexBuffer, uint32_t zIndex,
                           const FontAtlasUtility& fontAtlas)
{
    float x1 = boundingRect.x;
    float y1 = boundingRect.y;
    float x2 = boundingRect.x + boundingRect.w;
    float y2 = boundingRect.y + boundingRect.h;

    float r = ColorLibrary::Gray.r / 255.0f;
    float g = ColorLibrary::Gray.g / 255.0f;
    float b = ColorLibrary::Gray.b / 255.0f;
    float a = ColorLibrary::Gray.a / 255.0f;
    uint32_t mode = 0;
    float u = 0.0f, v = 0.0f;

    vertexBuffer.push_back({x1, y1, u, v, r, g, b, a, zIndex, mode});
    vertexBuffer.push_back({x2, y1, u, v, r, g, b, a, zIndex, mode});
    vertexBuffer.push_back({x1, y2, u, v, r, g, b, a, zIndex, mode});
    vertexBuffer.push_back({x2, y1, u, v, r, g, b, a, zIndex, mode});
    vertexBuffer.push_back({x2, y2, u, v, r, g, b, a, zIndex, mode});
    vertexBuffer.push_back({x1, y2, u, v, r, g, b, a, zIndex, mode});

    std::string button_text = getDisplayText();

    // Calculate centering start positions
    float totalTextWidth = fontAtlas.CalculateTextWidth(button_text);
    float glyphFontHeight = fontAtlas.GetFontHeight();

    float cursorX = boundingRect.x + (boundingRect.w - totalTextWidth) / 2.0f;
    float cursorY = boundingRect.y + (boundingRect.h - glyphFontHeight) / 2.0f;

    uint32_t textMode = 1;
    for (char c : button_text)
    {
        // Get UV coordinates and pixel dimensions for this character from our pre-baked atlas
        GlyphMetrics metrics = fontAtlas.GetGlyph(c);

        float tx1 = cursorX + metrics.offsetX;
        float ty1 = cursorY + metrics.offsetY;
        float tx2 = tx1 + metrics.width;
        float ty2 = ty1 + metrics.height;

        // Push 6 vertices for this character, marked as Mode 1 (Textured)
        vertexBuffer.push_back({tx1, ty1, metrics.u1, metrics.v1, 1.0f, 1.0f, 1.0f, 1.0f, zIndex, textMode});
        vertexBuffer.push_back({tx2, ty1, metrics.u2, metrics.v1, 1.0f, 1.0f, 1.0f, 1.0f, zIndex, textMode});
        vertexBuffer.push_back({tx1, ty2, metrics.u1, metrics.v2, 1.0f, 1.0f, 1.0f, 1.0f, zIndex, textMode});
        vertexBuffer.push_back({tx2, ty1, metrics.u2, metrics.v1, 1.0f, 1.0f, 1.0f, 1.0f, zIndex, textMode});
        vertexBuffer.push_back({tx2, ty2, metrics.u2, metrics.v2, 1.0f, 1.0f, 1.0f, 1.0f, zIndex, textMode});
        vertexBuffer.push_back({tx1, ty2, metrics.u1, metrics.v2, 1.0f, 1.0f, 1.0f, 1.0f, zIndex, textMode});

        // Advance the cursor for the next character
        cursorX += metrics.advanceX;
    }
}