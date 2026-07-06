// File: Transfer/src/Entities/UIElements/Overlay/FPSCounter.cpp

#include "Entities/UIElements/Overlay/FPSCounter.hpp"

FPSCounter::FPSCounter()
{
    updateLayout(SCREEN_WIDTH, SCREEN_HEIGHT);
    setVisibility(true);
    UIElementID = UIElementIdentifier::FPS_COUNTER_INDEX;
}

void FPSCounter::updateMe(UIState& UIState)
{
    float fps_local = UIState.getFPS();
    fps = static_cast<int>(fps_local);
}

// Local pushtext helper
static void pushText(std::vector<UIElementVertex>& vertexBuffer, const std::string& text, float startX, float startY,
                     const FontAtlasUtility& fontAtlas, uint32_t zIndex)
{
    float cursorX = startX;
    uint32_t textMode = 1;
    for (char c : text)
    {
        GlyphMetrics metrics = fontAtlas.GetGlyph(c);

        float tx1 = cursorX + metrics.offsetX;
        float ty1 = startY + metrics.offsetY;
        float tx2 = tx1 + metrics.width;
        float ty2 = ty1 + metrics.height;

        vertexBuffer.push_back({tx1, ty1, metrics.u1, metrics.v1, 1.0f, 1.0f, 1.0f, 1.0f, zIndex, textMode});
        vertexBuffer.push_back({tx2, ty1, metrics.u2, metrics.v1, 1.0f, 1.0f, 1.0f, 1.0f, zIndex, textMode});
        vertexBuffer.push_back({tx1, ty2, metrics.u1, metrics.v2, 1.0f, 1.0f, 1.0f, 1.0f, zIndex, textMode});
        vertexBuffer.push_back({tx2, ty1, metrics.u2, metrics.v1, 1.0f, 1.0f, 1.0f, 1.0f, zIndex, textMode});
        vertexBuffer.push_back({tx2, ty2, metrics.u2, metrics.v2, 1.0f, 1.0f, 1.0f, 1.0f, zIndex, textMode});
        vertexBuffer.push_back({tx1, ty2, metrics.u1, metrics.v2, 1.0f, 1.0f, 1.0f, 1.0f, zIndex, textMode});

        cursorX += metrics.advanceX;
    }
}

void FPSCounter::buildGeometry(std::vector<UIElementVertex>& vertexBuffer, uint32_t zIndex,
                               const FontAtlasUtility& fontAtlas)
{
    std::string fps_root = getDisplayText();
    std::string fps_text = "FPS: " + fps_root;
    pushText(vertexBuffer, fps_text, getX(), getY(), fontAtlas, zIndex);
}

void FPSCounter::updateLayout(float windowWidth, float windowHeight)
{
    setPosition(windowWidth / SCREEN_WIDTH * 10, windowHeight / SCREEN_HEIGHT * 10);
}