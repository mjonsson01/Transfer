// File: Transfer/Entities/src/UIElements/Sliders/Slider.cpp

#include "Entities/UIElements/Sliders/Slider.hpp"

Slider::Slider()
{
    orientation = Orientation::Horizontal;
    trackRect = SDL_FRect{0, 0, 0, 0};
    knobRect = {0, 0, 0, 0};
    hotZoneRect = {0, 0, 0, 0};
    sliderValue = 0.0;
    minValue = 0.0;
    max_value = 0.0;
}

void Slider::slideMe(DynamoEngine::Vector2D positionOfEvent, double& returnedElementValue, UIState& uiState)
{

    // Track start positions
    float track_start_x = trackRect.x;
    float track_start_y = trackRect.y;

    // Usable track lengths (accounting for knob size)
    float track_length_x = trackRect.w - knobRect.w;
    float track_length_y = trackRect.h - knobRect.h;

    if (orientation == Orientation::Horizontal)
    {
        float new_x = positionOfEvent.x_val - (knobRect.w / 2.0f);

        // Clamp the new centered position
        if (new_x < track_start_x)
            new_x = track_start_x;
        if (new_x > track_start_x + track_length_x)
            new_x = track_start_x + track_length_x;

        // Map knob position to slider value (handles negative minValue)
        sliderValue = minValue + ((new_x - track_start_x) / track_length_x) * (max_value - minValue);

        // Update knob position to reflect sliderValue
        knobRect.x = track_start_x + ((sliderValue - minValue) / (max_value - minValue)) * track_length_x;
    }
    else // Vertical
    {
        float new_y = positionOfEvent.y_val - (knobRect.h / 2.0f);
        if (new_y < track_start_y)
            new_y = track_start_y;
        if (new_y > track_start_y + track_length_y)
            new_y = track_start_y + track_length_y;

        // Vertical sliders usually invert direction (top = max, bottom = min)
        sliderValue = max_value - ((new_y - track_start_y) / track_length_y) * (max_value - minValue);

        // Update knob position to match sliderValue
        knobRect.y = track_start_y + ((max_value - sliderValue) / (max_value - minValue)) * track_length_y;
    }

    // Return updated value
    returnedElementValue = sliderValue;
    playTickSoundIfMoved(uiState);
    return;
}

static void pushQuad(std::vector<DynamoEngine::UIVertex>& vertexBuffer, const SDL_FRect& rect, SDL_Color color,
                     uint32_t z_index)
{
    float x1 = rect.x;
    float y1 = rect.y;
    float x2 = rect.x + rect.w;
    float y2 = rect.y + rect.h;

    float r = color.r / 255.0f;
    float g = color.g / 255.0f;
    float b = color.b / 255.0f;
    float a = color.a / 255.0f;
    uint32_t mode = static_cast<uint32_t>(DynamoEngine::UIVertexMode::Solid);
    float u = 0.0f, v = 0.0f;

    vertexBuffer.push_back({x1, y1, u, v, r, g, b, a, z_index, mode});
    vertexBuffer.push_back({x2, y1, u, v, r, g, b, a, z_index, mode});
    vertexBuffer.push_back({x1, y2, u, v, r, g, b, a, z_index, mode});
    vertexBuffer.push_back({x2, y1, u, v, r, g, b, a, z_index, mode});
    vertexBuffer.push_back({x2, y2, u, v, r, g, b, a, z_index, mode});
    vertexBuffer.push_back({x1, y2, u, v, r, g, b, a, z_index, mode});
}

static void pushText(std::vector<DynamoEngine::UIVertex>& vertexBuffer, const std::string& text, float startX,
                     float startY, const FontAtlasUtility& fontAtlas, uint32_t z_index)
{
    float cursorX = startX;
    uint32_t textMode = static_cast<uint32_t>(DynamoEngine::UIVertexMode::Textured);
    for (char c : text)
    {
        GlyphMetrics metrics = fontAtlas.GetGlyph(c);

        float tx1 = cursorX + metrics.offsetX;
        float ty1 = startY + metrics.offsetY;
        float tx2 = tx1 + metrics.width;
        float ty2 = ty1 + metrics.height;

        vertexBuffer.push_back({tx1, ty1, metrics.u1, metrics.v1, 1.0f, 1.0f, 1.0f, 1.0f, z_index, textMode});
        vertexBuffer.push_back({tx2, ty1, metrics.u2, metrics.v1, 1.0f, 1.0f, 1.0f, 1.0f, z_index, textMode});
        vertexBuffer.push_back({tx1, ty2, metrics.u1, metrics.v2, 1.0f, 1.0f, 1.0f, 1.0f, z_index, textMode});
        vertexBuffer.push_back({tx2, ty1, metrics.u2, metrics.v1, 1.0f, 1.0f, 1.0f, 1.0f, z_index, textMode});
        vertexBuffer.push_back({tx2, ty2, metrics.u2, metrics.v2, 1.0f, 1.0f, 1.0f, 1.0f, z_index, textMode});
        vertexBuffer.push_back({tx1, ty2, metrics.u1, metrics.v2, 1.0f, 1.0f, 1.0f, 1.0f, z_index, textMode});

        cursorX += metrics.advanceX;
    }
}

void Slider::buildGeometry(std::vector<DynamoEngine::UIVertex>& vertexBuffer, uint32_t z_index,
                           const FontAtlasUtility& fontAtlas)
{
    pushQuad(vertexBuffer, trackRect, ColorLibrary::Gray, z_index);
    pushQuad(vertexBuffer, knobRect, ColorLibrary::White, z_index);

    std::string slider_text = getDisplayText();
    pushText(vertexBuffer, slider_text, getX(), getY() + knobRect.h, fontAtlas, z_index);
}

void Slider::playTickSoundIfMoved(UIState& uiState)
{
    if (max_value == minValue)
        return; // avoid division by zero on an uninitialized/degenerate slider

    int currentTick =
        static_cast<int>(std::round((sliderValue - minValue) / (max_value - minValue) * NUM_SLIDER_TICKS));

    if (currentTick != lastTickIndex)
    {
        uiState.QueueSoundEffect("SliderTick");
        lastTickIndex = currentTick;
    }
}