// File: Transfer/src/DynamoEngine/Rendering/UIVertex.hpp

#pragma once

// Standard Library Imports
#include <cstdint>

namespace DynamoEngine
{
// How the UI fragment shader colors a vertex (must match UIElement.frag.hlsl)
enum class UIVertexMode : uint32_t
{
    None = 0,
    Solid = 1,   // flat color
    Textured = 2 // color * font-atlas alpha (text)
};

struct UIVertex
{
    float x = 0.0f, y = 0.0f;                     // screen position
    float u = 0.0f, v = 0.0f;                     // font atlas texture coordinates
    float r = 0.0f, g = 0.0f, b = 0.0f, a = 0.0f; // color and opacity 0-1
    uint32_t z_index = 0;                         // cpu-side drawing order index
    uint32_t mode = static_cast<uint32_t>(UIVertexMode::None);
};

} // namespace DynamoEngine