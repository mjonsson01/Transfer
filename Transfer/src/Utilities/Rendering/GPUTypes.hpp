// File: Transfer/src/Utilities/Rendering/GPUTypes.hpp

#pragma once

// This vertex will be used and batched over the SDL_GPU hlsl calls.
struct UnifiedBodyVertex
{
    float x, y;         // 8 bytes positions
    float prevX, prevY; // 8 bytes prev positions to help interpolation
    float radius;       // 4 bytes

    // Custom View Attributes
    float logMass; // logged mass to allow double range mass to fit within a float.
    float temperature;
    float charge;

    // Identification
    uint32_t flags; // gets all the property bools. room for 16 property bools.
    uint32_t seed;  // used to generate procedural patterns later.
};

struct TwinklingStarVertex
{
    float x, y;
    float radius;

    float alpha;
    float twinkleSpeed;

    uint32_t seed;
};

struct UIElementVertex
{
    float x, y;       // Screen coords
    float u, v;       // Texture coordinates
    float r, g, b, a; // Color and Opacity
    uint32_t zIndex;
    uint32_t mode;
};