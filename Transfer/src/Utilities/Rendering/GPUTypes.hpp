// File: Transfer/src/Utilities/Rendering/GPUTypes.hpp

#pragma once

// This vertex will be used and batched over the SDL_GPU hlsl calls.
struct UnifiedBodyVertex
{
    float x, y;         // 8 bytes positions
    float prevX, prevY; // 8 bytes prev positions to help interpolation
    float radius;       // 4 bytes

    // Custom View Attributes
    float log_mass; // logged mass to allow double range mass to fit within a float.
    float temperature;
    float charge;

    // Identification
    uint32_t flags; // gets all the property bools. room for 16 property bools.
    uint32_t seed;  // used to generate procedural patterns later.
};