// File: Transfer/src/Systems/Rendering/GPUTypes.hpp

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

struct VelocityVectorVertex
{
    float x, y;
    float r, g, b, a;
};

// One end of a planning-mode trajectory line segment, in world space (TrajectoryLine shaders). Written once and never
// changed: the fade-in is computed in the vertex shader from drawn_at and the trace clock, so each frame only the
// new segments need uploading.
struct TrajectoryLineVertex
{
    float x, y;
    float prevX, prevY; // only differs from x, y on a line's first point: it follows the body as drawn
    float r, g, b, a;   // a = the settled alpha
    float drawn_at;     // Trajectories::clock when the point was recorded
    float fades_in;     // 1 = fades in from drawn_at, 0 = always alpha a (the bright tips, the placement preview)
};

// One end of a debug-overlay line, in world space. Like the ship sprite it also carries its previous-tick position,
// so the overlay is interpolated exactly like the things it outlines and sits right on top of them.
struct DebugLineVertex
{
    float x, y;
    float prevX, prevY;
    float r, g, b, a;
};

struct StarshipVertex
{
    float x, y;
    float prevX, prevY;
    float square_size;
    float u, v;
    float r, g, b, a;
};