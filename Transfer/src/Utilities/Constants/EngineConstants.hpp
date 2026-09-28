// File: Transfer/src/Utilities/Constants/EngineConstants.hpp

#pragma once

// Load balancing max renderable bodies on screen at once
constexpr uint32_t MAX_UNIFIED_BODIES = 12000; // This only balances the rendering system, which is not actually the
                                               // bottleneck. Need to fix the physics system load balancing.

// Arbitrary limit to number of UI vertices
constexpr uint32_t MAX_UI_VERTICES = 65536;

// Load balancing to prevent too many particles from being instantiated
constexpr uint32_t MAX_LIVE_PARTICLES = 20000;

constexpr uint32_t MAX_STARSHIP_VERTICES = 32; // Unknown if needed
// Grav body max/mins
constexpr double MAX_MASS = 1e11;
constexpr double MAX_RADIUS = 300;

constexpr double MIN_PARTICLE_RADIUS = 1.0;
