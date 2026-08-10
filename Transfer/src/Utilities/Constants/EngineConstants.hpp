// File: Transfer/src/Utilities/Constants/EngineConstants.hpp

#pragma once

// Hardcoded Pi
constexpr double PI = 3.14159265358979323846;

// Near-zero comparison epsilon.
constexpr double EPSILON = 1e-8;

// Load balancing max renderable bodies on screen at once
constexpr uint32_t MAX_UNIFIED_BODIES = 20000;

// Arbitrary limit to number of UI vertices
constexpr uint32_t MAX_UI_VERTICES = 65536;

// Load balancing to prevent too many particles from being instantiated
static const uint32_t MAX_LIVE_PARTICLES = 20000;

// Load balancing to prevent too many shatters from occuring in a single frame.
static const uint32_t MAX_SIMULTANEOUS_SHATTERS_PER_TICK = 20; // headroom heuristic, see handleCollisions reserve()

// Grav body max/mins
<<<<<<< HEAD
constexpr double MAX_MASS = 1e11;
constexpr double MAX_RADIUS = 300;
=======
constexpr double MAX_MASS = 1e12;
constexpr double MAX_RADIUS = 600;
>>>>>>> 35bb4b8cd1894f2d0615f49b56a79c7c5cbf9702
