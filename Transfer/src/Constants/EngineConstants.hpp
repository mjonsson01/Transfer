// File: Transfer/src/Constants/EngineConstants.hpp

#pragma once

// Particle budget. Spawns are budgeted for the worst case (every macro body shattering into DEFAULT_FRAGMENT_COUNT
// particles), so physics can run in full and still never have more than this many particles.
// See PhysicsSystem::potentialParticleCount.
constexpr uint32_t MAX_LIVE_PARTICLES = 20000;

// The most see-through preview bodies drawn at once: a previewed cluster shows all of its particles, at most
// DEFAULT_FRAGMENT_COUNT (PhysicsConstants.hpp checks that it fits); a previewed planet is just 1.
constexpr uint32_t MAX_PREVIEW_BODIES = 800;

// Size of the body renderer's GPU buffers. Each macro body costs DEFAULT_FRAGMENT_COUNT of the budget but is only
// one body, so the budget also caps ALL bodies at MAX_LIVE_PARTICLES; on top of that, room for the spawn preview.
// If that's ever exceeded anyway, RenderSystem::uploadUnifiedBodies grows the buffers and logs a warning.
constexpr uint32_t INITIAL_UNIFIED_BODY_CAPACITY = MAX_LIVE_PARTICLES + MAX_PREVIEW_BODIES;

// Arbitrary limit to number of UI vertices
constexpr uint32_t MAX_UI_VERTICES = 65536;

constexpr uint32_t MAX_STARSHIP_VERTICES = 32; // Unknown if needed
// Debug overlay (F3): 2 vertices per line. Each planet outline is 32 lines (the budget allows at most 25 planets),
// the ship's hitbox 6 lines, so this leaves plenty of room for future debug drawings
constexpr uint32_t MAX_DEBUG_LINE_VERTICES = 8192;

// Planning mode's trajectory lines, 2 vertices per segment. PhysicsConstants.hpp checks the worst cases fit, so
// tuning the horizon can't silently cut lines off.
// Traced lines: written once, only new segments uploaded each frame (26 lines x 180 s x 30 points/s = ~281k, 11 MB)
constexpr uint32_t MAX_TRAJECTORY_LINE_VERTICES = 294912;
// Rebuilt every frame: the placement preview (one line), the time marks (sized on screen, so they follow the zoom)
// and the bright tips
constexpr uint32_t MAX_TRAJECTORY_FRAME_VERTICES = 16384;

// Grav body max/mins
constexpr double MAX_MASS = 1e11;
constexpr double MAX_RADIUS = 700;

constexpr double MIN_PARTICLE_RADIUS = 1.0;
