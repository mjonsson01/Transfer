// File: Transfer/src/Utilities/Constants/EngineConstants.hpp

#pragma once

// Particle budget. Spawns are budgeted for the worst case (every macro body shattering into DEFAULT_FRAGMENT_COUNT
// particles), so physics can run in full and still never have more than this many particles.
// See PhysicsSystem::potentialParticleCount.
constexpr uint32_t MAX_LIVE_PARTICLES = 20000;

// Size of the body renderer's GPU buffers. Each macro body costs DEFAULT_FRAGMENT_COUNT of the budget but is only
// one body, so the budget also caps ALL bodies at MAX_LIVE_PARTICLES; + 1 is the spawn-preview body.
// If that's ever exceeded anyway, RenderSystem::uploadUnifiedBodies grows the buffers and logs a warning.
constexpr uint32_t INITIAL_UNIFIED_BODY_CAPACITY = MAX_LIVE_PARTICLES + 1;

// Arbitrary limit to number of UI vertices
constexpr uint32_t MAX_UI_VERTICES = 65536;

constexpr uint32_t MAX_STARSHIP_VERTICES = 32; // Unknown if needed
// Grav body max/mins
constexpr double MAX_MASS = 1e11;
constexpr double MAX_RADIUS = 300;

constexpr double MIN_PARTICLE_RADIUS = 1.0;
