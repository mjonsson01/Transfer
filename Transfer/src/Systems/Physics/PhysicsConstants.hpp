// File: Transfer/src/Systems/Physics/PhysicsConstants.hpp

#pragma once
#include "Constants/EngineConstants.hpp"

// Global Physics DT
constexpr float PHYSICS_TIME_STEP = 1.0f / 120.0f; // 120Hz physics update rate
constexpr float INV_PHYSICS_TIME_STEP = 120.0f;
// Most physics steps one frame may run to catch up. Normal play needs 2 per frame (60 FPS), 4 at 2x speed.
// Time beyond this is dropped (the simulation briefly runs slower) instead of piling up into a spiral of death.
constexpr int MAX_PHYSICS_STEPS_PER_FRAME = 8;

// --- Trajectory prediction (planning mode) ---
constexpr float TRAJECTORY_TRACE_SPEED = 15.0f;          // game seconds predicted per real second: how fast lines grow
constexpr int MAX_TRACE_TICKS_PER_FRAME = 64;            // never more prediction ticks than this in one frame
constexpr double MAX_TRACE_MILLISECONDS_PER_FRAME = 4.0; // ...nor more time: heavy scenes trace slower, not choppier
constexpr double MAX_PREVIEW_MILLISECONDS_PER_FRAME = 4.0; // the placement preview's own share (it has no speed limit)
constexpr int TRAJECTORY_TICKS_PER_POINT = 4;            // record a point every 4 ticks (30 points per game second)
constexpr float TRAJECTORY_HORIZON_SECONDS = 180.0f;     // how far ahead the lines reach
constexpr float TRAJECTORY_TIME_MARK_SECONDS = 5.0f;     // a time mark on every line every 5 game seconds

// Grav constant, default to 1, scale other values
constexpr double GRAVITATIONAL_CONSTANT = 1.0;

// Collision Behavior Management constants
constexpr double MIN_BODY_BODY_ACCRETION_THRESHOLD_RATIO = 10.0;     // 10x mass ratio allows for accretion
constexpr double MIN_BODY_PARTICLE_ACCRETION_THRESHOLD_RATIO = 10.0; // 10x mass ratio allows for accretion
constexpr double ELASTIC_LOSS_FACTOR = 0.994;                        // Keep 99.4% of elastic collision energy.
constexpr double MUTUAL_SHATTER_MASS_RATIO_THRESHOLD = 10.0; // below this heavy/light ratio, both bodies shatter
constexpr double MAX_ACCRETION_COLLISION_SPEED = 160000.0;   // in px/s (in world space) need to rescale later to m/s
constexpr double MIN_SHATTER_SPEED = 400.0;                  // in px/s (in world space) need to rescale later to m/s
// How bouncy the player's ship is: the share of its closing speed it keeps after hitting something (0 = stops dead,
// 1 = no energy lost). 0.2 keeps a fifth of the speed, 4% of the impact energy: a dull thud. The rest is the
// impact energy that damage will be calculated from.
constexpr double SHIP_RESTITUTION = 0.2;

// Inactive for now, used if we reenable particle particle accretion and promotion
constexpr double MIN_PARTICLE_PARTICLE_ACCRETION_THRESHOLD_RATIO = 8.0;
constexpr double PARTICLE_PROMOTION_RADIUS_THRESHOLD =
    24.0; // Globe frame radius at which a particle promotes into a shatterable macro body

// Fragment generation constants (density-factor + sunflower-spiral substitution)
constexpr double GOLDEN_ANGLE = 2.399963229728653; // pi * (3 - sqrt(5)), the Vogel/sunflower spiral angle
constexpr double OVERLAP_MARGIN = 0.75;            // fragment radius vs. exact area-match
constexpr double IMPACT_SKEW_GROWTH_FACTOR = 3.0;  // max radius growth multiple far from impact point
constexpr uint32_t DEFAULT_FRAGMENT_COUNT = 800;   // fragments per shattered body, regardless of body size
// A cluster's spawn preview draws up to DEFAULT_FRAGMENT_COUNT particles, so the renderer reserves that much room
static_assert(DEFAULT_FRAGMENT_COUNT <= MAX_PREVIEW_BODIES, "the cluster preview wouldn't fit MAX_PREVIEW_BODIES");

// The trajectory line buffers must hold the worst case (2 vertices per segment or mark). Traced lines: every planet
// the particle budget allows, plus the ship, each to the full horizon. The per-frame buffer: the placement preview
// (one full line with its marks), the traced lines' time marks, and the bright tips (allowing up to 8 segments each).
// Raising the horizon or the point rate without raising the limits in EngineConstants.hpp stops the build here.
constexpr uint32_t MAX_TRACED_TRAJECTORY_PATHS = MAX_LIVE_PARTICLES / DEFAULT_FRAGMENT_COUNT + 1;
constexpr uint32_t TRAJECTORY_SEGMENTS_PER_PATH =
    static_cast<uint32_t>(TRAJECTORY_HORIZON_SECONDS * INV_PHYSICS_TIME_STEP) / TRAJECTORY_TICKS_PER_POINT;
constexpr uint32_t TRAJECTORY_MARKS_PER_PATH =
    static_cast<uint32_t>(TRAJECTORY_HORIZON_SECONDS / TRAJECTORY_TIME_MARK_SECONDS) + 1;
static_assert(MAX_TRACED_TRAJECTORY_PATHS * 2 * TRAJECTORY_SEGMENTS_PER_PATH <= MAX_TRAJECTORY_LINE_VERTICES,
              "traced trajectory lines wouldn't fit MAX_TRAJECTORY_LINE_VERTICES: raise it (EngineConstants.hpp)");
static_assert(2 * (TRAJECTORY_SEGMENTS_PER_PATH + TRAJECTORY_MARKS_PER_PATH) +
                      MAX_TRACED_TRAJECTORY_PATHS * 2 * TRAJECTORY_MARKS_PER_PATH +
                      (MAX_TRACED_TRAJECTORY_PATHS + 1) * 2 * 8 <=
                  MAX_TRAJECTORY_FRAME_VERTICES,
              "preview + time marks + tips wouldn't fit MAX_TRAJECTORY_FRAME_VERTICES: raise it (EngineConstants.hpp)");
