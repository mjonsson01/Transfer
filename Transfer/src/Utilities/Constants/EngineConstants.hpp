// File: Transfer/src/Utilities/Constants/EngineConstants.h

#pragma once

// Global Physics dt
constexpr float PHYSICS_TIME_STEP = 1.0f / 120.0f; // 120Hz physics update rate

// Newtonian Gravitational constant (scaled)
static const double GRAVITATIONAL_CONSTANT = 6.67430e-7 / 5; // in m^3 kg^-1 s^-2

static const double MAX_MASS = 1e16;
static const double MAX_RADIUS = 400;

static const double MIN_BODY_BODY_ACCRETION_THRESHOLD_RATIO = 10.0;
static const double MIN_BODY_PARTICLE_ACCRETION_THRESHOLD_RATIO = 10.0;

static const double MAX_ACCRETION_COLLISION_SPEED = 90.0; // in px/s need to rescale later to m/s // not necessarily the
                                                          // same as the min shatter speed
static const double MIN_SHATTER_SPEED = 90.0;             // in px/s need to rescale later to m/s

// Scaling factors for converting between simulation units and screen units
// TBI

static const double PI = 3.14159265358979323846;

static const double EPSILON = 1e-8;

static const uint32_t MAX_UNIFIED_BODIES = 200000;

static const double ELASTIC_LOSS_FACTOR = 0.85;       // 85% of Velocity retained after elastic collision
static const double STATIC_ELASTIC_LOSS_FACTOR = 1.0; // 85% of Velocity retained after elastic collision

static const uint32_t MAX_UI_VERTICES = 65536;

static const double GOLDEN_ANGLE = 2.399963229728653; // pi(3-sqrt(5))

static const double OVERLAP_MARGIN = 1.2; // test for now

static const double CHUNK_MASS_FRACTION = 0.6;    // fraction of the body's mass/area given to big chunks vs dust
static const double CHUNK_DENSITY_FACTOR = 300.0; // coarse divisor -> yields "a handful" of chunks
static const double CHUNK_OVERLAP_MARGIN = 1.0;   // chunks are allowed visible gaps between them, that's the look
static const double DUST_OVERLAP_MARGIN = 1.2;    // dust must fully tile to hide gaps (this is your old OVERLAP_MARGIN)
static const double IMPACT_SKEW_GROWTH_FACTOR =
    4.5; // far-side fragments can grow up to this multiple of their base radius