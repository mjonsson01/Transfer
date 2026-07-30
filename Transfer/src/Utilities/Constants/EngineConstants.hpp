// File: Transfer/src/Utilities/Constants/EngineConstants.h

#pragma once

// Global Physics dt
constexpr float PHYSICS_TIME_STEP = 1.0f / 120.0f; // 120Hz physics update rate
constexpr float INV_PHYSICS_TIME_STEP = 120.0f;
// Newtonian Gravitational constant (scaled)
// const double GRAVITATIONAL_CONSTANT = 6.67430e-7 / 2; // in m^3 kg^-1 s^-2
const double GRAVITATIONAL_CONSTANT = 1.0; //

const double MAX_MASS = 1e7;
const double MAX_RADIUS = 200;

const double MIN_BODY_BODY_ACCRETION_THRESHOLD_RATIO = 10.0;
const double MIN_BODY_PARTICLE_ACCRETION_THRESHOLD_RATIO = 10.0;

const double MAX_ACCRETION_COLLISION_SPEED = 45.0; // in px/s need to rescale later to m/s // not necessarily the
                                                   // same as the min shatter speed
const double MIN_SHATTER_SPEED = 90.0;             // in px/s need to rescale later to m/s

// Scaling factors for converting between simulation units and screen units

const double PI = 3.14159265358979323846;

const double EPSILON = 1e-8;

const uint32_t MAX_UNIFIED_BODIES = 200000;

const double ELASTIC_LOSS_FACTOR = 0.995; // need to tweak if no longer having particle particle gravity

const uint32_t MAX_UI_VERTICES = 65536;

// Fragment generation (density-factor + sunflower-spiral substitution)
static const double GOLDEN_ANGLE = 2.399963229728653;           // pi * (3 - sqrt(5)), the Vogel/sunflower spiral angle
static const double OVERLAP_MARGIN = 0.75;                      // tune by eye: fragment radius vs. exact area-match
static const double IMPACT_SKEW_GROWTH_FACTOR = 3.0;            // tune by eye: max radius growth far from impact point
static const uint32_t DEFAULT_FRAGMENT_COUNT = 800;             // fragments per shattered body, regardless of body size
static const double MUTUAL_SHATTER_MASS_RATIO_THRESHOLD = 10.0; // below this heavy/light ratio, both bodies shatter
static const double MIN_PARTICLE_PARTICLE_ACCRETION_THRESHOLD_RATIO =
    8.0; // deliberately above the ~27x (IMPACT_SKEW_GROWTH_FACTOR^3) spread skew alone can produce in one
         // explosion
static const uint32_t MAX_SIMULTANEOUS_SHATTERS_PER_TICK = 20; // headroom heuristic, see handleCollisions reserve()
static const double PARTICLE_PROMOTION_RADIUS_THRESHOLD =
    24.0; // tune by eye: radius at which an accreting particle graduates into a macro body

static const uint32_t MAX_LIVE_PARTICLES = 20000; // hard cap, see point 2
