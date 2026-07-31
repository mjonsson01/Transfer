// File: Transfer/src/Utilities/Constants/PhysicsConstants.hpp

#pragma once
#include "Utilities/Constants/EngineConstants.hpp"


// Global Physics DT
constexpr float PHYSICS_TIME_STEP = 1.0f / 120.0f; // 120Hz physics update rate
constexpr float INV_PHYSICS_TIME_STEP = 120.0f;

// Grav constant, default to 1, scale other values
constexpr double GRAVITATIONAL_CONSTANT = 1.0;

// Collision Behavior Management constants
constexpr double MIN_BODY_BODY_ACCRETION_THRESHOLD_RATIO = 10.0; // 10x mass ratio allows for accretion
constexpr double MIN_BODY_PARTICLE_ACCRETION_THRESHOLD_RATIO = 10.0; // 10x mass ratio allows for accretion
constexpr double ELASTIC_LOSS_FACTOR = 0.994; // Keep 99.4% of elastic collision energy. 
constexpr double MUTUAL_SHATTER_MASS_RATIO_THRESHOLD = 10.0; // below this heavy/light ratio, both bodies shatter
constexpr double MAX_ACCRETION_COLLISION_SPEED = 1600000.0; // in px/s (in world space) need to rescale later to m/s
constexpr double MIN_SHATTER_SPEED = 120.0;             // in px/s (in world space) need to rescale later to m/s

// Inactive for now, used if we reenable particle particle accretion and promotion
constexpr double MIN_PARTICLE_PARTICLE_ACCRETION_THRESHOLD_RATIO = 8.0;
constexpr double PARTICLE_PROMOTION_RADIUS_THRESHOLD = 24.0; // Globe frame radius at which a particle promotes into a shatterable macro body


// Fragment generation constants (density-factor + sunflower-spiral substitution)
constexpr double GOLDEN_ANGLE = 2.399963229728653;           // pi * (3 - sqrt(5)), the Vogel/sunflower spiral angle
constexpr double OVERLAP_MARGIN = 0.75;                      // fragment radius vs. exact area-match
constexpr double IMPACT_SKEW_GROWTH_FACTOR = 3.0;            // max radius growth multiple far from impact point
constexpr uint32_t DEFAULT_FRAGMENT_COUNT = 800;             // fragments per shattered body, regardless of body size

