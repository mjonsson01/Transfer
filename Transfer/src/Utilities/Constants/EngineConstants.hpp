// File: Transfer/src/Utilities/Constants/EngineConstants.h

#pragma once

// Global Physics dt
constexpr float PHYSICS_TIME_STEP = 1.0f / 120.0f; // 120Hz physics update rate
constexpr float INV_PHYSICS_TIME_STEP = 120.0f;
// Newtonian Gravitational constant (scaled)
// const double GRAVITATIONAL_CONSTANT = 6.67430e-7 / 2; // in m^3 kg^-1 s^-2
const double GRAVITATIONAL_CONSTANT = 1.0; //

const double MAX_MASS = 1e4;
const double MAX_RADIUS = 200;

const double MIN_BODY_BODY_ACCRETION_THRESHOLD_RATIO = 10.0;
const double MIN_BODY_PARTICLE_ACCRETION_THRESHOLD_RATIO = 10.0;

const double MAX_ACCRETION_COLLISION_SPEED = 90.0; // in px/s need to rescale later to m/s // not necessarily the
                                                   // same as the min shatter speed
const double MIN_SHATTER_SPEED = 90.0;             // in px/s need to rescale later to m/s

// Scaling factors for converting between simulation units and screen units

const double PI = 3.14159265358979323846;

const double EPSILON = 1e-8;

const uint32_t MAX_UNIFIED_BODIES = 200000;

const double ELASTIC_LOSS_FACTOR = 0.75; // 75% of Velocity retained after elastic collision

const uint32_t MAX_UI_VERTICES = 65536;