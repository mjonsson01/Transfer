// File: Transfer/src/Utilities/EngineConstants.h

#pragma once

// Global Physics dt
constexpr float PHYSICS_TIME_STEP = 1.0f/120.0f; // 120Hz physics update rate
// constexpr float PHYSICS_TIME_STEP = 1.0f/165.0f; // 165Hz physics update rate
// Newtonian Gravitational constant
static const double GRAVITATIONAL_CONSTANT = 6.67430e-8; // in m^3 kg^-1 s^-2

static const double MAX_MASS = 1e15;

static const double ACCRETION_THRESHOLD_RATIO = 100;

// static const uint32_t MERGE_COLLISION_TIME_THRESHOLD = 2000; // in milliseconds // switch if necessary

// Max Speed to NOT cause dynamic collision.
static const double MAX_ELASTIC_COLLISION_SPEED = 200.0; // in px/s need to rescale later to m/s
// Scaling factors for converting between simulation units and screen units

static const double PI = 3.14159265358979323846;

static const int MAX_BODIES_ALLOWED = 500; // will be defined by system performance?


// Add later