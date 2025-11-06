#pragma once

// Global Physics dt
constexpr float PHYSICS_TIME_STEP = 1.0f/120.0f; // 120Hz physics update rate
// Newtonian Gravitational constant
static const double GRAVITATIONAL_CONSTANT = 6.67430e-11; // in m^3 kg^-1 s^-2

static const double MAX_MASS = 1e15;
// Scaling factors for converting between simulation units and screen units

// Add later