// File: Transfer/src/Utilities/Physics/FragmentLayout.hpp

#pragma once

// Custom Imports
#include "Entities/Physics/GravitationalBody.hpp"

// Standard Library Imports
#include <cstdint>
#include <vector>

// How many particles `body` breaks into: as many as fit while each stays at least MIN_PARTICLE_RADIUS, but no more
// than max_count (and always at least 1)
uint32_t survivableFragmentCount(const GravitationalBody& body, uint32_t max_count);

// Appends the particles `body` breaks into to fragments_out: a sunflower (Vogel) spiral filling the body's circle, its
// mass split evenly between them. Only the layout: velocities are left at 0, and `body` itself isn't changed.
// Real clusters and debris (PhysicsSystem::substituteWithParticles) and the cluster spawn preview both use this, so
// the preview always shows exactly where the particles will appear.
void layOutFragments(const GravitationalBody& body, uint32_t max_count, std::vector<GravitationalBody>& fragments_out);