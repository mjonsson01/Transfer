// File: Transfer/src/Utilities/Physics/FragmentLayout.cpp

#include "Utilities/Physics/FragmentLayout.hpp"

// Custom Imports
#include "DynamoEngine/Constants/GlobalConstants.hpp"
#include "Utilities/Constants/EngineConstants.hpp"
#include "Utilities/Constants/PhysicsConstants.hpp"

// Standard Library Imports
#include <algorithm>
#include <cmath>

uint32_t survivableFragmentCount(const GravitationalBody& body, uint32_t max_count)
{
    double max_survivable = (OVERLAP_MARGIN * body.radius) / MIN_PARTICLE_RADIUS;
    return std::clamp(static_cast<uint32_t>(max_survivable * max_survivable), 1u, max_count);
}

void layOutAndPopulateFragments(const GravitationalBody& body, uint32_t max_count,
                                std::vector<GravitationalBody>& fragments_out)
{
    uint32_t num_particles = survivableFragmentCount(body, max_count);

    const double R = body.radius;
    const DynamoEngine::Vector2D center = body.position;

    double density_factor = (DynamoEngine::PI * R * R) / num_particles;
    const double fragment_radius = OVERLAP_MARGIN * sqrt(density_factor / DynamoEngine::PI);
    const double particle_mass = body.mass / num_particles;

    for (uint32_t k = 0; k < num_particles; ++k)
    {
        double r_k = R * sqrt((k + 0.5) / num_particles);
        double theta_k = k * GOLDEN_ANGLE;
        DynamoEngine::Vector2D pos_k = center + DynamoEngine::Vector2D{r_k * cos(theta_k), r_k * sin(theta_k)};

        GravitationalBody p;
        p.mass = particle_mass;
        p.invMass = 1.0 / particle_mass;
        p.radius = fragment_radius;
        p.position = pos_k;
        p.previousPosition = pos_k;
        p.isFragment = true;
        p.isAccretable = body.isAccretable; // inherit body's accretable flags
        p.isCollidable = body.isCollidable; // inherit body's collision flags.
        // p.isForceState = body.isForceState; // gameplay decision, tbd.
        p.macroIdentifier = body.macroIdentifier;
        p.isParticle = true;

        fragments_out.push_back(p);
    }
}