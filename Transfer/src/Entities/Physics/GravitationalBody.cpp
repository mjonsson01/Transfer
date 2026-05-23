#include "Entities/Physics/GravitationalBody.hpp"

UnifiedBodyVertex GravitationalBody::toUnifiedVertex() const
{
    UnifiedBodyVertex v;
    v.x = static_cast<float>(position.xVal);
    v.y = static_cast<float>(position.yVal);
    v.prevX = static_cast<float>(previousPosition.xVal);
    v.prevY = static_cast<float>(previousPosition.yVal);
    v.radius = static_cast<float>(radius);

    // Symmetric Log Mass logic
    float logMag = static_cast<float>(std::log10(std::max(1.0, std::abs(mass))));
    v.log_mass = (mass < 0) ? -logMag : logMag;

    v.temperature = static_cast<float>(temperature);
    v.charge = 0.0f; // Placeholder

    // Pack flags
    uint32_t f = 0;
    if (isMacro)
        f |= (1 << 0);
    if (isParticle)
        f |= (1 << 1);
    if (visible)
        f |= (1 << 2);
    // ... add others as needed
    v.flags = f;

    v.seed = static_cast<uint32_t>(macroIdentifier != -1 ? macroIdentifier : 0);
    return v;
}