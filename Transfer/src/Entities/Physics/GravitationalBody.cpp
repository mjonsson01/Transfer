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
    v.logMass = (mass < 0) ? -logMag : logMag;

    v.temperature = static_cast<float>(temperature);
    v.charge = 0.0f; // Placeholder

    // Pack flags
    uint32_t f = 0;
    if (visible)
        f |= (1 << 0);
    if (isMacro)
        f |= (1 << 1);
    if (isParticle)
        f |= (1 << 2);
    if (isMacroGhost)
        f |= (1 << 3);
    if (isAccretable)
        f |= (1 << 4);
    if (isCollidable)
        f |= (1 << 5);
    if (isShatterable)
        f |= (1 << 6);
    if (isBounce)
        f |= (1 << 7);
    if (isDust)
        f |= (1 << 8);
    if (isStatic)
        f |= (1 << 9);
    if (isGas)
        f |= (1 << 10);
    if (isFragment)
        f |= (1 << 11);
    if (isGravStar)
        f |= (1 << 12);
    if (isPlanet)
        f |= (1 << 13);
    if (isMoon)
        f |= (1 << 14);
    if (isPreview)
        f |= (1 << 15);

    v.flags = f;

    v.seed = static_cast<uint32_t>(macroIdentifier != -1 ? macroIdentifier : 0);
    return v;
}