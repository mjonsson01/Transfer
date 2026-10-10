// File: Transfer/src/Player/Starship.hpp

#pragma once

// Standard includes
#include <iostream>
#include <vector>

// Custom includes
#include "Core/DEPRECATED_InputState.hpp"
#include "Core/UIState.hpp"
#include "DynamoEngine/Math/Vector2.hpp"
#include "DynamoEngine/Rendering/UIVertex.hpp"
#include "Utilities/Constants/PhysicsConstants.hpp"
#include "Utilities/Rendering/GPUTypes.hpp"

class Starship
{
  public:
    Starship();
    ~Starship();

    // --- Physics: PhysicsSystem calls these every tick, in this order ---
    void applyRotation(UIState& uiState);
    // Reads W/S and sets the thrust acceleration (along the nose) for this tick
    void applyThrust(UIState& uiState);
    // Velocity Verlet, the same two halves the bodies use. Phase 1: half kick with last tick's acceleration, then drift
    void applyVelocityVerletPhase1();
    // Called between the two phases, once the planets' pull at the ship's NEW position is known
    void setGravityForce(const DynamoEngine::Vector2D& gravitational_force);
    // Phase 2: the other half kick, with the new acceleration (gravity + thrust)
    void applyVelocityVerletPhase2();

    // The ship's centre: its centre of mass, and the middle of its sprite and hitbox
    DynamoEngine::Vector2D center() const;
    // Half the side of the ship's square; gravity uses it as the ship's size when softening the pull up close
    double halfSize() const { return m_ship_size * 0.5; }
    double mass() const { return m_mass; }
    // --- Collisions: PhysicsSystem::handleShipCollisions uses these ---
    DynamoEngine::Vector2D velocity() const { return m_velocity; }
    // A sudden change in momentum: the velocity changes by impulse / mass (a heavier ship changes less)
    void applyImpulse(const DynamoEngine::Vector2D& impulse);
    // Moves the ship without touching its velocity (pushing it out of something it overlaps)
    void moveBy(const DynamoEngine::Vector2D& offset);
    // Distance from the centre to the hitbox's farthest corner: nothing farther than this (plus its own radius)
    // can be touching the ship
    double boundingRadius() const;
    // The kinetic energy the ship's collisions absorbed during the last tick (0 if it hit nothing). Damage will be
    // calculated from this.
    void setLastImpactEnergy(double energy) { m_last_impact_energy = energy; }
    double lastImpactEnergy() const { return m_last_impact_energy; }

    void buildGeometry(std::vector<StarshipVertex>& starshipVertexBuffer);
    DynamoEngine::Vector2D getPointingVector();
    // The hitbox: a convex polygon in world space at the ship's current position and rotation, corners in outline
    // order (what DynamoEngine::circleVsConvexPolygon expects)
    std::vector<DynamoEngine::Vector2D> collisionPolygon() const;
    // How far the ship moved during the last physics tick (position - previous position). Subtracting it from a
    // world point gives that point one tick ago, which the renderer needs to interpolate it smoothly.
    DynamoEngine::Vector2D movementThisTick() const { return m_position - m_prev_position; }

  private:
    // Turns an offset from the ship's centre by the ship's rotation. The sprite and the hitbox both use this,
    // so they can never disagree about which way the ship is facing.
    DynamoEngine::Vector2D rotatedOffset(DynamoEngine::Vector2D offset) const;

    // The hitbox traced from FutureShip.png (its outline's convex hull simplified to 6 corners, covering 96.4% of it),
    // in the image's pixels: (0, 0) = top-left of the 2048 x 2048 image, y down, clockwise on screen.
    // Re-trace these if the sprite is redrawn.
    static constexpr double HITBOX_IMAGE_SIZE_PX = 2048.0;
    static constexpr DynamoEngine::Vector2D HITBOX_CORNERS_PX[6] = {
        {32.0, 1464.0},   // left wingtip
        {959.0, 366.0},   // nose, left side
        {1094.0, 366.0},  // nose, right side
        {2021.0, 1464.0}, // right wingtip
        {1242.0, 1744.0}, // tail fin, right corner
        {811.0, 1744.0},  // tail fin, left corner
    };
    // The ship's mass without fuel. Light next to planets (the mass slider goes up to 10^10), comparable to debris.
    static constexpr double DRY_MASS = 1000.0;
    double m_mass = DRY_MASS; // just the dry mass for now; fuel will be added on top (rocket equation)
    double m_last_impact_energy = 0.0;

    // Speed gained per second of full thrust (px/s^2). 1200 indicates "+10 px/s every tick" at 120 ticks/s.
    static constexpr double THRUST_ACCELERATION = 1200.0;
    DynamoEngine::Vector2D m_thrust_acceleration = {}; // from the controls, set every tick by applyThrust
    DynamoEngine::Vector2D m_acceleration = {};        // gravity + thrust, as of the end of the last tick
    DynamoEngine::Vector2D m_velocity = {};
    DynamoEngine::Vector2D m_position = {};
    double m_rotation = {};
    DynamoEngine::Vector2D m_prev_position = {};
    float m_ship_size = {};
};