// File: Transfer/src/Entities/GravitationalBody.h

#pragma once

// Custom Imports
#include "Entities/PhysicsStructures.h"

// Standard Library Imports
#include <cstdint>
#include <iostream>
#include <vector>

// To be implemented in the future for non-circular gravitational bodies

// enum class ShapeType { Circle, Polygon };

// struct PolygonShape {
//     std::vector<Vector2D> vertices;
// };

// class GravitationalBody
// {
// public:
//     ShapeType shape = ShapeType::Circle;
//     PolygonShape polygon;   // only used when shape == Polygon

//     // everything else stays the same...
// };



// struct Particle {
//     Vector2D position {0.0, 0.0};
//     Vector2D prevPosition = {0.0, 0.0};
//     Vector2D velocity = {0.0, 0.0};
//     Vector2D netForce = {0.0, 0.0};
//     Vector2D prevForce = {0.0, 0.0};
//     double mass = 0.0;
//     double radius = 0.0;
//     int clusterID; // ID of the cluster this particle belongs to. -1 if belongs to none.
//     bool collisionEnabled = true;
//     bool markedForDeletion = false;
// };

// struct GravitationalCluster {
//     // std::vector<Particle> clusterParticles;
//     std::vector<int> clusterParticleIndices; // Indices of particles in the GameState's particle vector that belong to this cluster
//     double initialRadius = 0.0;
//     double totalMass = 0.0;
//     Vector2D centerOfMassPos {0.0, 0.0};
//     Vector2D bulkVelocity {0.0, 0.0}; // Velocity of bulk
//     bool isDebrisField = false;
// };



struct GravitationalBody
{
    Vector2D position = {0.0, 0.0};
    Vector2D prevPosition = {0.0, 0.0};
    Vector2D velocity = {0.0, 0.0};
    Vector2D netForce = {0.0, 0.0};
    Vector2D prevForce = {0.0, 0.0};
    double mass = 0.0;
    double radius = 0.0;
    // int clusterID; // ID of the cluster this particle belongs to. -1 if belongs to none.
    bool collisionEnabled = true;
    bool markedForDeletion = false;

    bool isStatic = false;
    bool isPlanet = false;
    bool isFragment = false;
    bool isDust = false;
};
