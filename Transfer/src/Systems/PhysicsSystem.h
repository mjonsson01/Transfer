// File: Transfer/src/Systems/PhysicsSystem.h

#pragma once

// Custom Imports
#include "Core/GameState.h"
#include "Utilities/EngineConstants.h"
#include "Utilities/GameSystemConstants.h"
#include "Entities/GravitationalBody.h"
#include "Entities/PhysicsStructures.h"

// Standard Library Imports
#include <cmath>
#include <numeric>

class PhysicsSystem
{
    public:
        // Constructor and Destructor
        PhysicsSystem();
        ~PhysicsSystem();

    public:
        // Methods to update physics. Essentially integrates one physics frame worth of information
        void CleanUp();
        void UpdateSystemFrame(GameState& state);

    private:
        // Internal state variables for physics calculations can be added here

        // Function to calculate the resultant force vectors between two bodies
        void calculateGravForceBetweenBodies(GravitationalBody& bodyA, GravitationalBody& bodyB);
        void updateBodyVectors(GravitationalBody& body);
        void updateBoundingBox(GravitationalBody& body);
        void handleCollisions(GameState& state);
        void handleElasticCollision(GravitationalBody& bodyA, GravitationalBody& bodyB, Vector2D& unit_direction_vector);
};