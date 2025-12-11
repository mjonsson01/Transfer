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
#include <algorithm>

#include <iostream>

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
        void integrateForwards(GravitationalBody& body);
        void updateBoundingBox(GravitationalBody& body);
        
        // Top-level collision handling
        void handleCollisions(GameState& state);
        
        // Simple Collision handling helpers
        void handleElasticCollision(GravitationalBody& bodyA, GravitationalBody& bodyB);
        
        // Collision handling helpers
        // void handleDynamicCollision(GravitationalBody& bodyA, GravitationalBody& bodyB, GameState& state, std::vector<GravitationalBody>& newFragments);
        void handleDynamicCollision(GravitationalBody& bodyA, GravitationalBody& bodyB, GameState& state);
        // bool isOverlappingAny(const GravitationalBody& bodyA, const GameState& state);

        std::vector<double> generateRandomFragmentMasses(double totalMass, int numFragments);
        // void resolvePenetration(GravitationalBody& bodyA, GravitationalBody& bodyB);
        Vector2D randomDirectionVector();

        // Helper to make sure we don't kill performance with too many bodies.
        void manageLoad(GameState& state);

        // void handleAccretion(GravitationalBody& smallerM, GravitationalBody& largerM );
};