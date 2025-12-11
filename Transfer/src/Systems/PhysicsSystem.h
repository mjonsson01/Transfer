// File: Transfer/src/Systems/PhysicsSystem.h

#pragma once

// Custom Imports
#include "Core/GameState.h"
#include "Core/UIState.h"
#include "Utilities/EngineConstants.h"
#include "Utilities/GameSystemConstants.h"
#include "Entities/GravitationalBody.h"
#include "Entities/PhysicsStructures.h"

// Standard Library Imports
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>

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
        void UpdateSystemFrame(GameState& state, UIState& UIState);

    private:
        // Internal state variables for physics calculations can be added here

        // void calculateGravForceBetweenParticles(Particle& particleA, Particle& particleB);
        // void applyInterClusterGravity(GameState& state);

        // void integrateParticleForwards(Particle& particle);
        // void integrateParticleVelocityVerlet(Particle& p);
        // void integrateParticleLeapfrog(Particle& p);

        
        // // Top-level collision handling
        // void handleCollisions(GameState& state);
        
        // void handleElasticParticleCollision(Particle& particleA, Particle& particleB);
        

        // std::vector<double> generateRandomFragmentMasses(double totalMass, int numFragments);

        // Vector2D randomDirectionVector();

        // // Helper to make sure we don't kill performance with too many bodies.
        // void manageLoad(GameState& state);


        // void createNewGravitationalCluster(UIState& UIState, GameState& state);
        // void recomputeClusterProperties(GameState& state, int clusterIndex);

        void applyGravity(PhysicsData& physicsData);
        void predictPositions(PhysicsData& physicsData);
        void buildSpatialGrid(PhysicsData& physicsData);
        void solveXPBDConstraints(PhysicsData& physicsData);
        void createParticleCluster(PhysicsData& physicsData, InputState& inputState);

};