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
        void calculateGravity(GravitationalBody& body1, GravitationalBody& body2);
        void calculateGravityForSmallFragments(GravitationalBody& particle, GravitationalBody& body);
        // void integrateForwards(GameState& state);
        void integrateForwards_Phase1(GameState& state);
        void integrateForwards_Phase2(GameState& state);
        void handleCollisions(GameState& state);
        void handleElasticCollisions(GravitationalBody& particle, GravitationalBody& body);
        // void handleDynamicCollision(GameState& state);
        void createPlanet(GameState& state, InputState& inputState);
        // void createPlanetaryCluster(Particle& originalBody, GameState& state);
        void calculateTotalEnergy(GameState& state);
        void substituteWithParticles(GravitationalBody& originalBody, GameState& state);
        void cleanupParticles(GameState& state);
        void cleanupMacroBodies(GameState& state);
};