// File: Transfer/src/Systems/PhysicsSystem.h

#pragma once

// Custom Imports
#include "Core/GameState.h"
#include "Core/UIState.h"
#include "Utilities/EngineConstants.h"
#include "Utilities/GameSystemConstants.h"
#include "Utilities/CustomMathUtilities.h"
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
        // Method to update Physics System. Handles all physics interactions and body instantiations for one physics frame
        void UpdateSystemFrame(GameState& gameState, UIState& UIState);
        // Helper method called in the destructor to clear up physics-related contents
        void CleanUp();
    
    private:
        // Add/remove any user-requested new Gravitational bodies;
        void updateGravBodyInstantiations(GameState& gameState, UIState& UIState);
        
        // Gravity Methods
        void updateGravityForSystem(GameState& gameState);                         // Gravity calculation dispatch helper
        void calculateGravity(GravitationalBody& body1, GravitationalBody& body2); // Calculate and apply gravity between two gravitational bodies
        
        // Two-step integration (a la Verlet with half-steps for stability)
        void integrateForwardsPhase1(GameState& gameState); // First half step of integration. Kicks velocity halfway and drifts position. Occurs with leftover forces from previous physics frame
        void integrateForwardsPhase2(GameState& gameState); // Second half step of integration. Kicks Velocity halfway. Occurs with newly calculated forces from current physics frame

        // Top-level collision handler. Makes decisions about the kinds of collisions encountered and dispatches to the subhandlers
        void handleCollisions(GameState& gameState);
        
        // Sub-level collision handlers
        void handleElasticCollisions(GravitationalBody& smallerBody, GravitationalBody& largerBody);                    // Handles 'bouncy' collision if the collision satisfies Engine-Constant-defined constraints
        void handleDynamicExplosionCollision(GravitationalBody& body1, GravitationalBody& body2, GameState& gameState); // Handles 'explosive' collisions that shatter the pieces if the collision satisfies Engine-Constant-defined constraints
        void handleAccretion(GravitationalBody& particle, GravitationalBody& body);                                     // Handles accretion collision events for bodies absorbing particles

        // Gravitational Body Creation Mechanisms
        void createMacroBody(GameState& gameState, InputState& inputState);          // Creates a Macro Gravitational Body with the user-defined attributes
        void createParticle(GameState& gameState, InputState& inputState);           // Creates a Particle Gravitational Body with the user-defined attriubutes
        // void createParticleCluster(Particle& originalBody, GameState& gameState); // Creates a cluster of Particle Gravitational Bodies with the user-defined attributes
        
        // Utility Functions
        void calculateTotalEnergy(GameState& gameState); // Calculates total energy of all Macro Bodies and Particles on Screen.
        void substituteWithParticles(GravitationalBody& originalBody, GameState& gameState); // Replaces a Macro Body with Particles in place to aid accretion

        // Cleanup Functions
        void cleanupParticles(GameState& gameState);   // Clears any Particles from the screen flagged as marked for deletion
        void cleanupMacroBodies(GameState& gameState); // Clears any Macro Bodies from the screen flagged as marked for deletion
};