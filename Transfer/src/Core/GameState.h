// File: Transfer/src/Core/GameState.h

#pragma once


// Custom Imports
#include "Entities/GravitationalBody.h"
#include "Entities/PhysicsStructures.h"
#include "Utilities/GameSystemConstants.h"
// Standard Library Imports
#include <vector>




class GameState
{
    public:
        // Constructor and Destructor
        GameState();
        ~GameState();

    public:
        // Getters and Setters for Game State
        bool IsPlaying() const { return isPlaying; }
        void SetPlaying(bool playing) { isPlaying = playing; }

        float getAlpha() const {return alpha;}
        void setAlpha(float alphaIn) {alpha = alphaIn;}

        float getTimeScaleFactor() const {if (toggleSlow) return SLOW_TIME_SCALE_FACTOR; else return REGULAR_TIME_SCALE_FACTOR;}
        // void setTimeScaleFactor(float scale) {timeScaleFactor = scale;}

        // bool getToggleSlow() const {return toggleSlow;}
        void invertToggleSlow() {toggleSlow = !toggleSlow;}

        PhysicsData physicsData; // physics data for the game
    private:
        // State variables
        bool isPlaying = false;
        
        // Frame helper vars
        float alpha = 0.0f;
        
        // float timeScaleFactor = REGULAR_TIME_SCALE_FACTOR; // default to 1.0 for standard scaling
        bool toggleSlow = false; // default to false for regular speed
        // Collection of bodies in the game 
        std::vector<GravitationalBody> bodies; 



        // New particle shit
        // std::vector<Particle> particles; // All particles in the sim
        // std::vector<GravitationalCluster> clusters; // All clusters in the sim
        // std::vector<XPBDSoftParticle> global_particles;
        // std::vector<XPBDConstraints> global_constraints;
};

