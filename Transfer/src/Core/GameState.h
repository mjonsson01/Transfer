// File: Transfer/src/Core/GameState.h

#pragma once


// Custom Imports
#include "Entities/GravitationalBody.h"
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

        // Getters for Particles with a mutable and nonmutable version to respect System hierarchies
        const std::vector<GravitationalBody>& getParticles() const {return particles;}
        std::vector<GravitationalBody>& getParticlesMutable() {return particles;}

        // Getters for Bodies with a mutable and nonmutable version to respect System hierarchies
        const std::vector<GravitationalBody>& getMacroBodies() const {return macroBodies;}
        std::vector<GravitationalBody>& getMacroBodiesMutable() {return macroBodies;}
        
        // Getter and setter for the alpha rendering variable
        float getAlpha() const {return alpha;}
        void setAlpha(float alphaIn) {alpha = alphaIn;}

        // Getter for the current time scale factor
        float getTimeScaleFactor() const {if (toggleSlow && !toggleFast) return SLOW_TIME_SCALE_FACTOR; else if (!toggleSlow && toggleFast) return FAST_TIME_SCALE_FACTOR; else return REGULAR_TIME_SCALE_FACTOR;}

        bool getToggleSlow() const {return toggleSlow;}
        void invertToggleSlow() {toggleSlow = !toggleSlow;}

        bool getToggleFast() const {return toggleFast;}
        void invertToggleFast() {toggleFast = !toggleFast;}


    private:
        // State variables
        bool isPlaying = false;
        

        // Frame helper vars
        float alpha = 0.0f;
        
        bool toggleSlow = false; // default to false for regular speed
        bool toggleFast = false; // default to false for regular speed

        // Database for all the Macro Bodies
        std::vector<GravitationalBody> macroBodies;
        // Database for all the Particle Bodies
        std::vector<GravitationalBody> particles;
};

