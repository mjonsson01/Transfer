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

        // void updateSelectedRadius(float radius) { inputState.selectedRadius = radius; }
        // float getSelectedRadius() const { return inputState.selectedRadius; }
        
        // void updateSelectedMass(double mass) { inputState.selectedMass = mass; }
        // double getSelectedMass() const { return inputState.selectedMass; }
        
        // void setBodySelectionValidity(bool isValid) { inputState.bodySelectionValidity = isValid; }
        // bool isBodySelectionValid() const { return inputState.bodySelectionValidity; }



        const std::vector<GravitationalBody>& getParticles() const {return particles;}
        std::vector<GravitationalBody>& getParticlesMutable() {return particles;}

        const std::vector<GravitationalBody>& getMacroBodies() const {return macroBodies;}
        std::vector<GravitationalBody>& getMacroBodiesMutable() {return macroBodies;}
        

        float getAlpha() const {return alpha;}
        void setAlpha(float alphaIn) {alpha = alphaIn;}

        float getTimeScaleFactor() const {if (toggleSlow) return SLOW_TIME_SCALE_FACTOR; else return REGULAR_TIME_SCALE_FACTOR;}
        // void setTimeScaleFactor(float scale) {timeScaleFactor = scale;}

        bool getToggleSlow() const {return toggleSlow;}
        void invertToggleSlow() {toggleSlow = !toggleSlow;}


    private:
        // State variables
        bool isPlaying = false;
        

        // Frame helper vars
        float alpha = 0.0f;
        
        // float timeScaleFactor = REGULAR_TIME_SCALE_FACTOR; // default to 1.0 for standard scaling
        bool toggleSlow = false; // default to false for regular speed

        std::vector<GravitationalBody> macroBodies;
        std::vector<GravitationalBody> particles;
};

