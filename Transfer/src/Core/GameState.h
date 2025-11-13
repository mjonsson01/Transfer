// File: Transfer/src/Core/GameState.h

#pragma once


// Custom Imports
#include "Entities/GravitationalBody.h"

// Standard Library Imports
#include <vector>

struct InputState
{
	// Body instantiation control vars
	float selectedRadius = 0.0;
	double selectedMass = 0.0;
	bool bodySelectionValidity = false;
};



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

        void updateSelectedRadius(float radius) { inputState.selectedRadius = radius; }
        float getSelectedRadius() const { return inputState.selectedRadius; }
        
        void updateSelectedMass(double mass) { inputState.selectedMass = mass; }
        double getSelectedMass() const { return inputState.selectedMass; }
        
        void setBodySelectionValidity(bool isValid) { inputState.bodySelectionValidity = isValid; }
        bool isBodySelectionValid() const { return inputState.bodySelectionValidity; }

        const std::vector<GravitationalBody>& getBodies() const { return bodies; }
        std::vector<GravitationalBody>& getBodiesMutable() {return bodies; }
        void addBody(const GravitationalBody& body) { bodies.push_back(body); }
        void clearBodies() { bodies.clear(); }

        float getAlpha() const {return alpha;}
        void setAlpha(float alphaIn) {alpha = alphaIn;}

    private:
        // State variables
        bool isPlaying = false;
        

        // Frame helper vars
        float alpha = 0.0f;
        
        // Collection of bodies in the game 
        std::vector<GravitationalBody> bodies;
        // Input related state variables
        InputState inputState;
};

