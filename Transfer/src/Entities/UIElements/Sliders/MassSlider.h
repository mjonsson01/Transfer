// File: Transfer/src/Entities/UIElements/Sliders/MassSlider.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Entities/UIElements/UIElement.h"
#include "Entities/UIElements/Sliders/Slider.h"
#include "Core/UIState.h"
#include "Utilities/EngineConstants.h"
#include "Utilities/GameSystemConstants.h"
#include "Utilities/Colors.h"

// Standard Library Imports
#include <string>

class MassSlider : public Slider
{
    public:
        MassSlider();
        ~MassSlider() = default;
        std::string getDisplayText() const override
        {
            return "Mass: " + std::to_string(sliderValue); // add units?
        }
            
};


// class MassSlider : public UIElement {
//     public:
//         MassSlider(); // constructor
//         virtual ~MassSlider() = default; // destructor
//         virtual void renderMe(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont) override;

// };

// class MassSlider : public UIElement {
//     public:
//         MassSlider(); // constructor
//         virtual ~MassSlider(); // destructor

//     public:
//         virtual void renderElement(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont) override; // Mass Slider specialty render method

//     public:
//         // Special Getter and Setter for Mass Value
//         float getSelectedMassValue() const { return selectedMassValue; }
//         void setSelectedMassValue(float mass) { selectedMassValue = mass; }
//         void getTrackAndKnobPositions();
//         virtual void setKnob(SDL_FRect newKnobRect) override { knobRect = constrainKnobToTrack(newKnobRect); }
//         SDL_FRect getTrackRect() const { return trackRect; }
//     private:
//         // Element State variables
//         float selectedMassValue = 0.0f; // Current mass value selected by the slider
//         float minMassValue = 0.0f;      // Minimum mass value -- revisit to allow negative masses?
//         float maxMassValue = MAX_MASS;    // Maximum mass value on slider
//         bool isSelected = false;        // Slider Interaction Status
        
//         // GUI Elements (copied from the master location in UIState  SDL_FRect massTrackRect and SDL_FRect massKnobRect)

//         const SDL_FRect trackRect = {50.0f, 3*SCREEN_HEIGHT/4, 200.0f, 20.0f}; // default position and size for the mass track, will
//         SDL_FRect knobRect;

//         // Helper method to constrain knob to track bounds
//         SDL_FRect constrainKnobToTrack(SDL_FRect proposedKnob) const {
//             SDL_FRect constrained = proposedKnob;
            
//             // Constrain left edge: knob.x >= trackRect.x
//             if (constrained.x < trackRect.x) {
//                 constrained.x = trackRect.x;
//             }
            
//             // Constrain right edge: knob.x + knob.w <= trackRect.x + trackRect.w
//             float knobRightEdge = constrained.x + constrained.w;
//             float trackRightEdge = trackRect.x + trackRect.w;
//             if (knobRightEdge > trackRightEdge) {
//                 constrained.x = trackRightEdge - constrained.w;
//             }
            
//             return constrained;
//         }
// };