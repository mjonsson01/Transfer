// File: Transfer/src/Entities/UIElements/MassSlider.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Entities/UIElement.h"
#include "Core/UIState.h"
#include "Utilities/EngineConstants.h"
#include "Utilities/Colors.h"

class MassSlider : public UIElement {
    public:
        MassSlider();
        virtual ~MassSlider();

    public:
        virtual void renderElement(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont) override;

        public:
        // Special Getter and Setter for Mass Value
        float getSelectedMassValue() const { return selectedMassValue; }
        void setSelectedMassValue(float mass) { selectedMassValue = mass; }
        // void handleInputEvent(const SDL_Event& event, UIState& UIState); // Handled in InputSystem. Events should be clear from this.
    private:
        // state variables
        float selectedMassValue = 0.0f; // Current mass value selected by the slider
        float minMassValue = 0.0f;      // Minimum mass value -- revisit to allow negative masses?
        float maxMassValue = MAX_MASS;    // Maximum mass value
        bool isSelected = false;        // Slider Interaction Status
        // GUI Elements
        SDL_FRect trackRect;
        SDL_FRect knobRect;
};