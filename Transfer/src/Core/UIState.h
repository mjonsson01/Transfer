// File: Transfer/src/Core/UIState.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>

// Custom Imports
#include "Core/InputState.h"
#include "Utilities/GameSystemConstants.h"

// Standard Library Imports
#include <vector>

// class UIElement;
// struct UIElementVisibility
// {
//     bool allElementsVisible = false;
//     bool FPSVisibile = false;
// };

class UIState
{
    public:
        UIState();
        ~UIState(); 
        InputState& getMutableInputState() { return inputState; }
        const InputState& getInputState() const { return inputState; }
        float getFPS() { return framesPerSecond;}
        void setFPS(float fps) { framesPerSecond = fps;}
        bool getAllVisibility() {return allUIElementsVisible;}
        // bool getFPSVisibility() { return allUIElementsVisibility.FPSVisibile;}

    private:
        InputState inputState;
        float framesPerSecond = TARGET_FPS;
        bool allUIElementsVisible = true; // Default to true
        // UIElementVisibility allUIElementsVisibility;
        
};

// class UIState
// {
//     public:
//         UIState();
//         ~UIState();

//     public:
//         float getFPS() const { return fps; }
//         void setFPS(float framesPerSecond) { fps = framesPerSecond; }


//         // Add UI state management methods and members here
//         const std::vector<UIElement*>& getUIElements() const { return UIElements; }
        
//         // Initializing helper.
//         void addUIElement(UIElement* uielement) { UIElements.push_back(uielement); }

//         // Cleanup helper.
//         void clearUIElements() { UIElements.clear(); }

//         bool getShowFPSCounter() const { return showFPSCounter; }
//         void setShowFPSCounter(bool show) { showFPSCounter = show; }

//         InputState& getMutableInputState() { return inputState; }
//         const InputState& getInputState() const { return inputState; }

//         bool getUIElementsVisible() const { return UIElementsVisible; }
//         void invertUIElementsVisibility() { UIElementsVisible = !UIElementsVisible; }

//         bool inRect(const SDL_FRect& rect) const {
//             bool inRect = inputState.mouseCurrPosition.xVal >= rect.x &&
//                     inputState.mouseCurrPosition.xVal <= rect.x + rect.w &&
//                     inputState.mouseCurrPosition.yVal >= rect.y &&
//                     inputState.mouseCurrPosition.yVal <= rect.y + rect.h;
//             return inRect;
//         }
//         SDL_FRect getMassKnobRect() const {
//             return massKnobRect;
//         }
//         SDL_FRect getMassElementRect() const {
//             return massElementRect;
//         }
//         void setMassKnobRectPosition(Vector2D position) {
//             // Only update X coordinate, center knob on click position (subtract half knob width)
//             massKnobRect.x = position.xVal - (massKnobRect.w / 2.0f);
//             // Y coordinate stays fixed for horizontal slider
//             if (UIElements.size() > UIElementType::MASS_KNOB_INDEX) {
//                 UIElements[UIElementType::MASS_KNOB_INDEX]->setKnob(massKnobRect);
//                 std::cout<< "Updated mass knob rect to: " << massKnobRect.x << ", " << massKnobRect.y << ", " << massKnobRect.w << ", " << massKnobRect.h << std::endl;
//             }
//         };
        
//     private:
//         // FPS counter state
//         float fps = 0.0f;
//         bool showFPSCounter = false;

//         // Rectangle locations for UI element targets.
//         // SDL_FRect massTrackRect;
//         SDL_FRect massElementRect;
//         SDL_FRect massKnobRect;

//         // Input state for UI interactions
//         InputState inputState;
        
//         bool UIElementsVisible = true; // Flag to toggle visibility of UI elements

//         // Collection of UI elements
//         std::vector<UIElement*> UIElements;
// };