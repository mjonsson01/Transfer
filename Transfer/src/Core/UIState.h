// File: Transfer/src/Core/UIState.h

#pragma once

// No SDL here. Will just include data to manage the UI state.

// Custom Imports
// #include "Entities/UIElement.h"

// Standard Library Imports
#include <vector>

// Forward declaration of UI Element class. Full include leads to circular dependency, causing failures in compilation.
class UIElement;

class UIState
{
    public:
        UIState();
        ~UIState();

    public:
        float getFPS() const { return fps; }
        void setFPS(float framesPerSecond) { fps = framesPerSecond; }


        // Add UI state management methods and members here
        std::vector<UIElement*> getUIElements() const { return elements; }
        // Initializing helper.
        void addUIElement(UIElement* element) { elements.push_back(element); }
        // Cleanup helper.
        void clearUIElements() { elements.clear(); }

        bool getShowFPSCounter() const { return showFPSCounter; }
        void setShowFPSCounter(bool show) { showFPSCounter = show; }
    private:
        // FPS counter state
        float fps = 0.0f;
        bool showFPSCounter = false;

        std::vector<UIElement*> elements;
};