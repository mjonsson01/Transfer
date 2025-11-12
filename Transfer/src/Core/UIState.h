#pragma once
// No SDL here. Will just include data to manage the UI state.
#include <vector>
#include "Entities/UIElement.h"
class UIState
{
    public:
        UIState();
        ~UIState();

    public:
        float getFPS() const { return fps; }
        void setFPS(float framesPerSecond) { fps = framesPerSecond; }
        // Add UI state management methods and members here

        bool getShowFPSCounter() const { return showFPSCounter; }
        void setShowFPSCounter(bool show) { showFPSCounter = show; }
    private:
        // FPS counter state
        float fps = 0.0f;
        bool showFPSCounter = false;

        std::vector<UIElement*> elements;
};