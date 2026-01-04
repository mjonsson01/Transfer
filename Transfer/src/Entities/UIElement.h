// File: Transfer/src/Entities/UIElement.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Utilities/Colors.h"
#include "Core/UIState.h"

// Standard Library Imports
#include <string>


class UIElement{
    public:
        UIElement();
        virtual ~UIElement();

        virtual void renderElement(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont) // Default renderElement implementation
        {
            // Default implementation (can be empty)
            return;
        }; // will be overloaded by derived classes

    public:
        // Getters and Setters
        void setPosition(float x, float y) { posX = x; posY = y; }
        void setSize(float w, float h) { width = w; height = h; }
        float getX() const { return posX; }
        float getY() const { return posY; }
        float getWidth() const { return width; }
        float getHeight() const { return height; }
        float getScaleFactor() const { return scaleFactor; }
        void setScaleFactor(float scale) { scaleFactor = scale; }

    private:
        float posX = 0;
        float posY = 0;
        float width = 0;
        float height = 0;
        float scaleFactor = 1.0f;
};

// Derived UI Element Classes
// (Each derived class should have its own header file)

// class SimulationSpeedSlider : public UIElement
// {
// };
// class VelocityVectorToggle : public UIElement
// {
// };

// class GravityToggle : public UIElement
// {
// };

// class MassSlider : public UIElement
// {
// };

// class PauseMenu : public UIElement
// {
// };

// class SelectionCheckbox : public UIElement
// {
// };