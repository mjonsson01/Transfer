// File: Transfer/src/Entities/UIElements/UIElement.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>


// Custom Imports
#include "Utilities/Colors.h"
#include "Entities/UIElements/UIElementTypeEnum.h"
#include "Core/UIState.h"
#include "Entities/MathStructures.h"

// Standard Library Imports
#include <string>


class UIElement
{
    public:
        UIElement();
        virtual ~UIElement();
        virtual void renderMe(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont)
        {
            return;
        }
        virtual void updateMe(Vector2D positionOfEvent, double& returnedElementValue) {}; //Default does nothing
        void setPosition(float x, float y)
        {
            posX = x;
            posY = y;
        }
        float getX() const {return posX;}
        float getY() const {return posY;}
        void setVisibility(bool desiredVisibility) { visible = desiredVisibility;}
        UIElementType isInDeadZone(const Vector2D& positionToCheck);
    private: 
        float posX = 0;
        float posY = 0;
        bool visible = false;
    protected:
        SDL_FRect hotZoneRect;
        UIElementType UIElementTypeIdentifier;
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