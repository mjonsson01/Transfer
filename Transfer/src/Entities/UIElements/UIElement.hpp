// File: Transfer/src/Entities/UIElements/UIElement.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Core/UIState.hpp"
#include "Entities/UIElements/UIElementIdentifierEnum.hpp"
#include "Utilities/Math/Vector2D.hpp"
#include "Utilities/Rendering/Colors.hpp"
#include "Utilities/Rendering/FontAtlasUtility.hpp"
#include "Utilities/Rendering/GPUTypes.hpp"

// Standard Library Imports
#include <string>
#include <vector>

class UIElement
{
  public:
    UIElement();
    virtual ~UIElement();
    virtual void slideMe(Vector2D positionOfEvent, double& returnedElementValue, UIState& UIState) {
    }; // Default does nothing
    virtual void clickMe(Vector2D positionOfEvent, UIState& UIState) {}; // Default does nothing
    void setPosition(float x, float y)
    {
        posX = x;
        posY = y;
    }
    float getX() const { return posX; }
    float getY() const { return posY; }
    void setVisibility(bool desiredVisibility) { visible = desiredVisibility; }
    bool isVisible() const { return visible; }
    UIElementIdentifier checkAndReturnIfHit(const Vector2D& positionToCheck);
    UIElementIdentifier getUIElementID() const { return UIElementID; }
    virtual void buildGeometry(std::vector<UIElementVertex>& vertexBuffer, uint32_t zIndex,
                               const FontAtlasUtility& fontAtlas) {};    // Default does nothing
    virtual void updateMe(UIState& UIState) {};                          // Default does nothing
    virtual void updateLayout(float windowWidth, float windowHeight) {}; // Default does nothing

  private:
    float posX = 0;
    float posY = 0;
    bool visible = false;

  protected:
    SDL_FRect hotZoneRect;
    UIElementIdentifier UIElementID;
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