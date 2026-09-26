// File: Transfer/src/Entities/UIElements/Buttons/Button.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "DynamoEngine/Rendering/UIVertex.hpp"
#include "Entities/UIElements/UIElement.hpp"
#include "Entities/UIElements/UIElementIdentifierEnum.hpp"
#include "Utilities/Rendering/GPUTypes.hpp"

// Standard Library Imports
#include <string>
#include <vector>

class Button : public UIElement
{
  public:
    Button();
    ~Button() = default;
    void buildGeometry(DynamoEngine::UIGeometryBuilder& builder) override;

    virtual std::string getDisplayText() const { return displayText; }
    virtual void clickMe(DynamoEngine::Vector2D positionOfEvent, UIState& uiState) override;
    double getButtonState() { return buttonSelected; }

  protected:
    SDL_FRect boundingRect;
    bool buttonSelected;
    std::string displayText;
    std::string altText;
};