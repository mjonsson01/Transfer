// File: Transfer/src/Entities/UIElements/Buttons/Button.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Entities/UIElements/UIElement.h"
#include "Entities/UIElements/UIElementIdentifierEnum.h"

// Standard Library Imports
#include <string>

class Button : public UIElement
{
  public:
    Button();
    ~Button() = default;
    void renderMe(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont) override;
    virtual std::string getDisplayText() const { return displayText; }
    virtual void clickMe(Vector2D positionOfEvent, UIState& UIState) override;
    double getButtonState() { return buttonSelected; }

  protected:
    SDL_FRect boundingRect;
    bool buttonSelected;
    std::string displayText;
    std::string altText;
};