// File: Transfer/src/Entities/UIElements/Overlay/FPSCounter.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Core/UIState.h"
#include "Entities/UIElements/UIElement.h"
#include "Entities/UIElements/UIElementIdentifierEnum.h"
#include "Utilities/Rendering/Colors.h"

// Standard Library Imports
#include <string>

class FPSCounter : public UIElement
{
  public:
    FPSCounter();
    ~FPSCounter() = default;
    virtual void renderMe(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont) override;
};