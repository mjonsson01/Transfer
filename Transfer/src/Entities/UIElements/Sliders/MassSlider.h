// File: Transfer/src/Entities/UIElements/Sliders/MassSlider.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
// #include "Core/UIState.h"
#include "Entities/UIElements/Sliders/Slider.h"
#include "Entities/UIElements/UIElement.h"
#include "Entities/UIElements/UIElementIdentifierEnum.h"
#include "Utilities/Constants/EngineConstants.h"
#include "Utilities/Constants/GameSystemConstants.h"
#include "Utilities/Math/Vector2D.h"
#include "Utilities/Rendering/Colors.h"

// Standard Library Imports
#include <algorithm>
#include <string>

class MassSlider : public Slider
{
  public:
    MassSlider();
    ~MassSlider() = default;
    std::string getDisplayText() const override
    {
        return "Mass: " + std::to_string(sliderValue); // add units?
    }
    void slideMe(Vector2D positionOfEvent, double& returnedElementValue) override;
};
