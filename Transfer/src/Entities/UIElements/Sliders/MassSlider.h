// File: Transfer/src/Entities/UIElements/Sliders/MassSlider.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Core/UIState.h"
#include "Entities/MathStructures.h"
#include "Entities/UIElements/Sliders/Slider.h"
#include "Entities/UIElements/UIElement.h"
#include "Entities/UIElements/UIElementTypeEnum.h"
#include "Utilities/Colors.h"
#include "Utilities/EngineConstants.h"
#include "Utilities/GameSystemConstants.h"

// Standard Library Imports
#include <string>

class MassSlider : public Slider {
  public:
    MassSlider();
    ~MassSlider() = default;
    std::string getDisplayText() const override {
        return "Mass: " + std::to_string(sliderValue); // add units?
    }
};
