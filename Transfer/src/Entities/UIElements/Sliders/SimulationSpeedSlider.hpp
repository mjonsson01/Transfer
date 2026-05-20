// File: Transfer/src/Entities/UIElements/Sliders/SimulationSpeedSlider.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Core/UIState.hpp"
#include "Entities/UIElements/Sliders/Slider.hpp"
#include "Entities/UIElements/UIElement.hpp"
#include "Entities/UIElements/UIElementIdentifierEnum.hpp"
#include "Utilities/Constants/EngineConstants.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"
#include "Utilities/Math/Vector2D.hpp"
#include "Utilities/Rendering/Colors.hpp"

// Standard Library Imports
#include <string>

class SimulationSpeedSlider : public Slider
{
  public:
    SimulationSpeedSlider();
    ~SimulationSpeedSlider() = default;
    std::string getDisplayText() const override
    {
        return "Simulation Speed: " + std::to_string(sliderValue); // add units?
    }
};
