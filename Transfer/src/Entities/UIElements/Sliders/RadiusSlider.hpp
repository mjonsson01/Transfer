// File: Transfer/src/Entities/UIElements/Sliders/RadiusSlider.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Core/UIState.hpp"
#include "DynamoEngine/Math/Vector2.hpp"
#include "Entities/UIElements/Sliders/Slider.hpp"
#include "Entities/UIElements/UIElement.hpp"
#include "Entities/UIElements/UIElementIdentifierEnum.hpp"
#include "Utilities/Constants/EngineConstants.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"
#include "Utilities/Rendering/Colors.hpp"

// Standard Library Imports
#include <string>

class RadiusSlider : public Slider
{
  public:
    RadiusSlider();
    ~RadiusSlider() = default;
    std::string getDisplayText() const override
    {
        return "Radius: " + std::to_string(sliderValue); // add units?
    }
    void updateLayout(float window_width, float window_height) override;
};
