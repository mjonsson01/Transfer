// File: Transfer/src/Entities/UIElements/Overlay/FPSCounter.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Core/UIState.hpp"
#include "Entities/UIElements/UIElement.hpp"
#include "Entities/UIElements/UIElementIdentifierEnum.hpp"
#include "Utilities/Rendering/Colors.hpp"
#include "Utilities/Rendering/FontAtlasUtility.hpp"

// Standard Library Imports
#include <string>

class FPSCounter : public UIElement
{
  public:
    FPSCounter();
    ~FPSCounter() = default;
    virtual void buildGeometry(std::vector<UIElementVertex>& vertexBuffer, uint32_t zIndex,
                               const FontAtlasUtility& fontAtlas) override;
    void updateMe(UIState& UIState) override;
    void updateLayout(float windowWidth, float windowHeight) override;
    std::string getDisplayText() const
    {
        return std::to_string(fps);
    } // Maybe make virtual later but not seeing a good reason rn

  private:
    int fps; // Simply local from the UIState call for updateMe;
};