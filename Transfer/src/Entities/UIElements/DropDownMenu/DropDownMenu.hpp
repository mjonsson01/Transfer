// File: Transfer/src/Entities/UIElements/DropDownMenu.hpp

#pragma once 

#include "Entities/UIElements/Buttons/VisorButton/VisorButton.hpp"

class DropDownMenu : public UIElement
{
    DropDownMenu();
    ~DropDownMenu();

    void buildGeometry(std::vector<UIElementVertex>& vertexBuffer, uint32_t zIndex,
                       const FontAtlasUtility& fontAtlas) override; 

    private: 
        std::vector<VisorButton&> visorButtons;
};