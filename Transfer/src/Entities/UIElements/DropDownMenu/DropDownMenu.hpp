// File: Transfer/src/Entities/UIElements/DropDownMenu.hpp

#pragma once

#include "Core/UIState.hpp"
#include "DynamoEngine/Math/Vector2.hpp"
#include "Entities/UIElements/Buttons/VisorButton/VisorButton.hpp"
#include <vector>

class DropDownMenu : public UIElement
{
  public:
    DropDownMenu();
    ~DropDownMenu();

    void buildGeometry(DynamoEngine::UIGeometryBuilder& builder) override;

    void clickMe(DynamoEngine::Vector2D positionOfEvent, UIState& uiState) override;

  private:
    // std::vector<VisorButton*> visorButtons;
};