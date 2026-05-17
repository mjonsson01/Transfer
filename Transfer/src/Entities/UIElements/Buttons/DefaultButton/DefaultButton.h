// File: Transfer/src/Entities/UIElements/Buttons/DefaultButton/DefaultButton.h

#pragma once

#include "Entities/UIElements/Buttons/Button.h"

class DefaultButton : public Button
{
  public:
    DefaultButton();
    ~DefaultButton() = default;
    void clickMe(Vector2D positionOfEvent, UIState& UIState) override;
};