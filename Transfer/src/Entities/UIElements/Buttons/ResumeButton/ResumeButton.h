// File: Transfer/src/Entities/UIElements/Buttons/ResumeButton/ResumeButton.h

#pragma once

#include "Entities/UIElements/Buttons/Button.h"

class ResumeButton : public Button
{
  public:
    ResumeButton();
    ~ResumeButton() = default;
    void clickMe(Vector2D positionOfEvent, UIState& UIState) override;
};