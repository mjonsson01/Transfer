// File: Transfer/src/Entities/UIElements/Buttons/PlayGameButton/PlayGameButton.h

#pragma once

#include "Entities/UIElements/Buttons/Button.h"

class PlayGameButton : public Button
{
  public:
    PlayGameButton();
    ~PlayGameButton() = default;
    void clickMe(Vector2D positionOfEvent, UIState& UIState) override;
};