// File: Transfer/src/Entities/UIElements/Buttons/DefaultButton/DefaultButton.cpp

#include "Entities/UIElements/Buttons/DefaultButton/DefaultButton.h"

DefaultButton::DefaultButton() : Button()
{
    boundingRect = SDL_FRect{100.0f, 100.0f, 300.0f, 200.0f};
    hotZoneRect = boundingRect;
    setPosition(boundingRect.x, boundingRect.y);
    displayText = "Default Button";
    altText = "Clicked";
    bool buttonSelected = false;
    UIElementID = UIElementIdentifier::DEFAULT_BUTTON_INDEX;
}

void DefaultButton::clickMe(Vector2D positionOfEvent, UIState& UIState)
{
    // UIState.setCurrentScene(GAME_SCENE);
    std::string temp = altText;
    altText = displayText;
    displayText = temp;
    return;
}