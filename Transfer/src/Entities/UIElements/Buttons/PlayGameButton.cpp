// File: Transfer/src/Entities/UIElements/Buttons/PlayGameButton.cpp

#include "Entities/UIElements/Buttons/PlayGameButton.h"

PlayGameButton::PlayGameButton() : Button()
{
    boundingRect = SDL_FRect{100.0f, 100.0f, 300.0f, 200.0f};
    hotZoneRect = boundingRect;
    setPosition(boundingRect.x, boundingRect.y);
    displayText = "Play Game";
    altText = "Clicked";
    bool buttonSelected = false;
    UIElementTypeIdentifier = UIElementType::PLAY_GAME_BUTTON_INDEX;
}