// File: Transfer/src/Entities/UIElements/Buttons/PlayGameButton/PlayGameButton.cpp

#include "Entities/UIElements/Buttons/PlayGameButton/PlayGameButton.hpp"

PlayGameButton::PlayGameButton() : Button()
{
    updateLayout(SCREEN_HEIGHT, SCREEN_WIDTH);
    setVisibility(true);
    displayText = "Play Game";
    altText = "";
    bool buttonSelected = false;
    UIElementID = UIElementIdentifier::PLAY_GAME_BUTTON_INDEX;
}

void PlayGameButton::clickMe(Vector2D positionOfEvent, UIState& UIState)
{
    UIState.setCurrentScene(SceneIdentifier::GAME_SCENE);
    UIState.QueueSoundEffect("ButtonClick");
    return;
}

void PlayGameButton::updateLayout(float windowWidth, float windowHeight)
{
    float width = 300.0f;
    float height = 200.0f;
    boundingRect = SDL_FRect{windowWidth / 2 - width / 2, windowHeight / 2 - height / 2, width, height};
    setPosition(boundingRect.x, boundingRect.y);
    hotZoneRect = boundingRect;
}