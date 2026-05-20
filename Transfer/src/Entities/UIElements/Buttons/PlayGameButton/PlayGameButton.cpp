// File: Transfer/src/Entities/UIElements/Buttons/PlayGameButton/PlayGameButton.cpp

#include "Entities/UIElements/Buttons/PlayGameButton/PlayGameButton.hpp"

PlayGameButton::PlayGameButton() : Button()
{
    float width = 300.0f;
    float height = 200.0f;
    boundingRect = SDL_FRect{SCREEN_WIDTH / 2 - width / 2, SCREEN_HEIGHT / 2 - height / 2, width, height};
    hotZoneRect = boundingRect;
    setPosition(boundingRect.x, boundingRect.y);
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