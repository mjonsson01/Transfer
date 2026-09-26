// File: Transfer/src/Entities/UIElements/Buttons/PlayGameButton/PlayGameButton.cpp

#include "Entities/UIElements/Buttons/PlayGameButton/PlayGameButton.hpp"

PlayGameButton::PlayGameButton() : Button()
{
    updateLayout(SCREEN_WIDTH, SCREEN_HEIGHT);
    setVisibility(true);
    displayText = "Play Game";
    altText = "";
    bool buttonSelected = false;
    UIElementID = UIElementIdentifier::PLAY_GAME_BUTTON_INDEX;
}

void PlayGameButton::clickMe(DynamoEngine::Vector2D positionOfEvent, UIState& uiState)
{
    uiState.setCurrentScene(SceneIdentifier::GAME_SCENE);
    uiState.QueueSoundEffect("ButtonClick");
    return;
}

void PlayGameButton::updateLayout(float window_width, float window_height)
{
    float width = 300.0f;
    float height = 200.0f;
    boundingRect = SDL_FRect{window_width / 2 - width / 2, window_height / 2 - height / 2, width, height};
    setPosition(boundingRect.x, boundingRect.y);
    hotZoneRect = boundingRect;
}