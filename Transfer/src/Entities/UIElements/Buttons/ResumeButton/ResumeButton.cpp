// File: Transfer/src/Entities/UIElements/Buttons/ResumeButton/ResumeButton.cpp

#include "Entities/UIElements/Buttons/ResumeButton/ResumeButton.h"

ResumeButton::ResumeButton() : Button()
{
    float width = 200.0f;
    float height = 100.0f;
    boundingRect = SDL_FRect{SCREEN_WIDTH / 2 - width / 2, SCREEN_HEIGHT / 2 - height / 2, width, height};
    hotZoneRect = boundingRect;
    setPosition(boundingRect.x, boundingRect.y);
    displayText = "Resume";
    altText = "Clicked";
    bool buttonSelected = false;
    UIElementID = UIElementIdentifier::RESUME_BUTTON_INDEX;
}

void ResumeButton::clickMe(Vector2D positionOfEvent, UIState& UIState)
{
    // UIState.setCurrentScene(SceneIdentifier::GAME_SCENE);
    std::string temp = displayText;
    displayText = altText;
    altText = temp;
    UIState.QueueSoundEffect("ButtonClick");
    return;
}