// File: Transfer/src/Entities/UIElements/Buttons/Button.cpp

#include "Entities/UIElements/Buttons/Button.hpp"

Button::Button()
{
    boundingRect = SDL_FRect{0.0f, 0.0f, 0.0f, 0.0};
    buttonSelected = false;
    setPosition(boundingRect.x, boundingRect.y);
    hotZoneRect = boundingRect;
}

void Button::clickMe(DynamoEngine::Vector2D positionOfEvent, UIState& uiState)
{
    std::string temp = altText;
    altText = displayText;
    displayText = temp;
    return;
}

void Button::buildGeometry(DynamoEngine::UIGeometryBuilder& builder)
{
    builder.addRect(boundingRect, ColorLibrary::Gray);
    builder.addTextCentered(getDisplayText(), boundingRect, ColorLibrary::White);
}