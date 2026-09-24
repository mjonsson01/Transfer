// File: Transfer/src/Entities/UIElements/UIElement.cpp

#include "UIElement.hpp"

// Base UIElement implementation
UIElement::UIElement() {}

UIElement::~UIElement() {}

UIElementIdentifier UIElement::checkAndReturnIfHit(const DynamoEngine::Vector2D& positionToCheck)
{
    SDL_FPoint point = {static_cast<float>(positionToCheck.x_val), static_cast<float>(positionToCheck.y_val)};
    UIElementIdentifier type_hit = UIElementIdentifier::NONE;
    if (SDL_PointInRectFloat(&point, &hotZoneRect))
    {
        type_hit = UIElementID;
    }
    return type_hit;
}