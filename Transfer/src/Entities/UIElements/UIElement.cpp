// File: Transfer/src/Entities/UIElements/UIElement.cpp

#include "UIElement.h"

// Base UIElement implementation
UIElement::UIElement() {}

UIElement::~UIElement() {}

UIElementIdentifier UIElement::checkAndReturnIfHit(const Vector2D& positionToCheck)
{
    SDL_FPoint point = {static_cast<float>(positionToCheck.xVal), static_cast<float>(positionToCheck.yVal)};
    UIElementIdentifier type_hit = UIElementIdentifier::NONE;
    if (SDL_PointInRectFloat(&point, &hotZoneRect))
    {
        type_hit = UIElementIdentifierIdentifier;
    }
    return type_hit;
}