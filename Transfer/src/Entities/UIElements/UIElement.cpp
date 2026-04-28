// File: Transfer/src/Entities/UIElements/UIElement.cpp

#include "UIElement.h"

// Base UIElement implementation
UIElement::UIElement() {}

UIElement::~UIElement() {}

UIElementType UIElement::checkAndReturnIfHit(const Vector2D& positionToCheck)
{
    // std::cout<<"check and return if hit called"<<std::endl;
    SDL_FPoint point = {static_cast<float>(positionToCheck.xVal), static_cast<float>(positionToCheck.yVal)};
    UIElementType type_hit = UIElementType::NONE;
    if (SDL_PointInRectFloat(&point, &hotZoneRect))
    {
        type_hit = UIElementTypeIdentifier;
    }
    return type_hit;
}