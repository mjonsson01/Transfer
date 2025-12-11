// File: Transfer/src/Core/InputState.h

#pragma once
#include "Entities/PhysicsStructures.h"


struct InputState
{   
    bool dirty = false;
    bool isMouseDragging = false;
    bool isCreatingCluster = false;
    bool isCreatingPlanet = false;
    bool isCreatingDust = false;
    bool isCreatingStatic = false;
    Vector2D mouseCurrPosition = {0.0, 0.0};
    Vector2D mouseDragStartPosition = {0.0, 0.0};
    Vector2D mouseDragEndPosition = {0.0, 0.0};
    double selectedMass = 0.0;
    double selectedRadius = 0.0;
};