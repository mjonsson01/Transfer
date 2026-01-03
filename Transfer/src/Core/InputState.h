// File: Transfer/src/Core/InputState.h

#pragma once
#include "Entities/PhysicsStructures.h"


struct InputState
{   
    // Persistent Flags
    Vector2D mouseCurrPosition = {0.0, 0.0};
    Vector2D mouseDragStartPosition = {0.0, 0.0};
    Vector2D mouseDragEndPosition = {0.0, 0.0};
    double selectedMass = 0.0;
    double selectedRadius = 0.0;
    bool isMouseDragging = false;
    bool isHoldingLeftMouseButton = false;
    bool isHoldingRightMouseButton = false;

    // Transient Flags reset at end of processing an event
    bool dirty = false;

    bool isCreatingCluster = false;
    bool isCreatingPlanet = false;
    bool isCreatingDust = false;
    bool isCreatingStatic = false;
    double spawnAccumulator = 0.0;
    
    // Managed by player input
    bool clearAll = false;
    bool isPaused = false;


    InputState() = default;
    ~InputState() = default;
    InputState& resetTransientFlags()
    {
        isCreatingCluster = false;
        isCreatingPlanet = false;
        isCreatingDust = false;
        isCreatingStatic = false;
        dirty = false;
        spawnAccumulator = 0.0;
        return *this;
    }
};

