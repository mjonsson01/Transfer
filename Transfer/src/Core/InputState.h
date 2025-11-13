// File: Transfer/src/Core/InputState.h

#pragma once

class InputState
{
    public:
        // Constructor and Destructor
        InputState();
        ~InputState();
    public:
        bool isMouseDragging = false;
        float mouseDragStartX = 0.0f;
        float mouseDragStartY = 0.0f;
        float mouseDragEndX = 0.0f;
        float mouseDragEndY = 0.0f;
        // Body instantiation control vars
        // float selectedRadius = 0.0;
        // double selectedMass = 0.0;
        // bool bodySelectionValidity = false;
};