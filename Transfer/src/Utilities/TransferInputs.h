// File: Transfer/src/Utilities/TransferInputs.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>

// Custom Imports
#include "Entities/MathStructures.h"

// Standard Library Imports
#include <iostream>

struct TransferInputs
{
    TransferInputs() = default;
    TransferInputs(const TransferInputs&) = default;
    TransferInputs& operator=(const TransferInputs&) = default;
    // Standard Keyboard Controls
    bool wPressed = false;
    bool aPressed = false;
    bool sPressed = false;
    bool dPressed = false;
    bool spacePressed = false;
    bool clearParticlesPressed = false;

    // Instantiation controls
    bool leftMousePressed = false;
    bool rightMousePressed = false;
    bool middleMousePressed = false;
    bool leftMouseJustReleased = false;
    bool leftMouseJustPressed = false;

    // Menu and Media controls
    bool escPressed = false;
    bool zeroPressed = false;    // Music play/pause
    bool escJustPressed = false; // Pause menu

    // Modifier Controls
    // bool ctrlPressed = false;
    bool shiftPressed = false;
    bool altPressed = false;

    Vector2D mouseCurrPosition = {0, 0};
    Vector2D mouseDragStartPosition = {0, 0};
    bool isDragging = false;

    // Helper utilities
    void resetAllMousePressedVars()
    {
        this->leftMousePressed = false;
        this->rightMousePressed = false;
        this->middleMousePressed = false;
    }
    void resetJustPressed()
    {
        leftMouseJustPressed = false;
        leftMouseJustReleased = false;
        escJustPressed = false;
    }
    void resetAllKeyPressedVars()
    {
        this->wPressed = false;
        this->aPressed = false;
        this->sPressed = false;
        this->dPressed = false;
        this->spacePressed = false;
        this->zeroPressed = false;
        this->shiftPressed = false;
        this->altPressed = false;
    }
};

std::ostream& operator<<(std::ostream& os, const TransferInputs& input);