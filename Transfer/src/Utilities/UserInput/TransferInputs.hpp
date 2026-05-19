// File: Transfer/src/Utilities/UserInput/TransferInputs.h

#pragma once

// Custom Imports
#include "Utilities/Math/Vector2D.hpp"

// Standard Library Imports

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
    void resetAllInputsForSceneChange()
    {
        this->wPressed = false;
        this->aPressed = false;
        this->sPressed = false;
        this->dPressed = false;
        this->spacePressed = false;
        this->clearParticlesPressed = false;
        this->leftMousePressed = false;
        this->rightMousePressed = false;
        this->middleMousePressed = false;
        this->leftMouseJustReleased = false;
        this->leftMouseJustPressed = false;
        this->escPressed = false;
        this->zeroPressed = false;
        this->escJustPressed = false;
        this->shiftPressed = false;
        this->altPressed = false;
        this->isDragging = false;
    }
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