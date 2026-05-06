// File: Transfer/src/Core/InputState.h

#pragma once

#include "Utilities/Math/Vector2D.h"

struct InputState
{
    // Persistent flags and details
    Vector2D mouseCurrPosition = {0.0, 0.0}; // set by from event input router
    Vector2D mouseDragStartPosition = {0.0, 0.0};

    // Physics locations if instantiateDirty gets thrown. If the event is not
    // consumed by a UI event, then instantiate dirty will be set by the
    // UIState.

    bool isDragging = false;                  // set in Input System
    bool isClickingRightMouseButton = false;  // set by transfer inputs
    bool isClickingMiddleMouseButton = false; // set by transfer inputs
    bool isClickingLeftMouseButton = false;   // set by tranfer inputs
    bool leftMouseButtonJustPressed = false;
    bool leftMouseButtonJustReleased = false;
    bool isPressingShift = false;

    double selectedMass = 0.0;   // spit back up from UI element to UISystem to UIState
    double selectedRadius = 0.0; // spit back up from UI element to UISystem to UIState
    double selectedSimSpeedScale = 1.0; // spit back up from UI element to UISystem to UIState

    bool UIInputConsumed = false; // Does not get reset as a transient flag

    // Creation type flags
    bool isCreatingMacro = false;
    bool isCreatingParticle = false;
    bool isCreatingParticleCluster = false;

    // Creation subtype flags -- only one may be true
    // Need to be implemented and well-defined
    bool isCreatingPlanet = false;
    bool isCreatingMoon = false;
    bool isCreatingGravStar = false;
    bool isCreatingDust = false;
    bool isCreatingFragment = false;
    bool isCreatingGas = false;

    // Attribute enablement flags on creation
    bool isCreatingStatic = false;
    bool isCreatingShatterable = false;
    bool isCreatingAccretable = false;
    bool isCreatingCollidable = false;
    bool isCreatingBounce = false;
    bool isCreatingMacroGhost = false;
    bool isCreatingWithInitialVelocity = false;

    // Managed by player input
    bool clearAll = false;        // Toggled when screen wipe is requested
    bool isPhysicsPaused = false; // Toggled when rendering continues but Physics
                                  // System integration is completely paused
    bool isPaused = false;
    bool isPreviewingMacro = false;
    bool isPreviewingWithInitialVelocity = false;

    InputState& resetTransientFlags()
    {
        // Creation type flags
        isCreatingMacro = false;
        isCreatingParticle = false;
        isCreatingParticleCluster = false;

        // Creation subtype flags -- only one may be true
        // Need to be implemented and well-defined
        isCreatingPlanet = false;
        isCreatingMoon = false;
        isCreatingGravStar = false;
        isCreatingDust = false;
        isCreatingFragment = false;
        isCreatingGas = false;

        // Attribute enablement flags on creation
        isCreatingStatic = false;
        isCreatingShatterable = false;
        isCreatingAccretable = false;
        isCreatingCollidable = false;
        isCreatingMacroGhost = false;
        isCreatingWithInitialVelocity = false;

        leftMouseButtonJustPressed = false;
        leftMouseButtonJustReleased = false;
        return *this;
    }

    // Toggle the pausing of the integrator while still rendering background
    // elements
    InputState& togglePhysicsPause()
    {
        isPhysicsPaused = !isPhysicsPaused;
        return *this;
    }

    // Wipe all the bodies from the screen and internal Game State
    InputState& clearAllBodies()
    {
        clearAll = true;
        return *this;
    }
};
