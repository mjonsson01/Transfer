// File: Transfer/src/Core/InputState.h

#pragma once

#include "Entities/MathStructures.h"

struct InputState
{
    // Persistent flags and details
    Vector2D mouseCurrPosition = {0.0, 0.0}; // set by from event input router
    Vector2D mouseDragStartPosition = {0.0, 0.0};

    // Physics locations if instantiateDirty gets thrown. If the event is not
    // consumed by a UI event, then instantiate dirty will be set by the
    // UIState.
    Vector2D instantiatePosition = {0.0, 0.0}; // derived from event input router, set by UISystem if determined
                                               // that click event did not interact with a UI element.
    Vector2D instantiateDragStartPosition = {0.0, 0.0}; // only possible in game or level editor, derived from transfer
                                                        // inputs / editor inputs and set in the Input System

    bool isDragging = false;                 // set in Input System
    bool isHoldingRightMouseButton = false;  // set by transfer inputs
    bool isHoldingMiddleMouseButton = false; // set by transfer inputs
    bool isHoldingLeftMouseButton = false;   // set by tranfer inputs

    double selectedMass = 0.0;   // spit back up from UI element to UISystem to UIState
    double selectedRadius = 0.0; // spit back up from UI element to UISystem to UIState

    bool UIInputConsumed = false; // Does not get reset as a transient flag

    // Transient Flags reset at end of processing an event
    bool instantiateDirty = false; // Flag for body instantiation/rendering required. If dirty is
                                   // set to true, instantiate at the instantiate position then flip
                                   // back to false.

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

    // Managed by player input
    bool clearAll = false;        // Toggled when screen wipe is requested
    bool isPhysicsPaused = false; // Toggled when rendering continues but Physics
                                  // System integration is completely paused
    // bool isPaused = false;

    InputState& resetTransientFlags()
    {
        instantiateDirty = false; // Flag for body instantiation/rendering required

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
