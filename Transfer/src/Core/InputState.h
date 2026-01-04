// File: Transfer/src/Core/InputState.h

#pragma once
#include "Entities/PhysicsStructures.h"

struct InputState
{   
    // Persistent Flags
    Vector2D mouseCurrPosition = {0.0, 0.0};
    Vector2D mouseDragStartPosition = {0.0, 0.0};
    Vector2D mouseDragEndPosition = {0.0, 0.0};
    bool isMouseDragging = false;
    bool isHoldingRightMouseButton = false;
    bool isHoldingLeftMouseButton = false;
    double selectedMass = 0.0;
    double selectedRadius = 0.0;

    
    // Transient Flags reset at end of processing an event
    bool dirty = false;             // Flag for body instantiation/rendering required

    // Creation type flags
    bool isCreatingMacro = false;    
    bool isCreatingParticle = false;
    // bool isCreatingParticleCluster = false;

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

    double spawnAccumulator = 0.0;  // To track time elapsed between last spawn.

    
    // Managed by player input
    bool clearAll = false; // Toggled when screen wipe is requested
    bool isPhysicsPaused = false; // Toggled when rendering continues but Physics System integration is completely paused
    // bool isPaused = false;

    InputState& resetTransientFlags()
    {
        dirty = false;             // Flag for body instantiation/rendering required

        // Creation type flags
        isCreatingMacro = false;    
        isCreatingParticle = false;
        // isCreatingParticleCluster = false;
        
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

        double spawnAccumulator = 0.0;  // To track time elapsed between last spawn.
        return *this;
    }
    
    // Toggle the pausing of the integrator while still rendering background elements
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

