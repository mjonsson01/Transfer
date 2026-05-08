// File: Transfer/src/Entities/Physics/GravitationalBody.h

#pragma once

// Custom Imports
#include "Utilities/Math/Vector2D.h"

// Standard Library Imports
#include <SDL3/SDL_stdinc.h>
#include <iostream>
#include <vector>

struct GravitationalBody
{
    // State data members
    Vector2D position = {0.0, 0.0};         // Current position vector (measured from upper left)
    Vector2D previousPosition = {0.0, 0.0}; // Previous frame position vector (used
                                            // for alpha frame interpolation)
    Vector2D velocity = {0.0, 0.0};         // Current frame net velocity vector
    Vector2D netForce = {0.0, 0.0};         // Current frame net force vector
    Vector2D prevForce = {0.0, 0.0};        // Net force vector from previous frame
    double mass = 0.0;                      // Mass of Gravitational Body
    double radius = 0.0;                    // Radius of Gravitational Body (units TBD)
    double temperature = 0.0;               // Kelvin

    // Top-level type flags
    bool isMacro = false;    // Is body with macro-particle and macro-macro gravitational
                             // forces and with macro-particle and macro-macro collisions.
    bool isParticle = false; // Is body with only particle-macro gravitational
                             // interactions and particle-macro collisions.

    // Subtype flags -- only one may be true
    // Need to be implemented and well-defined
    bool isPlanet = false;
    bool isMoon = false;
    bool isGravStar = false;
    bool isDust = false;
    bool isFragment = false;
    bool isGas = false;

    bool isBounce = false; // Body only bounces, no accretion-type behavior. False makes it
                           // follow the standard thresholds. true makes it always bounce.
    // Attribute enablement flags
    bool isStatic = false;      // makes the particle not follow any forces, can only have a
                                // fixed initial velocity. Basically as if it has infinite mass.
    bool isShatterable = false; // allows for shattering into particles
    bool isAccretable = false;  // allows body to be absorbed into another
    bool isCollidable = false;  // makes the body a ghost
    bool isMacroGhost = false;  // makes the body a ghost to macro bodies, but
                                // does not ghost for micro particles.
    bool visible = true;        // particles and macros should be visible by default
                                // but if instantiated specifically for collision
                                // purposes, maybe toggle visibility off
    // Control flag to send to cleanup function
    bool isTransient = false; // for bodies that are bounce artifacts of collisions
    Uint32 lifetime = 0;
    bool isMarkedForDeletion = false;

    int macroIdentifier = -1; // Defaults to minus one on construction unless
                              // explicitly being instantiated with a parent that might
                              // disintegrate. so create macro body should
                              // pass an iterating id based on global state. Only macro
                              // bodies will have an id. Can't just be the number of macro bodies
                              // because then duplicates, so just select +1 on each instantiation.
};

// --------- I/O OPERATOR OVERLOAD --------- //
// Print overloader operator to print full information of a given gravitational
// body.
inline std::ostream& operator<<(std::ostream& os, const GravitationalBody& b)
{
    os << "GravitationalBody{"
       << "\n --- State: ---"
       << "\nposition = " << b.position << "\nprevPosition = " << b.previousPosition << "\nvelocity = " << b.velocity
       << "\nnetForce = " << b.netForce << "\nprevForce = " << b.prevForce << "\nmass = " << b.mass
       << "\nradius = " << b.radius << "\ntemperature = " << b.temperature

       << "\n --- Type Classification: ---"
       << "\nisMacro = " << std::boolalpha << b.isMacro << "\nisParticle = " << std::boolalpha << b.isParticle

       << "\n --- Subtype Classification: ---"
       << "\nisPlanet = " << std::boolalpha << b.isPlanet << "\nisMoon = " << std::boolalpha << b.isMoon
       << "\nisGravStar = " << std::boolalpha << b.isGravStar << "\nisDust = " << std::boolalpha << b.isDust
       << "\nisFragment = " << std::boolalpha << b.isFragment << "\nisGas = " << std::boolalpha << b.isGas

       << "\n --- Attributes: ---"
       << "\nisStatic = " << std::boolalpha << b.isStatic << "\nisShatterable = " << std::boolalpha << b.isShatterable
       << "\nisAccretable = " << std::boolalpha << b.isAccretable << "\nisCollidable = " << std::boolalpha
       << b.isCollidable << "\nisMacroGhost = " << std::boolalpha << b.isMacroGhost << "\nisBounce = " << std::boolalpha
       << b.isBounce

       << "\n --- Status: ---"
       << "\nisTransient = " << std::boolalpha << b.isTransient << "\nlifetime = " << b.lifetime
       << "\nisMarkedForDeletion = " << std::boolalpha << b.isMarkedForDeletion

       << "\n}";
    return os;
}