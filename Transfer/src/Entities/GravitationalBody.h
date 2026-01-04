// File: Transfer/src/Entities/GravitationalBody.h

#pragma once

// Custom Imports
#include "Entities/PhysicsStructures.h"

// Standard Library Imports
#include <cstdint>
#include <iostream>
#include <vector>

struct GravitationalBody
{
    // State data members
    Vector2D position = {0.0, 0.0};     // Current position vector (measured from upper left)
    Vector2D previousPosition = {0.0, 0.0}; // Previous frame position vector (used for alpha frame interpolation)
    Vector2D velocity = {0.0, 0.0};     // Current frame net velocity vector
    Vector2D netForce = {0.0, 0.0};     // Current frame net force vector
    Vector2D prevForce = {0.0, 0.0};    // Net force vector from previous frame
    double mass = 0.0;                  // Mass of Gravitational Body
    double radius = 0.0;                // Radius of Gravitational Body (units TBD)
    double temperature = 0.0;           // Kelvin
    
    // Top-level type flags
    bool isMacro = false;    // Is body with macro-particle and macro-macro gravitational forces and with macro-particle and macro-macro collisions.
    bool isParticle = false; // Is body with only particle-macro gravitational interactions and particle-macro collisions.
    
    // Subtype flags -- only one may be true
    // Need to be implemented and well-defined
    bool isPlanet = false;
    bool isMoon = false;
    bool isGravStar = false;
    bool isDust = false;
    bool isFragment = false;
    bool isGas = false;
    
    // Attribute enablement flags
    bool isStatic = false;      // fixes the position of the body (doesn't get integrated forward)
    bool isShatterable = false; // allows for shattering into particles
    bool isAccretable = false;  // allows body to be absorbed into another
    bool isCollidable = false;  // makes the particle a ghost
    
    // Control flag to send to cleanup function
    bool isMarkedForDeletion = false;
};


// --------- I/O OPERATOR OVERLOAD --------- //
// Print overloader operator to print full information of a given gravitational body.
inline std::ostream& operator<<(std::ostream& os, const GravitationalBody& b)
{
    os << "GravitationalBody{"
       << "\n --- State: ---"
       << "\nposition = " << b.position
       << "\nprevPosition = " << b.previousPosition
       << "\nvelocity = " << b.velocity
       << "\nnetForce = " << b.netForce
       << "\nprevForce = " << b.prevForce
       << "\nmass = " << b.mass
       << "\nradius = " << b.radius
       << "\n --- Type Classification: ---"
       << "\nisMacro = " << std::boolalpha << b.isMacro
       << "\nisParticle = " << std::boolalpha << b.isParticle
       << "\n --- Subtype Classification: ---"
       << "\nisPlanet = " << std::boolalpha << b.isPlanet
       << "\nisMoon = " << std::boolalpha << b.isMoon
       << "\nisGravStar = " << std::boolalpha << b.isGravStar
       << "\nisDust = " << std::boolalpha << b.isDust
       << "\nisFragment = " << std::boolalpha << b.isFragment
       << "\nisGas = " << std::boolalpha << b.isGas
       << "\n --- Attributes: ---"
       << "\nisStatic = " << std::boolalpha << b.isStatic
       << "\nisShatterable = " << std::boolalpha << b.isShatterable
       << "\nisAccretable = " << std::boolalpha << b.isAccretable
       << "\nisCollidable = " << std::boolalpha << b.isCollidable
       << "\n --- Status: ---"
       << "\nisMarkedForDeletion = " <<std::boolalpha << b.isMarkedForDeletion
       << "\n}";
    return os;
}
