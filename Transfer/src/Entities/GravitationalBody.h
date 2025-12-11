// File: Transfer/src/Entities/GravitationalBody.h

#pragma once

// Custom Imports
#include "Entities/PhysicsStructures.h"

// Standard Library Imports
#include <cstdint>
#include <iostream>

// To be implemented in the future for non-circular gravitational bodies

// enum class ShapeType { Circle, Polygon };

// struct PolygonShape {
//     std::vector<Vector2D> vertices;
// };

// class GravitationalBody
// {
// public:
//     ShapeType shape = ShapeType::Circle;
//     PolygonShape polygon;   // only used when shape == Polygon

//     // everything else stays the same...
// };

class GravitationalBody
{
    public:
        // Constructor and Destructor
        GravitationalBody();
        ~GravitationalBody();

    public:
    // Getters and Setters

    float getRadius() const { return radius; }
    void setRadius(float r) { radius = r; }

    float getRadiusSquared() const {return radius*radius;}
    
    double getMass() const { return mass; }
    void setMass(double m) { mass = m; }

    Vector2D getPosition() const { return bodyPosition; }
    void setPosition(const Vector2D& position) { bodyPosition = position; }

    Vector2D getPrevPosition() const {return prevPosition; }
    void setPrevPosition(const Vector2D& position) { prevPosition = position; }

    
    Vector2D getNetVelocity() const { return bodyNetVelocity; }
    void setNetVelocity(const Vector2D& velocity) { bodyNetVelocity = velocity;}

    Vector2D getNetForce() const { return bodyNetForce; }
    void setNetForce(const Vector2D& force) { bodyNetForce = force;}

    void addNetForce(const Vector2D& forceComp) {bodyNetForce = bodyNetForce + forceComp;}

    BoundingBox getBoundingBox() const { return boundingBox;}
    void setBoundingBox(BoundingBox updatedBoundingBox) {boundingBox = updatedBoundingBox;}

    // bool getIsColliding() const { return isColliding; }
    // void setIsColliding(bool colliding) { isColliding = colliding; }

    bool getCollisionEnabled() const { return collisionEnabled; }
    void setCollisionEnabled(bool enabled) { collisionEnabled = enabled; }

    bool getMarkedForDeletion() const { return markedForDeletion; }
    void setMarkedForDeletion(bool marked) { markedForDeletion = marked; }

    bool getIsFragment() const {return isFragment;}
    void setIsFragment(bool isFrag) {isFragment = isFrag;}

    friend std::ostream& operator<<(std::ostream& os, const GravitationalBody& b);
    // bool getIsGhost() const {return isGhost;}
    // int getGhostFramesRemaining() {return ghostFramesRemaining;}

    // void setIsGhost(bool is_ghost) {isGhost = is_ghost;}
    // void setGhostFramesRemaining(int num) {ghostFramesRemaining = num;}
    // void decrementGhostFrames() {ghostFramesRemaining -= 1;}
    // void setGhost(int frames) {
    //     isGhost = true;
    //     ghostFramesRemaining = frames;
    // }

    // bool getIsGhost() const { return isGhost; }

    // void updateGhostState() {
    //     if (isGhost) {
    //         if (ghostFramesRemaining > 0) {
    //             ghostFramesRemaining--;
    //         } else {
    //             isGhost = false; // re-enable after cooldown
    //         }
    //     }
    // }

    private:
        float radius = 0.0; // in pixels
        double mass = 0.0; // in kilograms
        Vector2D bodyNetVelocity = { 0.0, 0.0 }; // in pixels per second?
        Vector2D bodyNetForce = { 0.0, 0.0 }; // in Newtons (may have to shift to different unit later)
        Vector2D bodyPosition = { 0.0, 0.0 }; // position measured in pixel xy coordinates measured from top-left corner.
        Vector2D prevPosition = { 0.0, 0.0 }; // previous pixel position
        BoundingBox boundingBox = { 0.0, 0.0, 0.0, 0.0};
        // bool isStatic = false; // whether the body is static (immovable) or dynamic

        // uint32_t collisionCooldownMs = 0;
        bool collisionEnabled = false;

        bool markedForDeletion = false;

        bool isFragment = false;

        // bool isGhost = false;
        // int ghostFramesRemaining = 0;

};
