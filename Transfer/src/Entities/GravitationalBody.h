// File: Transfer/src/Entities/GravitationalBody.h

#pragma once

// Custom Imports
#include "Entities/PhysicsStructures.h"

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

    BoundingBox getBoundingBox() const { return boundingBox;}
    void setBoundingBox(BoundingBox updatedBoundingBox) {boundingBox = updatedBoundingBox;}

    private:
        float radius = 0.0; // in pixels
        double mass = 0.0; // in kilograms
        Vector2D bodyNetVelocity = { 0.0, 0.0 }; // in pixels per second?
        Vector2D bodyNetForce = { 0.0, 0.0 }; // in Newtons (may have to shift to different unit later)
        Vector2D bodyPosition = { 0.0, 0.0 }; // position measured in pixel xy coordinates measured from top-left corner.
        Vector2D prevPosition = { 0.0, 0.0 }; // previous pixel position
        BoundingBox boundingBox = { 0.0, 0.0, 0.0, 0.0};
};