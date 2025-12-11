// File: Transfer/src/Entities/GravitationalBody.h

#pragma once

// Custom Imports
#include "Entities/PhysicsStructures.h"

// Standard Library Imports
#include <cstdint>
#include <iostream>
#include <vector>

// To be implemented in the future for non-circular gravitational bodies
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


struct Particle {
    Vector2D position {0.0, 0.0};
    Vector2D prevPosition = {0.0, 0.0};
    Vector2D velocity = {0.0, 0.0};
    Vector2D netForce = {0.0, 0.0};
    Vector2D prevForce = {0.0, 0.0};
    double mass = 0.0;
    double radius = 0.0;
    int clusterID; // ID of the cluster this particle belongs to. -1 if belongs to none.
    bool collisionEnabled = true;
    bool markedForDeletion = false;
};

struct GravitationalCluster {
    // std::vector<Particle> clusterParticles;
    std::vector<int> clusterParticleIndices; // Indices of particles in the GameState's particle vector that belong to this cluster
    double initialRadius = 0.0;
    double totalMass = 0.0;
    Vector2D centerOfMassPos {0.0, 0.0};
    Vector2D bulkVelocity {0.0, 0.0}; // Velocity of bulk
    bool isDebrisField = false;
};


