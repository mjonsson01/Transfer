// File: Transfer/src/Entities/PhysicsStructures.h

#pragma once

#include "Utilities/GameSystemConstants.h"
#include "Utilities/EngineConstants.h"

// Standard Library Imports
#include <cmath>
#include <iostream>

struct Vector2D
{
	double x_val;
	double y_val;

	// Constructor
	Vector2D(): x_val(0.0), y_val(0.0){};
	Vector2D(double x, double y): x_val(x), y_val(y){};

	 // --- Arithmetic operators ---
    Vector2D operator+(const Vector2D& other) const {
        return {x_val + other.x_val, y_val + other.y_val};
    }

    Vector2D operator-(const Vector2D& other) const {
        return {x_val - other.x_val, y_val - other.y_val};
    }

    Vector2D operator*(double scalar) const {
        return {x_val * scalar, y_val * scalar};
    }

    Vector2D operator/(double scalar) const {
        return {x_val / scalar, y_val / scalar};
    }
	 // --- Compound assignment (optional but handy) ---
    Vector2D& operator+=(const Vector2D& other) {
        x_val += other.x_val;
        y_val += other.y_val;
        return *this;
    }

    Vector2D& operator-=(const Vector2D& other) {
        x_val -= other.x_val;
        y_val -= other.y_val;
        return *this;
    }

    Vector2D& operator*=(double scalar) {
        x_val *= scalar;
        y_val *= scalar;
        return *this;
    }

    Vector2D& operator/=(double scalar) {
        x_val /= scalar;
        y_val /= scalar;
        return *this;
    }
	
    // custom utilities
    double magnitude() const {
        return sqrt(x_val * x_val + y_val * y_val);
    }
	double square_magnitude() const {
		return (x_val * x_val + y_val * y_val);
	}

    double dot(const Vector2D& other){
        return x_val*other.x_val + y_val*other.y_val;
    }
	Vector2D& normalize() {
    	double mag = magnitude();
    	if (mag != 0.0) {
        	x_val /= mag;
        	y_val /= mag;
   		}
		return *this;
	}
};
// **DECLARATION ONLY:** The compiler just needs to know this function exists.
std::ostream& operator<<(std::ostream& os, const Vector2D& Vec);


struct BoundingBox
{
	double min_X;
	double max_X;
	double min_Y;
	double max_Y;

	// checks for overlaps with another bounding box
	bool overlaps(const BoundingBox& other) const {
		return !(max_X < other.min_X ||
        min_X > other.max_X ||
        max_Y < other.min_Y ||
        min_Y > other.max_Y);
	}
};


struct XPBDSoftParticle
{
    Vector2D position = {0.0, 0.0};
    Vector2D prevPosition = {0.0, 0.0};
    Vector2D velocity = {0.0, 0.0};
    double radius = PARTICLE_RADIUS;
    double mass = 0;
    double inv_mass = 0;
};
struct XPBDConstraints
{
    int idx1, idx2; // Particle Indexes
    double restLength;
    double compliance; //XPBD softness
    double lambda; // Lagrange multiplier from prev frame.
    bool isBroken;
};

struct GridCell
{
    std::vector<int> particleIndices;
};


struct PhysicsData {
    double cellSize = PARTICLE_RADIUS;
    int gridWidth = int(SCREEN_WIDTH/cellSize) + 1;
    int gridHeight = int(SCREEN_HEIGHT/cellSize) + 1;
    
    std::vector<XPBDSoftParticle> particles;
    std::vector<XPBDConstraints> constraints;
    std::vector<GridCell> grid;

};