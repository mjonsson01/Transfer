// File: Transfer/src/Entities/PhysicsStructures.h

#pragma once

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
static inline Vector2D lerp(const Vector2D& a, const Vector2D& b, double t)
{
    return a + (b - a) * t;
}

// **DECLARATION ONLY:** The compiler just needs to know this function exists.
std::ostream& operator<<(std::ostream& os, const Vector2D& Vec);

