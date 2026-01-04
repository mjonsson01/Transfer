// File: Transfer/src/Entities/PhysicsStructures.h

#pragma once

// Standard Library Imports
#include <cmath>
#include <iostream>

struct Vector2D
{
	double xVal;
	double yVal;

	// Constructor
	Vector2D(): xVal(0.0), yVal(0.0){};
	Vector2D(double x, double y): xVal(x), yVal(y){};

	// --- Vector-Vector Arithmetic operators --- //
    Vector2D operator+(const Vector2D& other) const {
        return {xVal + other.xVal, yVal + other.yVal};
    }
    Vector2D operator-(const Vector2D& other) const {
        return {xVal - other.xVal, yVal - other.yVal};
    }
    
    // --- Vector-Scalar Arithmetic operators --- //
    Vector2D operator*(double scalar) const {
        return {xVal * scalar, yVal * scalar};
    }
    Vector2D operator/(double scalar) const {
        return {xVal / scalar, yVal / scalar};
    }
	
    // --- Vector-Vector Assignment operations --- //
    Vector2D& operator+=(const Vector2D& other) {
        xVal += other.xVal;
        yVal += other.yVal;
        return *this;
    }
    Vector2D& operator-=(const Vector2D& other) {
        xVal -= other.xVal;
        yVal -= other.yVal;
        return *this;
    }

    // --- Vector-Scalar Assignment operations --- //
    Vector2D& operator*=(double scalar) {
        xVal *= scalar;
        yVal *= scalar;
        return *this;
    }
    Vector2D& operator/=(double scalar) {
        xVal /= scalar;
        yVal /= scalar;
        return *this;
    }
	
    // --- Special Vector Utilities --- //
    double magnitude() const {
        return sqrt(xVal * xVal + yVal * yVal);
    }
	double square_magnitude() const {
		return (xVal * xVal + yVal * yVal);
	}
    double dot(const Vector2D& other){
        return xVal*other.xVal + yVal*other.yVal;
    }
	Vector2D& normalizeInPlace() {
    	double mag = magnitude();
    	if (mag != 0.0) {
        	xVal /= mag;
        	yVal /= mag;
   		}
		return *this;
	}
    Vector2D normalize() {
        double mag = magnitude();
    	if (mag != 0.0) {
        	xVal /= mag;
        	yVal /= mag;
   		}
		return Vector2D(xVal, yVal);
    }
    
};

// --------- SPECIALTY METHODS --------- //

// Linear interpolation betweeen two vectors
static inline Vector2D lerp(const Vector2D& a, const Vector2D& b, double t)
{
    return a + (b - a) * t;
}

// --------- I/O OPERATOR OVERLOAD --------- //

// Print operator overload
inline std::ostream& operator<<(std::ostream& os, const Vector2D& Vec) {
    os << "{ "<< Vec.xVal << ", " << Vec.yVal << " }";
    return os;
}