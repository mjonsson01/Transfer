#pragma once
#include <cmath>

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

	Vector2D& normalize() {
    	double mag = magnitude();
    	if (mag != 0.0) {
        	x_val /= mag;
        	y_val /= mag;
   		}
		return *this;
	}
};

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


// struct InitializerVelocities
// {
// 	int x_init;
// 	int y_init;
// 	int x_end;
// 	int y_end;
// 	VelocityVector2D velocity;
// };