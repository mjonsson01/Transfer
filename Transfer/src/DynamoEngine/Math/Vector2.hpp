// File: Transfer/src/DynamoEngine/Math/Vector2.hpp

#pragma once

// Standard Library Imports
#include <cmath>
#include <concepts>
#include <iostream>
#include <type_traits>

namespace DynamoEngine
{
// 2D vector over floating-point types. Aliases for DynamoEngine::Vector2D (double) and Vector2F (float) exist

template <std::floating_point T> struct Vector2
{
    T xVal = T(0);
    T yVal = T(0);

    // --- Constructors --- //
    constexpr Vector2() = default;
    constexpr Vector2(T x, T y) : xVal(x), yVal(y) {}

    // Conversion for double -> float Vector2F(Vector2D)
    template <std::floating_point U>
    constexpr explicit(sizeof(U) > sizeof(T)) Vector2(const Vector2<U>& other)
        : xVal(static_cast<T>(other.xVal)), yVal(static_cast<T>(other.yVal))
    {
    }
    // Conversion for float -> double is implicit DynamoEngine::Vector2D(Vector2F)

    // --- Vector-Vector Arithmetic operators --- //
    constexpr Vector2 operator+(const Vector2& other) const { return {xVal + other.xVal, yVal + other.yVal}; }
    constexpr Vector2 operator-(const Vector2& other) const { return {xVal - other.xVal, yVal - other.yVal}; }

    // --- Vector-Scalar Arithmetic operators --- //

    // V * scalar
    constexpr Vector2 operator*(T scalar) const { return {xVal * scalar, yVal * scalar}; }
    // scalar * V
    friend constexpr Vector2 operator*(T scalar, Vector2& vec) { return vec * scalar; }
    // V / scalar
    constexpr Vector2 operator/(T scalar) const { return {xVal / scalar, yVal / scalar}; }

    // --- Vector-Vector Assignment operators --- //

    constexpr Vector2& operator+=(const Vector2& other)
    {
        xVal += other.xVal;
        yVal += other.yVal;
        return *this;
    }
    constexpr Vector2& operator-=(const Vector2& other)
    {
        xVal -= other.xVal;
        yVal -= other.yVal;
        return *this;
    }

    // --- Vector-Scalar Assignment operators --- //
    constexpr Vector2& operator*=(T scalar)
    {
        xVal *= scalar;
        yVal *= scalar;
        return *this;
    }
    constexpr Vector2& operator/=(T scalar)
    {
        xVal /= scalar;
        yVal /= scalar;
        return *this;
    }

    // --- Special Vector Utilities --- //

    constexpr T magnitude() const { return std::sqrt(xVal * xVal + yVal * yVal); }
    constexpr T square_magnitude() const { return (xVal * xVal + yVal * yVal); }
    constexpr T dot(const Vector2& other) const { return xVal * other.xVal + yVal * other.yVal; }

    // Normalize to unit vector with zero vector left unchanged
    Vector2& normalize_in_place()
    {
        T mag = magnitude();
        if (mag != T(0))
        {
            xVal /= mag;
            yVal /= mag;
        }
        return *this;
    }
    Vector2& normalize() const
    {
        Vector2 copy = *this;
        return copy.normalize_in_place();
    }
};

using Vector2D = Vector2<double>;
using Vector2F = Vector2<float>;

// --- Specialty methods --- //

// Linear interpolation between two vectors
template <std::floating_point T>
constexpr Vector2<T> lerp(const Vector2<T>& a, const Vector2<T>& b, std::type_identity_t<T> t)
{
    return a + (b - a) * t;
}

// --- I/O Operator Overload --- //
template <std::floating_point T> std::ostream& operator<<(std::ostream& os, const Vector2<T>& vec)
{
    os << "{ " << vec.xVal << ", " << vec.yVal << " }";
    return os;
}
} // namespace DynamoEngine
