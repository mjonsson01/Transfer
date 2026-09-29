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
    T x_val = T(0);
    T y_val = T(0);

    // --- Constructors --- //
    constexpr Vector2() = default;
    constexpr Vector2(T x, T y) : x_val(x), y_val(y) {}

    // Conversion for double -> float Vector2F(Vector2D)
    template <std::floating_point U>
    constexpr explicit(sizeof(U) > sizeof(T)) Vector2(const Vector2<U>& other)
        : x_val(static_cast<T>(other.x_val)), y_val(static_cast<T>(other.y_val))
    {
    }
    // Conversion for float -> double is implicit DynamoEngine::Vector2D(Vector2F)

    // --- Vector-Vector Arithmetic operators --- //
    constexpr Vector2 operator+(const Vector2& other) const { return {x_val + other.x_val, y_val + other.y_val}; }
    constexpr Vector2 operator-(const Vector2& other) const { return {x_val - other.x_val, y_val - other.y_val}; }

    // --- Vector-Scalar Arithmetic operators --- //

    // V * scalar
    constexpr Vector2 operator*(T scalar) const { return {x_val * scalar, y_val * scalar}; }
    // scalar * V
    friend constexpr Vector2 operator*(T scalar, Vector2& vec) { return vec * scalar; }
    // V / scalar
    constexpr Vector2 operator/(T scalar) const { return {x_val / scalar, y_val / scalar}; }

    // --- Vector-Vector Assignment operators --- //

    constexpr Vector2& operator+=(const Vector2& other)
    {
        x_val += other.x_val;
        y_val += other.y_val;
        return *this;
    }
    constexpr Vector2& operator-=(const Vector2& other)
    {
        x_val -= other.x_val;
        y_val -= other.y_val;
        return *this;
    }

    // --- Vector-Scalar Assignment operators --- //
    constexpr Vector2& operator*=(T scalar)
    {
        x_val *= scalar;
        y_val *= scalar;
        return *this;
    }
    constexpr Vector2& operator/=(T scalar)
    {
        x_val /= scalar;
        y_val /= scalar;
        return *this;
    }

    // --- Special Vector Utilities --- //

    constexpr T magnitude() const { return std::sqrt(x_val * x_val + y_val * y_val); }
    constexpr T squareMagnitude() const { return (x_val * x_val + y_val * y_val); }
    constexpr T dot(const Vector2& other) const { return x_val * other.x_val + y_val * other.y_val; }

    // Normalize to unit vector with zero vector left unchanged
    Vector2& normalizeInPlace()
    {
        T mag = magnitude();
        if (mag != T(0))
        {
            x_val /= mag;
            y_val /= mag;
        }
        return *this;
    }
    Vector2 normalize() const
    {
        Vector2 copy = *this;
        return copy.normalizeInPlace();
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
    os << "{ " << vec.x_val << ", " << vec.y_val << " }";
    return os;
}
} // namespace DynamoEngine
