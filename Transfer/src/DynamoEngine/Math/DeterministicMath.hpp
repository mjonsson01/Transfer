// File: Transfer/src/DynamoEngine/Math/DeterministicMath.hpp

#pragma once

namespace DynamoEngine
{
// Math functions that give bit-identical results on every machine.
//
// The C++ standard only requires +, -, *, / and sqrt to be exactly rounded (IEEE 754), so those match everywhere.
// std::sin, std::cos, std::pow, std::cbrt are left to each platform's math library, and Apple's and Microsoft's can
// disagree in the last bit. In physics a last-bit difference grows tick by tick until the outcome differs, so physics
// code uses these instead: they are built only from +, -, *, /, and exact helpers (floor, frexp, ldexp), in a fixed
// order. (That also needs the compiler to keep that order: no fused multiply-add, see CMakeLists.txt.)
//
// Accuracy: within a couple of units in the last place of the standard functions. sin and cos are accurate for
// angles up to about a million radians; bigger angles still give repeatable results, just less accurate ones.

// sin(x), x in radians
double deterministicSin(double x);
// cos(x), x in radians
double deterministicCos(double x);
// The cube root (x^(1/3)), negative x included: deterministicCbrt(-8) = -2
double deterministicCbrt(double x);
} // namespace DynamoEngine
