// File: Transfer/src/DynamoEngine/Math/DeterministicMath.cpp
#include "DynamoEngine/Math/DeterministicMath.hpp"

// Standard Library Imports
#include <cmath>
#include <cstdint>

namespace DynamoEngine
{
namespace
{
// 2 / pi
constexpr double TWO_OVER_PI = 6.36619772367581382433e-01;
// pi / 2 split into three parts, each with few enough bits that (whole number) x (part) is EXACT for whole numbers up
// to about a million. Subtracting the parts one at a time removes whole multiples of pi / 2 from an angle without the
// rounding error a single multiply would leave. (The classic constants from fdlibm, Sun's math library.)
constexpr double PI_OVER_TWO_PART_1 = 1.57079632673412561417e+00;
constexpr double PI_OVER_TWO_PART_2 = 6.07710050630396597660e-11;
constexpr double PI_OVER_TWO_PART_3 = 2.02226624871116645580e-21;

// Taylor series of sin and cos around 0, used for |r| <= pi / 4. There the first term left out is below 1e-19, far
// under the last bit of a double:
//   sin(r) = r * (1 - r^2/3! + r^4/5! - ... + r^16/17!)
//   cos(r) =      1 - r^2/2! + r^4/4! - ... - r^18/18!
// Evaluated from the innermost bracket out (Horner's rule): series = series * r^2 + next coefficient. Always the same
// operations in the same order. (The coefficients are divisions of constants, worked out by the compiler, exactly
// rounded on every compiler.)
double sinNearZero(double r)
{
    double r2 = r * r;
    double series = 1.0 / 355687428096000.0;     // +1/17!
    series = series * r2 - 1.0 / 1307674368000.0; // -1/15!
    series = series * r2 + 1.0 / 6227020800.0;    // +1/13!
    series = series * r2 - 1.0 / 39916800.0;      // -1/11!
    series = series * r2 + 1.0 / 362880.0;        // +1/9!
    series = series * r2 - 1.0 / 5040.0;          // -1/7!
    series = series * r2 + 1.0 / 120.0;           // +1/5!
    series = series * r2 - 1.0 / 6.0;             // -1/3!
    series = series * r2 + 1.0;
    return r * series;
}

double cosNearZero(double r)
{
    double r2 = r * r;
    double series = -1.0 / 6402373705728000.0;     // -1/18!
    series = series * r2 + 1.0 / 20922789888000.0; // +1/16!
    series = series * r2 - 1.0 / 87178291200.0;    // -1/14!
    series = series * r2 + 1.0 / 479001600.0;      // +1/12!
    series = series * r2 - 1.0 / 3628800.0;        // -1/10!
    series = series * r2 + 1.0 / 40320.0;          // +1/8!
    series = series * r2 - 1.0 / 720.0;            // -1/6!
    series = series * r2 + 1.0 / 24.0;             // +1/4!
    series = series * r2 - 1.0 / 2.0;              // -1/2!
    series = series * r2 + 1.0;
    return series;
}

// Splits x into (quarter turns) and (what's left), x = quarter_turns * pi/2 + remainder, |remainder| <= pi/4.
// Returns the quarter turns modulo 4 (0..3): which quadrant the angle is in.
int reduceAngle(double x, double& remainder_out)
{
    double quarter_turns = std::floor(x * TWO_OVER_PI + 0.5); // the nearest whole number (floor is exact)
    remainder_out = ((x - quarter_turns * PI_OVER_TWO_PART_1) - quarter_turns * PI_OVER_TWO_PART_2) -
                    quarter_turns * PI_OVER_TWO_PART_3;
    // "& 3" keeps the last two bits: the quarter turns modulo 4, also for negative counts (-1 & 3 = 3)
    return static_cast<int>(static_cast<int64_t>(quarter_turns) & 3);
}
} // namespace

double deterministicSin(double x)
{
    if (!std::isfinite(x))
    {
        return x - x; // NaN for infinity and NaN, like std::sin
    }
    double r = 0.0;
    int quadrant = reduceAngle(x, r);
    // sin(r + 90 deg) = cos(r), sin(r + 180 deg) = -sin(r), sin(r + 270 deg) = -cos(r)
    switch (quadrant)
    {
    case 0:
        return sinNearZero(r);
    case 1:
        return cosNearZero(r);
    case 2:
        return -sinNearZero(r);
    default:
        return -cosNearZero(r);
    }
}

double deterministicCos(double x)
{
    if (!std::isfinite(x))
    {
        return x - x;
    }
    double r = 0.0;
    int quadrant = reduceAngle(x, r);
    // cos(r + 90 deg) = -sin(r), cos(r + 180 deg) = -cos(r), cos(r + 270 deg) = sin(r)
    switch (quadrant)
    {
    case 0:
        return cosNearZero(r);
    case 1:
        return -sinNearZero(r);
    case 2:
        return -cosNearZero(r);
    default:
        return sinNearZero(r);
    }
}

double deterministicCbrt(double x)
{
    if (x == 0.0 || !std::isfinite(x))
    {
        return x; // 0, -0, infinity and NaN are their own cube roots
    }
    double size = std::abs(x);

    // 1. Write size = fraction * 2^exponent with the exponent a multiple of 3, so the cube root of the 2^exponent part
    //    is exactly 2^(exponent / 3). frexp gives fraction in [0.5, 1); moving one or two powers of 2 into the
    //    fraction makes the exponent divisible by 3 and leaves the fraction in [0.5, 4). (frexp and ldexp only read
    //    and write the exponent bits: exact.)
    int exponent = 0;
    double fraction = std::frexp(size, &exponent);
    int leftover = ((exponent % 3) + 3) % 3; // 0, 1 or 2, also for negative exponents
    fraction = std::ldexp(fraction, leftover);
    exponent -= leftover;

    // 2. Newton's method for y^3 = fraction: y <- y - (y^3 - fraction) / (3 y^2) = (2y + fraction / y^2) / 3.
    //    A straight-line first guess is within ~25% on [0.5, 4); each step roughly squares the error, so 6 steps
    //    reach the last bit. Always 6 steps: a fixed amount of work, the same answer everywhere.
    double y = 0.6 + 0.4 * fraction;
    for (int step = 0; step < 6; ++step)
    {
        y = (2.0 * y + fraction / (y * y)) / 3.0;
    }

    double root = std::ldexp(y, exponent / 3);
    return (x < 0.0) ? -root : root;
}
} // namespace DynamoEngine
