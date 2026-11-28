// File: Tests/DynamoEngine/Math/Test_DeterministicMath.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/Constants/GlobalConstants.hpp"
#include "DynamoEngine/Math/DeterministicMath.hpp"

// Standard Library Imports
#include <cmath>
#include <limits>

using namespace DynamoEngine;

// --------- HELPERS --------- //

// How far apart two doubles are, in units of the last place (ulps) of the expected value
static double ulpsApart(double actual, double expected)
{
    double ulp = std::nextafter(std::abs(expected), std::numeric_limits<double>::infinity()) - std::abs(expected);
    return std::abs(actual - expected) / ulp;
}

// --------- TESTS --------- //

TEST(DeterministicMath, SinAndCosOfSpecialAngles)
{
    EXPECT_EQ(deterministicSin(0.0), 0.0);
    EXPECT_EQ(deterministicCos(0.0), 1.0);
    // PI / 3 etc. are only the nearest doubles to those angles, so the answers can be a bit off in the last place
    // (std::cos(PI / 3) is 0.50000000000000011 too)
    EXPECT_NEAR(deterministicSin(PI / 2.0), 1.0, 3e-16);
    EXPECT_NEAR(deterministicCos(PI), -1.0, 3e-16);
    EXPECT_NEAR(deterministicSin(PI / 6.0), 0.5, 3e-16);
    EXPECT_NEAR(deterministicCos(PI / 3.0), 0.5, 3e-16);
}

TEST(DeterministicMath, SinAndCosMatchTheStandardOnesClosely)
{
    // Every 0.001 rad over several turns, both signs: within 1e-15 of the standard library (absolute: near the zero
    // crossings an "ulp" comparison would be meaningless)
    for (double x = -20.0; x <= 20.0; x += 0.001)
    {
        ASSERT_NEAR(deterministicSin(x), std::sin(x), 1e-15) << "x = " << x;
        ASSERT_NEAR(deterministicCos(x), std::cos(x), 1e-15) << "x = " << x;
    }
}

TEST(DeterministicMath, SinAndCosStayAccurateForLargeAngles)
{
    // A ship that has turned a lot: up to a million radians, still within 1e-12
    for (double x = 1000.0; x <= 1.0e6; x *= 1.37)
    {
        EXPECT_NEAR(deterministicSin(x), std::sin(x), 1e-12) << "x = " << x;
        EXPECT_NEAR(deterministicCos(x), std::cos(x), 1e-12) << "x = " << x;
    }
}

TEST(DeterministicMath, SinSquaredPlusCosSquaredIsOne)
{
    for (double x = -7.0; x <= 7.0; x += 0.01)
    {
        double s = deterministicSin(x);
        double c = deterministicCos(x);
        EXPECT_NEAR(s * s + c * c, 1.0, 4e-16) << "x = " << x;
    }
}

TEST(DeterministicMath, SinAndCosOfInfinityAndNaNAreNaN)
{
    double infinity = std::numeric_limits<double>::infinity();
    EXPECT_TRUE(std::isnan(deterministicSin(infinity)));
    EXPECT_TRUE(std::isnan(deterministicCos(-infinity)));
    EXPECT_TRUE(std::isnan(deterministicSin(std::nan(""))));
}

TEST(DeterministicMath, CbrtOfPerfectCubes)
{
    EXPECT_DOUBLE_EQ(deterministicCbrt(8.0), 2.0);
    EXPECT_DOUBLE_EQ(deterministicCbrt(27.0), 3.0);
    EXPECT_DOUBLE_EQ(deterministicCbrt(-8.0), -2.0);
    EXPECT_DOUBLE_EQ(deterministicCbrt(1.0), 1.0);
    EXPECT_EQ(deterministicCbrt(0.0), 0.0);
}

TEST(DeterministicMath, CbrtMatchesTheStandardOneWithinTwoUlps)
{
    // From tiny to huge (handleAccretion uses it on mass ratios, 0.9 .. 2, but it must work everywhere)
    for (double x = 1e-300; x < 1e300; x *= 1.093)
    {
        ASSERT_LE(ulpsApart(deterministicCbrt(x), std::cbrt(x)), 2.0) << "x = " << x;
        ASSERT_LE(ulpsApart(deterministicCbrt(-x), std::cbrt(-x)), 2.0) << "x = " << -x;
    }
    for (double x = 0.5; x <= 4.0; x += 0.0001)
    {
        ASSERT_LE(ulpsApart(deterministicCbrt(x), std::cbrt(x)), 2.0) << "x = " << x;
    }
}
