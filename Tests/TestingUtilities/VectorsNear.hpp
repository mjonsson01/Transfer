// File: Tests/TestingUtilities/VectorsNear.hpp

#pragma once

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/Math/Vector2.hpp"

// Standard Library Includes
#include <cmath>
#include <concepts>
#include <type_traits>

template <std::floating_point T>
::testing::AssertionResult VectorsNear(const DynamoEngine::Vector2<T>& actual,
                                       const std::type_identity_t<DynamoEngine::Vector2<T>>& expected,
                                       std::type_identity_t<T> tolerance)
{
    if (std::abs(actual.xVal - expected.xVal) <= tolerance && std::abs(actual.yVal - expected.yVal) <= tolerance)
    {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure() << "actual " << actual << " vs expected " << expected << " (tolerance "
                                         << tolerance << ")";
}