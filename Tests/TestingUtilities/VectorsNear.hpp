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
    if (std::abs(actual.x_val - expected.x_val) <= tolerance && std::abs(actual.y_val - expected.y_val) <= tolerance)
    {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure() << "actual " << actual << " vs expected " << expected << " (tolerance "
                                         << tolerance << ")";
}