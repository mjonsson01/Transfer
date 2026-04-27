// File: Transfer/src/Utilities/CustomMathUtilities.h

// Standard Library Imports
#include <random>

// Custom Imports
#include "Utilities/Constants/EngineConstants.h"

#pragma once

static double randomDouble(double minVal, double maxVal)
{
    static thread_local std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<double> dist(minVal, maxVal);
    return dist(rng);
}

static bool firstWithinEpsilonOfSecond(double valueToCheck, double valueToCompareTo)
{
    return ((valueToCheck <= valueToCompareTo + EPSILON) && (valueToCheck >= valueToCompareTo - EPSILON));
}