// File: Transfer/src/Utilities/CustomMathUtilities.hpp

// Custom Imports
#include "DynamoEngine/Constants/GlobalConstants.hpp"
#include "Constants/EngineConstants.hpp"

#pragma once

// (randomDouble lived here until 2026-11-27: a time-seeded generator made every shatter different on every run.
// Physics randomness is now hashed from the universe's seed: DynamoEngine/Math/DeterministicRandom.hpp.)

static bool firstWithinEpsilonOfSecond(double valueToCheck, double valueToCompareTo)
{
    return ((valueToCheck <= valueToCompareTo + DynamoEngine::EPSILON) &&
            (valueToCheck >= valueToCompareTo - DynamoEngine::EPSILON));
}