// File: Tests/DynamoEngine/Math/Test_DeterministicRandom.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/Math/DeterministicRandom.hpp"

using namespace DynamoEngine;

// --------- TESTS --------- //

TEST(DeterministicRandom, SameInputsGiveTheSameNumber)
{
    uint64_t first = combineHash(combineHash(42, 7), 3);
    uint64_t second = combineHash(combineHash(42, 7), 3);
    EXPECT_EQ(first, second);
}

TEST(DeterministicRandom, KnownValuesNeverChange)
{
    // Pinned outputs: if these change, every saved level and replay would play out differently. Only change them on
    // purpose, together with the levels.
    EXPECT_EQ(mixBits(0), 0xE220A8397B1DCDAFULL);
    EXPECT_EQ(mixBits(1), 0x910A2DEC89025CC1ULL);
}

TEST(DeterministicRandom, InputOrderMatters)
{
    EXPECT_NE(combineHash(combineHash(1, 2), 3), combineHash(combineHash(1, 3), 2));
}

TEST(DeterministicRandom, UnitIntervalStaysInZeroToOne)
{
    EXPECT_EQ(hashToUnitInterval(0), 0.0);
    EXPECT_LT(hashToUnitInterval(0xFFFFFFFFFFFFFFFFULL), 1.0);
    for (uint64_t i = 0; i < 10000; ++i)
    {
        double value = hashToUnitInterval(mixBits(i));
        ASSERT_GE(value, 0.0);
        ASSERT_LT(value, 1.0);
    }
}

TEST(DeterministicRandom, RangeIsCoveredEvenly)
{
    // 10 equal buckets over [0.8, 1.1): each gets about a tenth of 100,000 draws
    int buckets[10] = {};
    for (uint64_t i = 0; i < 100000; ++i)
    {
        double value = hashToRange(combineHash(123, i), 0.8, 1.1);
        ASSERT_GE(value, 0.8);
        ASSERT_LT(value, 1.1);
        int bucket = static_cast<int>((value - 0.8) / 0.3 * 10.0);
        buckets[bucket < 10 ? bucket : 9] += 1;
    }
    for (int count : buckets)
    {
        EXPECT_NEAR(count, 10000, 400);
    }
}

TEST(DeterministicRandom, DoubleBitsTellApartValuesThatCompareEqual)
{
    EXPECT_NE(doubleBits(0.0), doubleBits(-0.0)); // equal as numbers, different bits
    EXPECT_EQ(doubleBits(1.5), doubleBits(1.5));
}
