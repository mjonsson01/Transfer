// File: Tests/TemplateTests.cpp
// Template / smoke test: proves googletest is wired up. Copy this file's structure for real tests,
// mirroring the src/ path (e.g. src/Engine/Input/InputDevices.cpp -> Tests/Engine/Input/InputDevicesTests.cpp).

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
// #include "Engine/Input/InputDevices.hpp"   // code under test goes here (include paths are relative to Transfer/src)

// Standard Library Imports
#include <vector>

// --------- PLAIN TESTS --------- //
// TEST(SuiteName, TestName): suite = the thing being tested, name = the behavior being checked.
// Avoid underscores in either name (googletest reserves them internally).

TEST(Template, BasicAssertionsPass)
{
    // EXPECT_* records a failure and KEEPS GOING (use by default, so you see every failure at once).
    EXPECT_EQ(2 + 2, 4);
    EXPECT_TRUE(true);
    EXPECT_NE(1, 2);

    // ASSERT_* records a failure and STOPS this test (use when later lines would crash if it failed).
    std::vector<int> values = {1, 2, 3};
    ASSERT_EQ(values.size(), 3u);
    EXPECT_EQ(values[2], 3); // safe: the ASSERT above guarantees index 2 exists
}

TEST(Template, FloatingPointComparisons)
{
    // Never EXPECT_EQ on doubles -- rounding makes exact equality fragile.
    double momentum_velocity = 500.0 / 110.0;
    EXPECT_DOUBLE_EQ(momentum_velocity, 500.0 / 110.0); // equal within 4 ULPs (float noise)
    EXPECT_NEAR(momentum_velocity, 4.545, 1e-3);        // equal within an explicit tolerance
}

// --------- FIXTURE TESTS --------- //
// A fixture gives every test in the suite a fresh copy of the same setup.
// Each TEST_F gets a brand new TemplateFixture object, so tests can't leak state into each other.

class TemplateFixture : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // Runs before EACH test. Build the object(s) under test here.
        numbers = {10, 20, 30};
    }

    void TearDown() override
    {
        // Runs after EACH test. Usually empty: members are destroyed automatically.
    }

    std::vector<int> numbers;
};

TEST_F(TemplateFixture, StartsWithThreeNumbers) { EXPECT_EQ(numbers.size(), 3u); }

TEST_F(TemplateFixture, ModificationsDoNotLeakBetweenTests)
{
    numbers.push_back(40);
    EXPECT_EQ(numbers.size(), 4u); // the next test still sees 3 -- SetUp runs fresh each time
}