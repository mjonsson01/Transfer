// File: Tests/DynamoEngine/Physics/Test_Collision2D.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/Math/Vector2.hpp"
#include "DynamoEngine/Physics/Collision2D.hpp"
#include "TestingUtilities/VectorsNear.hpp"

// Standard Library Imports
#include <cmath>
#include <vector>

using namespace DynamoEngine;

// --------- TEST HELPERS --------- //

// No SetUp needed: circleVsConvexPolygon is a free function with no state. The fixture only holds the shared shapes.
class Collision2DTest : public ::testing::Test
{
  protected:
    // A 10 x 10 square with its top-left corner at the origin (y grows downward, so y = 0 is the TOP edge)
    static std::vector<Vector2D> square() { return {{0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}}; }

    // The same square with its corners in the opposite order, to check the winding doesn't change any result
    static std::vector<Vector2D> squareReversed() { return {{0.0, 10.0}, {10.0, 10.0}, {10.0, 0.0}, {0.0, 0.0}}; }

    // A right triangle whose slanted edge runs from (10, 0) to (0, 10)
    static std::vector<Vector2D> triangle() { return {{0.0, 0.0}, {10.0, 0.0}, {0.0, 10.0}}; }
};

TEST_F(Collision2DTest, FarAwayIsNotTouching)
{
    CircleContact contact = circleVsConvexPolygon(Vector2D(30.0, 5.0), 2.0, square());

    EXPECT_FALSE(contact.touching);
}

TEST_F(Collision2DTest, SmallGapIsNotTouching)
{
    CircleContact contact = circleVsConvexPolygon(Vector2D(12.5, 5.0), 2.0, square());

    EXPECT_FALSE(contact.touching);
}

TEST_F(Collision2DTest, ExactlyTouchingIsNotOverlapping)
{
    CircleContact contact = circleVsConvexPolygon(Vector2D(12.0, 5.0), 2.0, square());

    EXPECT_FALSE(contact.touching);
}

TEST_F(Collision2DTest, OverlapOnRightEdge)
{
    CircleContact contact = circleVsConvexPolygon(Vector2D(11.0, 5.0), 2.0, square());

    EXPECT_TRUE(contact.touching);
    EXPECT_TRUE(VectorsNear(contact.normal, Vector2D(1.0, 0.0), 1e-9));
    EXPECT_NEAR(contact.depth, 1.0, 1e-9);
}

TEST_F(Collision2DTest, OverlapOnCornerPointsDiagonally)
{
    CircleContact contact = circleVsConvexPolygon(Vector2D(11.0, 11.0), 2.0, square());

    EXPECT_TRUE(contact.touching);
    EXPECT_TRUE(VectorsNear(contact.normal, Vector2D(std::sqrt(2.0) / 2.0, std::sqrt(2.0) / 2.0), 1e-9));
    EXPECT_NEAR(contact.depth, 2.0 - std::sqrt(2.0), 1e-9);
}

TEST_F(Collision2DTest, CenterInsideNearTopEdge)
{
    CircleContact contact = circleVsConvexPolygon(Vector2D(5.0, 1.0), 2.0, square());

    EXPECT_TRUE(contact.touching);
    EXPECT_TRUE(VectorsNear(contact.normal, Vector2D(0.0, -1.0), 1e-9));
    EXPECT_NEAR(contact.depth, 3.0, 1e-9);
}

TEST_F(Collision2DTest, WindingDoesNotMatter)
{
    // Same circles as OverlapOnRightEdge and CenterInsideNearTopEdge, corners in the opposite order: same answers
    CircleContact outside = circleVsConvexPolygon(Vector2D(11.0, 5.0), 2.0, squareReversed());

    EXPECT_TRUE(outside.touching);
    EXPECT_TRUE(VectorsNear(outside.normal, Vector2D(1.0, 0.0), 1e-9));
    EXPECT_NEAR(outside.depth, 1.0, 1e-9);

    // The inside case uses the edge's outward normal directly, so it's the one a wrong flip would break
    CircleContact inside = circleVsConvexPolygon(Vector2D(5.0, 1.0), 2.0, squareReversed());

    EXPECT_TRUE(inside.touching);
    EXPECT_TRUE(VectorsNear(inside.normal, Vector2D(0.0, -1.0), 1e-9));
    EXPECT_NEAR(inside.depth, 3.0, 1e-9);
}

TEST_F(Collision2DTest, TriangleSlantedEdge)
{
    // The slanted edge runs from (10, 0) to (0, 10): the line x + y = 10. Its outward direction is (1, 1) / sqrt(2).
    // The centre (6, 6) is (6 + 6 - 10) / sqrt(2) = sqrt(2) beyond it, so the radius-2 circle overlaps by 2 - sqrt(2).
    CircleContact contact = circleVsConvexPolygon(Vector2D(6.0, 6.0), 2.0, triangle());

    const double half_root_two = std::sqrt(2.0) / 2.0;
    EXPECT_TRUE(contact.touching);
    EXPECT_TRUE(VectorsNear(contact.normal, Vector2D(half_root_two, half_root_two), 1e-9));
    EXPECT_NEAR(contact.depth, 2.0 - std::sqrt(2.0), 1e-9);
}

// EXPECT_DEBUG_DEATH: in Debug builds the call must stop the program (the assert fires) with a message containing
// the given text; in Release builds asserts are compiled out, so the call just has to not crash the test
TEST_F(Collision2DTest, FewerThanThreeCornersAsserts)
{
    std::vector<Vector2D> just_a_line = {{0.0, 0.0}, {10.0, 0.0}};

    EXPECT_DEBUG_DEATH(circleVsConvexPolygon(Vector2D(5.0, 5.0), 2.0, just_a_line), "at least 3 corners");
}
