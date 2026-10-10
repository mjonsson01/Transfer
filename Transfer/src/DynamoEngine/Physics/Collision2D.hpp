// File: Transfer/src/DynamoEngine/Physics/Collision2D.hpp
#pragma once

// Custom Imports
#include "DynamoEngine/Math/Vector2.hpp"

// Standard Library Imports
#include <vector>

namespace DynamoEngine
{
// Result of testing a circle against a convex polygon.
struct CircleContact
{
    bool touching = false; // when false, normal and depth mean nothing
    Vector2D normal;       // unit vector from the polygon toward the circle: push the circle this way to separate them
    double depth = 0.0;    // how far they overlap along normal; moving them apart by this makes them just touch
};

// Tests a circle against a CONVEX polygon given by its corners in order around the outline (either direction),
// in the same space as the circle (e.g. world space). Needs at least 3 corners.
CircleContact circleVsConvexPolygon(const Vector2D& circle_center, double circle_radius,
                                    const std::vector<Vector2D>& polygon);
} // namespace DynamoEngine