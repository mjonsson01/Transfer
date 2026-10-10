// File: Transfer/src/DynamoEngine/Physics/Collision2D.cpp
#include "DynamoEngine/Physics/Collision2D.hpp"

// Custom Imports
#include "DynamoEngine/Math/Vector2.hpp"

// Standard Library Imports
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace DynamoEngine
{
namespace
{
// The point on the segment from a to b that is closest to point p
Vector2D closestPointOnSegment(const Vector2D& p, const Vector2D& a, const Vector2D& b)
{
    Vector2D edge = b - a;
    double edge_length_squared = edge.squareMagnitude();
    if (edge_length_squared == 0.0)
    {
        return a; // a zero-length edge is just a point
    }

    // How far along the edge the closest point lies: 0 = at a, 1 = at b. Clamped, so it never leaves the segment.
    double t = (p - a).dot(edge) / edge_length_squared;
    t = std::clamp(t, 0.0, 1.0);
    return a + edge * t;
}
} // namespace

CircleContact circleVsConvexPolygon(const Vector2D& circle_center, double circle_radius,
                                    const std::vector<Vector2D>& polygon)
{
    assert(polygon.size() >= 3 && "circleVsConvexPolygon needs at least 3 corners");

    // The average of the corners is always inside a convex polygon: it tells us which side of each edge is "inside"
    Vector2D polygon_center;
    for (const Vector2D& corner : polygon)
    {
        polygon_center += corner;
    }
    polygon_center /= static_cast<double>(polygon.size());

    bool center_is_inside = true; // until we find an edge it's on the outer side of
    double closest_distance_squared = std::numeric_limits<double>::infinity();
    Vector2D closest_point;
    Vector2D closest_edge_outward_normal;

    for (size_t i = 0; i < polygon.size(); ++i)
    {
        const Vector2D& a = polygon[i];
        const Vector2D& b = polygon[(i + 1) % polygon.size()]; // the last corner connects back to the first

        // This edge's outward normal: perpendicular to the edge, flipped if it points toward the polygon's centre.
        // Working it out from the centre is what makes either corner order (winding) work.
        Vector2D edge = b - a;
        Vector2D outward_normal = Vector2D(edge.y_val, -edge.x_val).normalize();
        if (outward_normal.dot(polygon_center - a) > 0.0)
        {
            outward_normal = outward_normal * -1.0;
        }

        // Convex polygon: the circle's centre is outside if it's on the outer side of ANY edge
        if (outward_normal.dot(circle_center - a) > 0.0)
        {
            center_is_inside = false;
        }

        Vector2D point_on_edge = closestPointOnSegment(circle_center, a, b);
        double distance_squared = (circle_center - point_on_edge).squareMagnitude();
        if (distance_squared < closest_distance_squared)
        {
            closest_distance_squared = distance_squared;
            closest_point = point_on_edge;
            closest_edge_outward_normal = outward_normal;
        }
    }

    double distance = std::sqrt(closest_distance_squared);
    CircleContact contact;

    if (center_is_inside)
    {
        // Deep hit: the centre has crossed into the polygon, so push it back out through the nearest edge.
        // It must travel back to that edge (distance) AND a full radius beyond it to stop overlapping.
        contact.touching = true;
        contact.normal = closest_edge_outward_normal;
        contact.depth = circle_radius + distance;
        return contact;
    }

    if (distance >= circle_radius)
    {
        return contact; // outside and not overlapping (exactly touching counts as not overlapping)
    }

    // Outside but overlapping: push away from the closest point. Near a corner this points diagonally,
    // so the polygon's corners behave as if they were rounded.
    contact.touching = true;
    contact.normal = (circle_center - closest_point) / distance;
    contact.depth = circle_radius - distance;
    return contact;
}
} // namespace DynamoEngine