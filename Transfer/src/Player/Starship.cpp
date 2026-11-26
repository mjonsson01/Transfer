// File: Transfer/src/Player/Starship.cpp
#include "Starship.hpp"

#include "DynamoEngine/Constants/GlobalConstants.hpp"

Starship::Starship() { m_ship_size = 150.0; }

Starship::~Starship() {}

DynamoEngine::Vector2D Starship::getPointingVector()
{
    // rotation = 0 → pointing up (-y), matches your current nose placement
    return DynamoEngine::Vector2D{std::sin(m_rotation), -std::cos(m_rotation)};
}

void Starship::buildGeometry(std::vector<StarshipVertex>& starshipVertexBuffer)
{
    // The sprite is one square of side shipSize (two triangles), centered where the old placeholder box was and
    // rotated around that center. Every corner is sent at its current AND its previous-tick position, so the vertex
    // shader can interpolate between them (render alpha) like everything else in the world.
    float half = m_ship_size * 0.5f;
    float center_x = float(m_position.x_val) + half;
    float center_y = float(m_position.y_val) + half;
    float prev_center_x = float(m_prev_position.x_val) + half;
    float prev_center_y = float(m_prev_position.y_val) + half;

    // Each corner: its offset from the center before rotating, and which point of the image it shows.
    // u runs left -> right and v top -> bottom (0..1), so v = 0 is the PNG's top row: the nose.
    struct Corner
    {
        float dx, dy, u, v;
    };
    const Corner corners[4] = {
        {-half, -half, 0.0f, 0.0f}, // top-left
        {half, -half, 1.0f, 0.0f},  // top-right
        {half, half, 1.0f, 1.0f},   // bottom-right
        {-half, half, 0.0f, 1.0f},  // bottom-left
    };

    // The square as two triangles: top-left, top-right, bottom-right, then top-left, bottom-right, bottom-left
    const int triangle_corners[6] = {0, 1, 2, 0, 2, 3};
    for (int corner_index : triangle_corners)
    {
        const Corner& corner = corners[corner_index];
        DynamoEngine::Vector2D rotated = rotatedOffset({corner.dx, corner.dy});
        float rotated_dx = float(rotated.x_val);
        float rotated_dy = float(rotated.y_val);
        starshipVertexBuffer.push_back({center_x + rotated_dx, center_y + rotated_dy, prev_center_x + rotated_dx,
                                        prev_center_y + rotated_dy, m_ship_size, corner.u, corner.v, 1.0f, 1.0f, 1.0f,
                                        1.0f});
    }
}

DynamoEngine::Vector2D Starship::rotatedOffset(DynamoEngine::Vector2D offset) const
{
    double cos_r = std::cos(m_rotation);
    double sin_r = std::sin(m_rotation);
    return {offset.x_val * cos_r - offset.y_val * sin_r, offset.x_val * sin_r + offset.y_val * cos_r};
}

std::vector<DynamoEngine::Vector2D> Starship::collisionPolygon() const
{
    std::vector<DynamoEngine::Vector2D> corners;
    corners.reserve(6);
    for (const DynamoEngine::Vector2D& corner_px : HITBOX_CORNERS_PX)
    {
        // Image pixels -> offset from the sprite's center, in world units: the image's center is (0, 0), and the
        // whole image spans shipSize, exactly like the sprite's square
        DynamoEngine::Vector2D offset =
            (corner_px / HITBOX_IMAGE_SIZE_PX - DynamoEngine::Vector2D(0.5, 0.5)) * double(m_ship_size);
        corners.push_back(center() + rotatedOffset(offset)); // same center as the sprite
    }
    return corners;
}

void Starship::applyThrust(const ShipControls& controls)
{
    adjustThrottle(controls.thrust, PHYSICS_TIME_STEP); // one tick of game time
}

void Starship::adjustThrottle(int direction, double seconds)
{
    // W pushes the throttle up, S pulls it down, at THROTTLE_CHANGE_PER_SECOND; with neither it stays where it is.
    // No reverse thrust: S only goes down to 0 (to brake, turn around and burn).
    m_throttle += direction * THROTTLE_CHANGE_PER_SECOND * seconds;
    m_throttle = std::clamp(m_throttle, 0.0, 1.0);

    // An ACCELERATION along the nose, not a jump in velocity: the speed it adds depends on how long it's held,
    // and it's integrated together with gravity
    m_thrust_acceleration = getPointingVector() * (m_throttle * THRUST_ACCELERATION);
}

void Starship::applyRotation(const ShipControls& controls)
{
    const double turn_speed = 0.05;

    // NOTE: signs are reversed on purpose because the y origin is in the upper left instead of lower left.
    if (controls.rotation > 0)
        m_rotation -= turn_speed;
    if (controls.rotation < 0)
        m_rotation += turn_speed;

    // Optional: keep rotation in a sane range to avoid float drift over time
    m_rotation = std::remainder(m_rotation, DynamoEngine::TWO_PI); // result in [-π, π]
}

void Starship::applyVelocityVerletFirstHalfStep()
{
    m_prev_position = m_position;
    m_velocity += m_acceleration * (PHYSICS_TIME_STEP / 2.0); // first half kick, with last tick's acceleration
    m_position += m_velocity * PHYSICS_TIME_STEP;             // drift
}

void Starship::setGravityForce(const DynamoEngine::Vector2D& gravitational_force)
{
    // a = F / m for gravity, plus the thrust's acceleration
    m_acceleration = gravitational_force / m_mass + m_thrust_acceleration;
}

void Starship::applyVelocityVerletSecondHalfStep()
{
    m_velocity += m_acceleration * (PHYSICS_TIME_STEP / 2.0); // second half kick, with the new acceleration
}
void Starship::applyImpulse(const DynamoEngine::Vector2D& impulse) { m_velocity += impulse / m_mass; }

void Starship::moveBy(const DynamoEngine::Vector2D& offset)
{
    // Only the current position: the renderer then glides from the previous position to the corrected one
    m_position += offset;
}

void Starship::bounceInsideWalls(const DynamoEngine::Vector2D& box_min, const DynamoEngine::Vector2D& box_max,
                                 double restitution)
{
    const double reach = boundingRadius();
    DynamoEngine::Vector2D center_now = center();
    DynamoEngine::Vector2D push = {0.0, 0.0};
    if (center_now.x_val - reach < box_min.x_val)
    {
        push.x_val = box_min.x_val + reach - center_now.x_val;
        if (m_velocity.x_val < 0.0)
            m_velocity.x_val = -m_velocity.x_val * restitution;
    }
    else if (center_now.x_val + reach > box_max.x_val)
    {
        push.x_val = box_max.x_val - reach - center_now.x_val;
        if (m_velocity.x_val > 0.0)
            m_velocity.x_val = -m_velocity.x_val * restitution;
    }
    if (center_now.y_val - reach < box_min.y_val)
    {
        push.y_val = box_min.y_val + reach - center_now.y_val;
        if (m_velocity.y_val < 0.0)
            m_velocity.y_val = -m_velocity.y_val * restitution;
    }
    else if (center_now.y_val + reach > box_max.y_val)
    {
        push.y_val = box_max.y_val - reach - center_now.y_val;
        if (m_velocity.y_val > 0.0)
            m_velocity.y_val = -m_velocity.y_val * restitution;
    }
    m_position += push;
}

void Starship::bounceOffBox(const DynamoEngine::Vector2D& box_min, const DynamoEngine::Vector2D& box_max,
                            double restitution)
{
    const DynamoEngine::Vector2D previous_center = center() - movementThisTick();
    const bool was_inside = previous_center.x_val >= box_min.x_val && previous_center.x_val <= box_max.x_val &&
                            previous_center.y_val >= box_min.y_val && previous_center.y_val <= box_max.y_val;
    if (was_inside)
    {
        bounceInsideWalls(box_min, box_max, restitution);
        return;
    }

    // From outside: if the bounding circle now overlaps the box, push it back out through the side it came from
    const double reach = boundingRadius();
    const DynamoEngine::Vector2D center_now = center();
    const bool overlaps = center_now.x_val + reach > box_min.x_val && center_now.x_val - reach < box_max.x_val &&
                          center_now.y_val + reach > box_min.y_val && center_now.y_val - reach < box_max.y_val;
    if (!overlaps)
    {
        return;
    }
    DynamoEngine::Vector2D push = {0.0, 0.0};
    if (previous_center.x_val < box_min.x_val)
    {
        push.x_val = box_min.x_val - reach - center_now.x_val;
        if (m_velocity.x_val > 0.0)
            m_velocity.x_val = -m_velocity.x_val * restitution;
    }
    else if (previous_center.x_val > box_max.x_val)
    {
        push.x_val = box_max.x_val + reach - center_now.x_val;
        if (m_velocity.x_val < 0.0)
            m_velocity.x_val = -m_velocity.x_val * restitution;
    }
    else if (previous_center.y_val < box_min.y_val)
    {
        push.y_val = box_min.y_val - reach - center_now.y_val;
        if (m_velocity.y_val > 0.0)
            m_velocity.y_val = -m_velocity.y_val * restitution;
    }
    else
    {
        push.y_val = box_max.y_val + reach - center_now.y_val;
        if (m_velocity.y_val < 0.0)
            m_velocity.y_val = -m_velocity.y_val * restitution;
    }
    m_position += push;
}

double Starship::boundingRadius() const
{
    double farthest = 0.0;
    for (const DynamoEngine::Vector2D& corner_px : HITBOX_CORNERS_PX)
    {
        // Same pixels -> world conversion as collisionPolygon (rotation doesn't change a distance from the center)
        DynamoEngine::Vector2D offset =
            (corner_px / HITBOX_IMAGE_SIZE_PX - DynamoEngine::Vector2D(0.5, 0.5)) * double(m_ship_size);
        farthest = std::max(farthest, offset.magnitude());
    }
    return farthest;
}
DynamoEngine::Vector2D Starship::center() const
{
    // The ship's position is the top-left corner of its (unrotated) square
    return m_position + DynamoEngine::Vector2D(halfSize(), halfSize());
}