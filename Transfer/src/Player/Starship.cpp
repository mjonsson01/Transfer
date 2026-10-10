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
    // The sprite is one square of side shipSize (two triangles), centred where the old placeholder box was and
    // rotated around that centre. Every corner is sent at its current AND its previous-tick position, so the vertex
    // shader can interpolate between them (render alpha) like everything else in the world.
    float half = m_ship_size * 0.5f;
    float center_x = float(m_position.x_val) + half;
    float center_y = float(m_position.y_val) + half;
    float prev_center_x = float(m_prev_position.x_val) + half;
    float prev_center_y = float(m_prev_position.y_val) + half;

    // Each corner: its offset from the centre before rotating, and which point of the image it shows.
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
        // Image pixels -> offset from the sprite's centre, in world units: the image's centre is (0, 0), and the
        // whole image spans shipSize, exactly like the sprite's square
        DynamoEngine::Vector2D offset =
            (corner_px / HITBOX_IMAGE_SIZE_PX - DynamoEngine::Vector2D(0.5, 0.5)) * double(m_ship_size);
        corners.push_back(center() + rotatedOffset(offset)); // same centre as the sprite
    }
    return corners;
}

void Starship::applyThrust(UIState& uiState)
{
    DEPRECATED_InputState& input_state = uiState.getMutableDEPRECATED_InputState();
    int sign = 0; // +1 = forward (W), -1 = backward (S), 0 = no thrust
    if (input_state.isRequestingThrust)
    {
        if (input_state.positiveThrust)
            sign = 1;
        if (input_state.negativeThrust)
            sign = -1;
    }

    // An ACCELERATION along the nose, not a jump in velocity: the speed it adds depends on how long it's held,
    // and it's integrated together with gravity
    m_thrust_acceleration = getPointingVector() * (sign * THRUST_ACCELERATION);
}

void Starship::applyRotation(UIState& uiState)
{
    DEPRECATED_InputState& input_state = uiState.getMutableDEPRECATED_InputState();
    const double turn_speed = 0.05;

    // NOTE: signs are reversed on purpose because the y origin is in the upper left instead of lower left.
    if (input_state.positiveRotation)
        m_rotation -= turn_speed;
    if (input_state.negativeRotation)
        m_rotation += turn_speed;

    // Optional: keep rotation in a sane range to avoid float drift over time
    m_rotation = std::remainder(m_rotation, DynamoEngine::TWO_PI); // result in [-π, π]
}

void Starship::applyVelocityVerletPhase1()
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

void Starship::applyVelocityVerletPhase2()
{
    m_velocity += m_acceleration * (PHYSICS_TIME_STEP / 2.0); // second half kick, with the new acceleration
}
void Starship::applyImpulse(const DynamoEngine::Vector2D& impulse) { m_velocity += impulse / m_mass; }

void Starship::moveBy(const DynamoEngine::Vector2D& offset)
{
    // Only the current position: the renderer then glides from the previous position to the corrected one
    m_position += offset;
}

double Starship::boundingRadius() const
{
    double farthest = 0.0;
    for (const DynamoEngine::Vector2D& corner_px : HITBOX_CORNERS_PX)
    {
        // Same pixels -> world conversion as collisionPolygon (rotation doesn't change a distance from the centre)
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