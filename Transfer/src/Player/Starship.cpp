// File: Transfer/src/Player/Starship.cpp
#include "Starship.hpp"

#include "DynamoEngine/Constants/GlobalConstants.hpp"

Starship::Starship() { shipSize = 150.0; }

Starship::~Starship() {}

DynamoEngine::Vector2D Starship::getPointingVector()
{
    // rotation = 0 → pointing up (-y), matches your current nose placement
    return DynamoEngine::Vector2D{std::sin(rotation), -std::cos(rotation)};
}

void Starship::buildGeometry(std::vector<StarshipVertex>& starshipVertexBuffer)
{
    // The sprite is one square of side shipSize (two triangles), centred where the old placeholder box was and
    // rotated around that centre. Every corner is sent at its current AND its previous-tick position, so the vertex
    // shader can interpolate between them (render alpha) like everything else in the world.
    float half = shipSize * 0.5f;
    float center_x = float(position.x_val) + half;
    float center_y = float(position.y_val) + half;
    float prev_center_x = float(prevPosition.x_val) + half;
    float prev_center_y = float(prevPosition.y_val) + half;

    float cosR = float(std::cos(rotation));
    float sinR = float(std::sin(rotation));

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
        float rotated_dx = corner.dx * cosR - corner.dy * sinR;
        float rotated_dy = corner.dx * sinR + corner.dy * cosR;
        starshipVertexBuffer.push_back({center_x + rotated_dx, center_y + rotated_dy, prev_center_x + rotated_dx,
                                        prev_center_y + rotated_dy, shipSize, corner.u, corner.v, 1.0f, 1.0f, 1.0f,
                                        1.0f});
    }
}
void Starship::applyVelocity(UIState& uiState)
{
    DEPRECATED_InputState& input_state = uiState.getMutableDEPRECATED_InputState();
    if (input_state.isRequestingThrust)
    {
        int sign = 0;
        if (input_state.positiveThrust)
            sign = 1;
        if (input_state.negativeThrust)
            sign = -1;

        DynamoEngine::Vector2D direction = getPointingVector(); // unit vector, nose direction
        const double thrustMagnitude = 10.0;

        velocity.x_val += sign * direction.x_val * thrustMagnitude;
        velocity.y_val += sign * direction.y_val * thrustMagnitude;
    }
};

void Starship::applyRotation(UIState& uiState)
{
    DEPRECATED_InputState& input_state = uiState.getMutableDEPRECATED_InputState();
    const double turnSpeed = 0.05;

    // NOTE: signs are reversed on purpose because the y origin is in the upper left instead of lower left.
    if (input_state.positiveRotation)
        rotation -= turnSpeed;
    if (input_state.negativeRotation)
        rotation += turnSpeed;

    // Optional: keep rotation in a sane range to avoid float drift over time
    rotation = std::remainder(rotation, DynamoEngine::TWO_PI); // result in [-π, π]
}

void Starship::integratePosition()
{
    prevPosition = position;
    position += velocity * PHYSICS_TIME_STEP;
}