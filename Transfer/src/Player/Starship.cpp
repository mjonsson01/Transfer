// File: Transfer/src/Player/Starship.cpp
#include "Starship.hpp"

Starship::Starship()
{
    velocity.xVal = 0;
    velocity.yVal = 0;
    shipSize = 50.0;
    position.xVal = 0.0;
    position.yVal = 0.0;
    prevPosition = position;
}

Starship::~Starship() {}

Vector2D Starship::getPointingVector()
{
    // rotation = 0 → pointing up (-y), matches your current nose placement
    return Vector2D{std::sin(rotation), -std::cos(rotation)};
}

void Starship::buildGeometry(std::vector<StarshipVertex>& starshipVertexBuffer)
{
    float x1 = float(position.xVal);
    float y1 = float(position.yVal);
    float x2 = x1 + shipSize;
    float y2 = y1 + shipSize;

    float prev_x1 = float(prevPosition.xVal);
    float prev_y1 = float(prevPosition.yVal);
    float prev_x2 = prev_x1 + shipSize;
    float prev_y2 = prev_y1 + shipSize;

    float noseHeight = shipSize * 0.5f;
    float x_nose = (x1 + x2) * 0.5f;
    float y_nose = y1 - noseHeight;
    float prev_x_nose = (prev_x1 + prev_x2) * 0.5f;
    float prev_y_nose = prev_y1 - noseHeight;

    // Rotate each point around its own center (current center for current points,
    // previous center for previous points, so interpolation stays sane)
    float cx = (x1 + x2) * 0.5f, cy = (y1 + y2) * 0.5f;
    float prev_cx = (prev_x1 + prev_x2) * 0.5f, prev_cy = (prev_y1 + prev_y2) * 0.5f;

    float cosR = float(std::cos(rotation));
    float sinR = float(std::sin(rotation));

    auto rotate = [](float px, float py, float ox, float oy, float c, float s) -> std::pair<float, float>
    {
        float dx = px - ox;
        float dy = py - oy;
        return {ox + dx * c - dy * s, oy + dx * s + dy * c};
    };

    auto [rx1, ry1] = rotate(x1, y1, cx, cy, cosR, sinR);
    auto [rx2y1_x, rx2y1_y] = rotate(x2, y1, cx, cy, cosR, sinR);
    auto [rx1y2_x, rx1y2_y] = rotate(x1, y2, cx, cy, cosR, sinR);
    auto [rx2, ry2] = rotate(x2, y2, cx, cy, cosR, sinR);
    auto [rnose_x, rnose_y] = rotate(x_nose, y_nose, cx, cy, cosR, sinR);

    auto [prx1, pry1] = rotate(prev_x1, prev_y1, prev_cx, prev_cy, cosR, sinR);
    auto [prx2y1_x, prx2y1_y] = rotate(prev_x2, prev_y1, prev_cx, prev_cy, cosR, sinR);
    auto [prx1y2_x, prx1y2_y] = rotate(prev_x1, prev_y2, prev_cx, prev_cy, cosR, sinR);
    auto [prx2, pry2] = rotate(prev_x2, prev_y2, prev_cx, prev_cy, cosR, sinR);
    auto [prnose_x, prnose_y] = rotate(prev_x_nose, prev_y_nose, prev_cx, prev_cy, cosR, sinR);

    float u = 0.0f, v = 0.0f, r = 1.0f, g = 1.0f, b = 1.0f, a = 1.0f;

    starshipVertexBuffer.push_back({rx1, ry1, prx1, pry1, shipSize, u, v, r, g, b, a});
    starshipVertexBuffer.push_back({rx2y1_x, rx2y1_y, prx2y1_x, prx2y1_y, shipSize, u, v, r, g, b, a});
    starshipVertexBuffer.push_back({rx1y2_x, rx1y2_y, prx1y2_x, prx1y2_y, shipSize, u, v, r, g, b, a});

    starshipVertexBuffer.push_back({rx2y1_x, rx2y1_y, prx2y1_x, prx2y1_y, shipSize, u, v, r, g, b, a});
    starshipVertexBuffer.push_back({rx2, ry2, prx2, pry2, shipSize, u, v, r, g, b, a});
    starshipVertexBuffer.push_back({rx1y2_x, rx1y2_y, prx1y2_x, prx1y2_y, shipSize, u, v, r, g, b, a});

    starshipVertexBuffer.push_back({rx1, ry1, prx1, pry1, shipSize, u, v, r, g, b, a});
    starshipVertexBuffer.push_back({rx2y1_x, rx2y1_y, prx2y1_x, prx2y1_y, shipSize, u, v, r, g, b, a});
    starshipVertexBuffer.push_back({rnose_x, rnose_y, prnose_x, prnose_y, shipSize, u, v, r, g, b, a});
}

void Starship::applyVelocity(UIState& uiState)
{
    InputState& input_state = uiState.getMutableInputState();
    if (input_state.isRequestingThrust)
    {
        int sign = 0;
        if (input_state.positiveThrust)
            sign = 1;
        if (input_state.negativeThrust)
            sign = -1;

        Vector2D direction = getPointingVector(); // unit vector, nose direction
        const double thrustMagnitude = 10.0;

        velocity.xVal += sign * direction.xVal * thrustMagnitude;
        velocity.yVal += sign * direction.yVal * thrustMagnitude;
    }
};

void Starship::applyRotation(UIState& uiState)
{
    InputState& input_state = uiState.getMutableInputState();
    const double turnSpeed = 0.05; // radians per tick — tune to taste

    // NOTE: signs are reversed on purpose because the y origin is in the upper left instead of lower left.
    if (input_state.positiveRotation)
        rotation -= turnSpeed;
    if (input_state.negativeRotation)
        rotation += turnSpeed;

    // Optional: keep rotation in a sane range to avoid float drift over time
    const double TWO_PI = 6.28318530718;
    if (rotation > TWO_PI)
        rotation -= TWO_PI;
    if (rotation < -TWO_PI)
        rotation += TWO_PI;
}

void Starship::integratePosition()
{
    prevPosition = position;
    position += velocity * PHYSICS_TIME_STEP;
}