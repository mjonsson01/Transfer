// File: Transfer/src/Player/Starship.cpp
#include "Starship.hpp"

Starship::Starship() {}
Starship::~Starship() {}

void Starship::buildGeometry(std::vector<StarshipVertex>& starshipVertexBuffer)
{
    // Just render a basic square for now;

    float x1 = 50.0f;
    float y1 = 50.0f;
    float x2 = x1 + 20.0f;
    float y2 = y1 + 20.0f;

    float u = 0.0f;
    float v = 0.0f;
    float r = 1.0f;
    float g = 1.0f;
    float b = 1.0f;
    float a = 1.0f;

    starshipVertexBuffer.push_back({x1, y1, u, v, r, g, b, a});
    starshipVertexBuffer.push_back({x2, y1, u, v, r, g, b, a});
    starshipVertexBuffer.push_back({x1, y2, u, v, r, g, b, a});
    starshipVertexBuffer.push_back({x2, y1, u, v, r, g, b, a});
    starshipVertexBuffer.push_back({x2, y2, u, v, r, g, b, a});
    starshipVertexBuffer.push_back({x1, y2, u, v, r, g, b, a});
}