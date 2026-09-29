// File: Transfer/src/Player/Starship.hpp

#pragma once

// Standard includes
#include <iostream>
#include <vector>

// Custom includes
#include "Core/DEPRECATED_InputState.hpp"
#include "Core/UIState.hpp"
#include "DynamoEngine/Math/Vector2.hpp"
#include "DynamoEngine/Rendering/UIVertex.hpp"
#include "Utilities/Constants/PhysicsConstants.hpp"
#include "Utilities/Rendering/GPUTypes.hpp"

class Starship
{
  public:
    Starship();
    ~Starship();

    // void UpdateStarshipState();
    // void applyDamage();
    // void applyGravity();
    // void updateAttitude();
    // void applyLinearThrust();
    // void applyAngularThrust();
    void integratePosition();

    void applyVelocity(UIState& uiState);
    void applyRotation(UIState& uiState);
    void buildGeometry(std::vector<StarshipVertex>& starshipVertexBuffer);
    DynamoEngine::Vector2D getPointingVector();

  private:
    DynamoEngine::Vector2D velocity = {};
    double rotation = {};
    DynamoEngine::Vector2D position = {};
    DynamoEngine::Vector2D prevPosition = {};
    float shipSize = {};
};