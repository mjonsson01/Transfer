// File: Transfer/src/Player/Starship.hpp

#pragma once

// Standard includes
#include <iostream>
#include <vector>

// Custom includes
#include "Core/InputState.hpp"
#include "Core/UIState.hpp"
#include "Utilities/Constants/PhysicsConstants.hpp"
#include "Utilities/Math/Vector2D.hpp"
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
    Vector2D getPointingVector();

  private:
    Vector2D velocity;
    Vector2D unitPointingVector;
    double rotation;
    Vector2D position;
    Vector2D prevPosition;
    float shipSize; // square magnitude
};