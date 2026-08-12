// File: Transfer/src/Player/Starship.hpp

#pragma once

// Standard includes
#include <vector>

// Custom includes
#include "Utilities/Rendering/GPUTypes.hpp"

class Starship
{
  public:
    Starship();
    ~Starship();

    void UpdateStarshipState();
    void applyDamage();
    void applyGravity();
    void updateAttitude();
    void applyLinearThrust();
    void applyAngularThrust();

    void buildGeometry(std::vector<StarshipVertex>& starshipVertexBuffer);

  private:
};