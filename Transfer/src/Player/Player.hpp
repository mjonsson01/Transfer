// File: Transfer/src/Player/Player.hpp
#pragma once

// Custom includes
#include "Player/Starship.hpp"

class Player
{
  public:
    Player();
    ~Player();

    Starship starship;
};