// File: Transfer/src/Core/GameState.cpp

#include "Core/Game.h"

GameState::GameState()
{
    // Initialize game state variables if needed
    physicsData.grid.resize(physicsData.gridWidth * physicsData.gridHeight);
}

GameState::~GameState()
{
    // Cleanup if necessary
}
