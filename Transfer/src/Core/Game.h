// File: Transfer/src/Core/Game.h

#pragma once

// Custom Imports
#include "Core/GameState.h"
#include "Core/InputState.h"
#include "Core/UIState.h"
#include "Systems/AudioSystem.h"
#include "Systems/InputSystem.h"
#include "Systems/PhysicsSystem.h"
#include "Systems/RenderSystem.h"
#include "Systems/UISystem.h"
#include "Utilities/Constants/EngineConstants.h"
#include "Utilities/Constants/GameSystemConstants.h"

// Standard Library Imports
#include <numeric>
#include <unordered_map>

class Game
{
  public:
    // Constructor and Destructor
    Game();
    ~Game();

    // Initializes SDL windows, renderer, and starts the main game loop by
    // calling Run()
    void StartGame();
    // Tears down the 'systems' and cleans up allocated resources.
    void EndGame();
    // Game Loop Entry Point
    void Run();

  public:
  private:
    // Core Game Loop Methods
    void ProcessInput();          // Handles User Input from keyboard and mouse events
    void IntegratePhysicsFrame(); // Integrate Game State Bodies forward by 1 time step
    void UpdateInstantiations();  // Updates Game State Bodies from User Input
    void RenderFrame();           // Renders the current frame to the screen including UI
    void PlayAudio();             // Plays audio based on the current game state and UI state
    // Helpers for Run() method
    void updateFPS(Uint32 renderEnd, Uint32 lastRender, float& fpsAccumulator,
                   float& currentFPS); // Rolling average frame calculation and update
    void limitFrameRate(Uint32 renderStart,
                        Uint32 renderEnd); // Soft Clamps frame rate based on
                                           // static Game System Constant value

  private:
    // Systems and State
    GameState gameState;         // Contains all game entities and their states
    UIState UIState;             // Contains all UI related states
    InputSystem inputSystem;     // Manages all user input
    PhysicsSystem physicsSystem; // Manages physics calculations and Frame Updates
    RenderSystem renderSystem;   // Manages all rendering operations
    AudioSystem audioSystem;     // Manages all audio operations
    UISystem UISystem;
};