// File: Transfer/src/Core/Game.h

#pragma once

// Custom Imports
#include "Core/GameState.hpp"
#include "Core/InputState.hpp"
#include "Core/UIState.hpp"
#include "Entities/Sound/MusicModeEnum.hpp"
#include "Systems/AudioSystem.hpp"
#include "Systems/InputSystem.hpp"
#include "Systems/PhysicsSystem.hpp"
#include "Systems/RenderSystem.hpp"
#include "Systems/UISystem.hpp"
#include "Utilities/Constants/EngineConstants.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"

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
    // void updateFPS(Uint32 renderEnd, Uint32 lastRender, float& fpsAccumulator,
    //                float& currentFPS); // Rolling average frame calculation and update
    void updateFPS(Uint64 renderEnd, Uint64 lastRender, float& fpsAccumulator, float& currentFPS);
    // void limitFrameRate(Uint32 renderStart,
    //                     Uint32 renderEnd); // Soft Clamps frame rate based on
    //                                        // static Game System Constant value
    void limitFrameRate(Uint64 renderStart, Uint64 renderEnd, Uint64 perfFreq);

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