// File: Transfer/src/Core/Game.h

#pragma once

// Standard Library Imports
#include <iostream>
#include <numeric>
#include <unordered_map>

// Custom Imports
#include "Core/GameState.hpp"
#include "Systems/Input/PlayerInput.hpp"
#include "Core/UIState.hpp"
#include "DynamoEngine/Debug/FrameProfiler.hpp"
#include "DynamoEngine/Scenes/SceneManager.hpp"
#include "Systems/Audio/MusicModeEnum.hpp"
#include "Scenes/TransferScenes.hpp"
#include "Systems/Audio/AudioSystem.hpp"
#include "Systems/Input/InputSystem.hpp"
#include "Systems/Physics/PhysicsSystem.hpp"
#include "Systems/Rendering/RenderSystem.hpp"
#include "Constants/EngineConstants.hpp"
#include "Constants/GameSystemConstants.hpp"
// #include "Utilities/SystemPathUtility.hpp"

class Game
{
  public:
    // Constructor and Destructor
    Game();
    ~Game();

    // Initializes SDL windows and starts the main game loop by
    // calling Run()
    void StartGame();
    // Tears down the 'systems' and cleans up allocated resources.
    void EndGame();
    // Game Loop Entry Point
    void Run();

  public:
  private:
    // Core Game Loop Methods
    // Handles User Input from keyboard and mouse events (the current scene's UI gets it first)
    void ProcessInput(float frame_seconds);
    void IntegratePhysicsFrame(); // Integrate Game State Bodies forward by 1 time step
    void UpdateInstantiations();  // Updates Game State Bodies from User Input
    void RenderFrame();           // Renders the current frame to the screen including UI
    void PlayAudio();             // Plays audio based on the current game state and UI state
    // Helpers for Run() method
    void updateFPS(Uint64 render_end, Uint64 last_render, float& fps_accumulator, float& current_fps);
    void limitFrameRate(Uint64 renderStart, Uint64 render_end, Uint64 perfFreq);
    // Turns the profiler's latest report into the F4 overlay's lines (and prints them while the overlay is shown)
    void publishProfilerReport();
    // The same for the physics profiler (F5: each part of the physics tick)
    void publishPhysicsProfilerReport();
    // Sets the universe's world edges to the star field's extent (no walls: the normal game)
    void setUpWorldEdges();

  private:
    // Systems and State
    GameState game_state;        // Contains all game entities and their states
    UIState ui_state;            // Contains all UI related states
    PlayerInput m_player_input;  // what the player is doing this frame: InputSystem writes it, the others read it
    InputSystem inputSystem;     // Manages all user input
    PhysicsSystem physicsSystem; // Manages physics calculations and Frame Updates
    RenderSystem renderSystem;   // Manages all rendering operations
    AudioSystem audioSystem;     // Manages all audio operations

    // Every scene (and its UI), and which one is active
    DynamoEngine::SceneManager scenes;

    // F4 profiler: how long each part of the frame takes. The section numbers come from m_profiler.addSection.
    DynamoEngine::FrameProfiler m_profiler;
    struct ProfilerSections
    {
        int input = 0;
        int spawns = 0;
        int audio = 0;
        int physics = 0;      // every physics tick this frame, summed
        int trajectories = 0; // planning mode's prediction
        int render = 0;       // CPU side only: building and submitting the frame (the GPU works on it afterwards)
        int frame_work = 0;   // everything above, before the frame limiter waits
        int frame_total = 0;  // the whole frame, including the wait (about 16.7 ms at 60 FPS)
    };
    ProfilerSections m_profile;
};