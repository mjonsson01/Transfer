// File: Transfer/src/Core/Game.h

#pragma once

// Standard Library Imports
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>

// Custom Imports
#include "Core/GameState.hpp"
#include "Core/PhysicsThread.hpp"
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
    // Tells the physics thread this frame's controls, speed and planning state, and hands it the queued changes:
    // spawns / undo / clear from the player, and edits to the universe from the scenes
    void sendRequestsToPhysics();
    // Swaps the physics thread's newest universe into game_state (if there is one) and picks up its numbers: ticks
    // run, whether it fell behind, the effective speed, the F5 report
    void takeUniverseFromPhysics();
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
    // The simulation, on its own thread (it owns the live universe and its own PhysicsSystem)
    PhysicsThread m_physics_thread;
    // A second PhysicsSystem, used only on the main thread for planning mode's trajectory lines and the placement
    // preview (they step copies of the universe), so the two threads never share one
    PhysicsSystem predictionPhysics;
    RenderSystem renderSystem;   // Manages all rendering operations
    AudioSystem audioSystem;     // Manages all audio operations

    // Every scene (and its UI), and which one is active
    DynamoEngine::SceneManager scenes;

    // F4 profiler: how long each part of the frame takes. The section numbers come from m_profiler.addSection.
    DynamoEngine::FrameProfiler m_profiler;
    struct ProfilerSections
    {
        int input = 0;
        int to_physics = 0; // handing this frame's requests and changes to the physics thread
        int audio = 0;
        int from_physics = 0; // taking the physics thread's newest universe (a swap: the ticks run on their own thread)
        // The physics thread's own work per frame (ticks, commands, the copy it publishes), as reported with each
        // universe. It runs AT THE SAME TIME as the main thread, so it isn't part of frame_work below.
        int physics_thread = 0;
        int trajectories = 0; // planning mode's prediction
        int render = 0;       // CPU side only: building and submitting the frame (the GPU works on it afterwards)
        int frame_work = 0;   // everything the MAIN thread did this frame, before the frame limiter waits
        int frame_total = 0;  // the whole frame, including the wait (about 16.7 ms at 60 FPS)
    };
    ProfilerSections m_profile;
    // Physics ticks per frame, for the profiler overlays: counted over the F4 window, then turned into m_ticks_line
    int m_ticks_in_window = 0;
    int m_most_ticks_in_a_frame = 0;
    int m_frames_behind_in_window = 0; // frames for which the physics thread dropped owed game time
    std::string m_ticks_line = "Ticks/frame: -";
    // When the universe in game_state was published (for interpolating between its last two positions)
    Uint64 m_universe_publish_time = 0;
    // The physics thread's latest F5 report (it arrives with a universe)
    std::vector<DynamoEngine::FrameProfiler::SectionReport> m_physics_report;
};