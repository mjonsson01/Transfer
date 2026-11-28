// File: Transfer/src/Core/PhysicsThread.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/Debug/FrameProfiler.hpp"
#include "Systems/Input/PlayerInput.hpp"
#include "Systems/Physics/PhysicsSystem.hpp"
#include "Systems/Physics/SpawnSettings.hpp"
#include "Systems/Physics/UniverseState.hpp"

// SDL Imports
#include <SDL3/SDL_stdinc.h>

// Standard Library Imports
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

// The simulation, running on its own thread (Notes/PHYSICS_REQUIREMENTS.md, requirement C3).
//
// The physics thread OWNS the live universe and its PhysicsSystem, and runs ticks on its own clock: as many as the
// speed slider asks for, back to back when it's behind. The main thread never touches the live universe. Instead:
//   - it SENDS: the ship's controls, the speed, planning on/off every frame (sendFrameRequest), and changes to the
//     universe as commands (spawns, undo, clear, scene edits), which the physics thread applies between two ticks
//   - it TAKES a copy: after a tick, the physics thread copies the universe into a spare buffer and "publishes" it;
//     each frame the main thread swaps the newest published copy into the universe it renders (takeNewestUniverse)
// So the threads only meet inside short locked blocks that swap buffers (vectors swap their insides: no copying
// under the lock). Rendering, input and the trajectory lines all read the main thread's copy.
//
// The ticks are the same ticks the game always ran, in the same order, so the results are identical (the
// determinism check runs PhysicsSystem directly and doesn't need this class).
//
// When a tick takes longer than the game time it covers, the owed game time grows until MAX_PHYSICS_BACKLOG_TICKS,
// and the rest is dropped: the simulation runs slower than requested ("effective speed"), but the main thread keeps
// rendering at full frame rate, because it never waits for physics.
class PhysicsThread
{
  public:
    // What the main thread tells the simulation every frame
    struct FrameRequest
    {
        ShipControls ship_controls;  // the keys held right now (used by every tick until the next request)
        double speed = 1.0;          // the speed slider
        bool is_planning = false;    // planning mode stops time (W/S still move the throttle)
        bool runs_simulation = true; // false in scenes without a simulation (menus): no ticks
    };

    // What came back with the newest published universe
    struct Handoff
    {
        int ticks = 0;                   // ticks run since the main thread last took a universe
        bool dropped_time = false;       // whether owed game time was dropped since then (couldn't keep up)
        double effective_speed = 1.0;    // game seconds per real second: measured when behind, else the slider
        Uint64 publish_time = 0;         // SDL_GetPerformanceCounter() when it was published
        double busy_ms = 0.0;            // how long the physics thread worked (ticks, commands, copies) for this one
        bool has_new_profiler_report = false;
        std::vector<DynamoEngine::FrameProfiler::SectionReport> profiler_report; // F5, when has_new_profiler_report
    };

    PhysicsThread() = default;
    ~PhysicsThread(); // stops the thread if it's still running

    // One thread, one owner: copying would make two threads share one universe
    PhysicsThread(const PhysicsThread&) = delete;
    PhysicsThread& operator=(const PhysicsThread&) = delete;

    // Starts the thread with this universe (copied), and publishes it right away so the first frame has one
    void start(const UniverseState& starting_universe);
    // Asks the thread to finish its current tick and waits for it to end
    void stop();

    // The latest controls, speed and planning state (replaces the previous request)
    void sendFrameRequest(const FrameRequest& request);
    // A spawn, undo or clear (the commands in `input`), with the spawn panel's settings at that moment
    void queueSpawnCommands(const PlayerInput& input, const SpawnSettings& settings);
    // Any other change to the universe (e.g. the test scene adding its collision box), applied between two ticks
    void queueUniverseEdit(std::function<void(UniverseState&)> edit);

    // If a newer universe was published since the last call: swaps it into `presented` (which the physics thread then
    // reuses as a spare buffer), fills `handoff_out` and returns true. Otherwise leaves both alone and returns false.
    bool takeNewestUniverse(UniverseState& presented, Handoff& handoff_out);

  private:
    // A change the main thread asked for: either spawn commands (with their settings) or an edit function
    struct Command
    {
        bool is_edit = false;
        PlayerInput spawn_input;
        SpawnSettings spawn_settings;
        std::function<void(UniverseState&)> edit;
    };

    // The thread's whole life: read requests, apply commands, run owed ticks, publish, wait
    void run();
    // Copies the live universe into the spare buffer and swaps it in as the newest publish. Only when the main thread
    // has taken the previous one (at most one copy per frame); returns false if it didn't publish.
    bool tryPublish(Uint64 now);

    // --- Only the physics thread touches these ---
    PhysicsSystem m_physics;
    UniverseState m_universe; // the live universe
    UniverseState m_staging;  // the spare buffer a publish is copied into
    int m_ticks_since_take = 0;
    double m_busy_ms_since_take = 0.0; // time spent working (not waiting) since the last publish
    bool m_dropped_time_since_take = false;
    bool m_has_new_profiler_report = false;
    std::vector<DynamoEngine::FrameProfiler::SectionReport> m_profiler_report;
    // The effective speed (see run()): measured over windows of EFFECTIVE_SPEED_WINDOW_SECONDS
    double m_requested_speed = 1.0; // the slider, as of the latest request
    double m_measured_speed = 1.0;  // game seconds per real second over the last finished window
    bool m_was_behind = false;      // whether game time was dropped in the last finished window
    double m_window_game_seconds = 0.0;
    double m_window_real_seconds = 0.0;
    bool m_window_dropped_time = false;

    // --- Shared: only touched while m_mutex is locked ---
    std::mutex m_mutex;
    std::condition_variable m_wake; // wakes the physics thread early: new commands, a taken publish, stop
    bool m_has_news = false;        // set with m_wake, so a wake-up that happens before the wait isn't lost
    bool m_stop_requested = false;
    FrameRequest m_request;
    std::vector<Command> m_pending_commands;
    UniverseState m_published; // the newest published universe (or the main thread's old buffer, once taken)
    bool m_publish_unread = false;
    Handoff m_published_handoff;

    std::thread m_thread;
};
