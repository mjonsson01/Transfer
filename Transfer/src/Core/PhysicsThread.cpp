// File: Transfer/src/Core/PhysicsThread.cpp
#include "Core/PhysicsThread.hpp"

// Custom Imports
#include "Constants/GameSystemConstants.hpp"
#include "Systems/Physics/PhysicsConstants.hpp"

// SDL Imports
#include <SDL3/SDL_timer.h>

// Standard Library Imports
#include <chrono>
#include <utility>

PhysicsThread::~PhysicsThread() { stop(); }

void PhysicsThread::start(const UniverseState& starting_universe)
{
    m_universe = starting_universe;
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_published = starting_universe;
        m_publish_unread = true;
        m_published_handoff = Handoff();
        m_published_handoff.publish_time = SDL_GetPerformanceCounter();
        m_stop_requested = false;
    }
    // std::thread starts running `run` on `this` object right away, on a new thread
    m_thread = std::thread(&PhysicsThread::run, this);
}

void PhysicsThread::stop()
{
    if (!m_thread.joinable())
    {
        return; // never started, or already stopped
    }
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_stop_requested = true;
        m_has_news = true;
    }
    m_wake.notify_one();
    m_thread.join(); // waits until run() has returned
}

void PhysicsThread::sendFrameRequest(const FrameRequest& request)
{
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_request = request;
        m_has_news = true;
    }
    m_wake.notify_one(); // a new speed or new controls take effect right away, not after the current wait
}

void PhysicsThread::queueSpawnCommands(const PlayerInput& input, const SpawnSettings& settings)
{
    Command command;
    command.is_edit = false;
    command.spawn_input = input;
    command.spawn_settings = settings;
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_pending_commands.push_back(command);
        m_has_news = true;
    }
    m_wake.notify_one();
}

void PhysicsThread::queueUniverseEdit(std::function<void(UniverseState&)> edit)
{
    Command command;
    command.is_edit = true;
    command.edit = std::move(edit);
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_pending_commands.push_back(std::move(command));
        m_has_news = true;
    }
    m_wake.notify_one();
}

bool PhysicsThread::takeNewestUniverse(UniverseState& presented, Handoff& handoff_out)
{
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (!m_publish_unread)
        {
            return false;
        }
        // Swap, don't copy: `presented` gets the published universe, and the main thread's old one becomes the buffer
        // the physics thread publishes into next time (std::swap exchanges the vectors' insides: nothing is copied)
        std::swap(presented, m_published);
        handoff_out = std::move(m_published_handoff);
        m_published_handoff = Handoff();
        m_publish_unread = false;
        m_has_news = true;
    }
    m_wake.notify_one(); // an idle physics thread may have a change waiting to be published
    return true;
}

bool PhysicsThread::tryPublish(Uint64 now)
{
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (m_publish_unread)
        {
            return false; // the main thread hasn't taken the previous one yet: publish later
        }
    }
    // The copy happens here, outside the lock (m_staging belongs to this thread), so the main thread never waits for it
    m_staging = m_universe;
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        std::swap(m_staging, m_published);
        m_publish_unread = true;
        m_published_handoff.ticks = m_ticks_since_take;
        m_published_handoff.dropped_time = m_dropped_time_since_take;
        // Behind: what was really achieved. Keeping up: exactly the slider, right away (so moving the slider doesn't
        // make the readout show the old speed until the next window ends)
        m_published_handoff.effective_speed = m_was_behind ? m_measured_speed : m_requested_speed;
        m_published_handoff.publish_time = now;
        m_published_handoff.busy_ms = m_busy_ms_since_take;
        m_published_handoff.has_new_profiler_report = m_has_new_profiler_report;
        if (m_has_new_profiler_report)
        {
            m_published_handoff.profiler_report = m_profiler_report;
        }
    }
    m_ticks_since_take = 0;
    m_busy_ms_since_take = 0.0;
    m_dropped_time_since_take = false;
    m_has_new_profiler_report = false;
    return true;
}

void PhysicsThread::run()
{
    const double clock_frequency = static_cast<double>(SDL_GetPerformanceFrequency());
    Uint64 last_clock = SDL_GetPerformanceCounter();
    double owed_game_seconds = 0.0; // game time the slider asked for that hasn't been simulated yet
    bool needs_publish = false;     // the live universe changed since the last publish

    while (true)
    {
        const Uint64 pass_start = SDL_GetPerformanceCounter(); // the work of this pass is timed, the wait isn't

        // 1. What the main thread sent (copied out under the lock, quickly)
        FrameRequest request;
        std::vector<Command> commands;
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            if (m_stop_requested)
            {
                return;
            }
            request = m_request;
            commands.swap(m_pending_commands);
        }
        m_requested_speed = request.speed;

        // 2. Apply the commands, in the order they were sent, between two ticks
        for (Command& command : commands)
        {
            if (command.is_edit)
            {
                command.edit(m_universe);
            }
            else
            {
                m_physics.updateGravBodyInstantiations(m_universe, command.spawn_input, command.spawn_settings);
            }
            needs_publish = true;
        }

        // 3. Time: how much game time the slider asks for since the last pass
        Uint64 now = SDL_GetPerformanceCounter();
        double real_seconds = static_cast<double>(now - last_clock) / clock_frequency;
        last_clock = now;
        const bool ticking = request.runs_simulation && !request.is_planning;
        if (ticking)
        {
            owed_game_seconds += real_seconds * request.speed;
            // Can't keep up: drop what's owed beyond the backlog cap (the simulation runs slower than requested)
            const double most_owed = MAX_PHYSICS_BACKLOG_TICKS * static_cast<double>(PHYSICS_TIME_STEP);
            if (owed_game_seconds > most_owed)
            {
                owed_game_seconds = most_owed;
                m_dropped_time_since_take = true;
                m_window_dropped_time = true;
            }
        }
        else
        {
            owed_game_seconds = 0.0;
            // Planning mode stops time, but W/S still set up the burn, in real time
            if (request.runs_simulation && request.is_planning && request.ship_controls.thrust != 0)
            {
                m_physics.adjustThrottleWhilePlanning(m_universe, request.ship_controls,
                                                      static_cast<float>(real_seconds));
                needs_publish = true;
            }
        }

        // 4. At most one tick per pass, so new controls and commands are picked up between every two ticks
        if (owed_game_seconds >= PHYSICS_TIME_STEP)
        {
            m_physics.advanceUniverseState(m_universe, request.ship_controls);
            owed_game_seconds -= PHYSICS_TIME_STEP;
            m_ticks_since_take += 1;
            m_window_game_seconds += PHYSICS_TIME_STEP;
            needs_publish = true;
            if (m_physics.profiler().takeNewReport())
            {
                m_profiler_report = m_physics.profiler().report();
                m_has_new_profiler_report = true;
            }
        }

        // 5. The effective speed, measured every EFFECTIVE_SPEED_WINDOW_SECONDS: game time simulated over real time
        //    passed. It only matters if game time was dropped in the window (m_was_behind); otherwise the simulation
        //    kept up and the effective speed is simply the slider (tryPublish)
        m_window_real_seconds += real_seconds;
        if (m_window_real_seconds >= EFFECTIVE_SPEED_WINDOW_SECONDS)
        {
            m_measured_speed = m_window_game_seconds / m_window_real_seconds;
            m_was_behind = m_window_dropped_time;
            m_window_game_seconds = 0.0;
            m_window_real_seconds = 0.0;
            m_window_dropped_time = false;
        }

        // 6. Publish (if the main thread took the previous one)
        if (needs_publish && tryPublish(now))
        {
            needs_publish = false;
        }

        // (This pass's work goes to the next publish: the copy above is counted there too)
        m_busy_ms_since_take +=
            static_cast<double>(SDL_GetPerformanceCounter() - pass_start) * 1000.0 / clock_frequency;

        // 7. Nothing owed right now: sleep until the next tick is due, or until the main thread sends something
        if (owed_game_seconds < PHYSICS_TIME_STEP)
        {
            double wait_seconds = 0.1; // idle (time stopped, or a menu): the main thread wakes us when needed
            if (ticking && request.speed > 0.0)
            {
                wait_seconds = (PHYSICS_TIME_STEP - owed_game_seconds) / request.speed;
            }
            else if (request.is_planning && request.ship_controls.thrust != 0)
            {
                wait_seconds = static_cast<double>(PHYSICS_TIME_STEP); // keep moving the throttle smoothly
            }
            std::unique_lock<std::mutex> lock(m_mutex);
            // Returns when the time is up, or earlier if m_has_news was set (checked first, so a notify that came in
            // before this line isn't missed)
            m_wake.wait_for(lock, std::chrono::duration<double>(wait_seconds), [this]() { return m_has_news; });
            m_has_news = false;
        }
    }
}
