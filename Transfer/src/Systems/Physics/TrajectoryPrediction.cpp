// File: Transfer/src/Systems/Physics/TrajectoryPrediction.cpp

// PhysicsSystem's trajectory prediction (planning mode). These are PhysicsSystem member functions, declared in
// PhysicsSystem.hpp; they live in their own file to keep PhysicsSystem.cpp readable.
//
// The prediction is a COPY of the universe, stepped with integrateUniversePhysics: exactly the forces and motion of
// the real tick, without collisions. So a line is exactly what the game will do, as long as nothing collides.

#include "Systems/PhysicsSystem.hpp"

// Custom Imports
#include "Utilities/Constants/PhysicsConstants.hpp"

// Standard Library Imports
#include <chrono>
#include <cmath>

// How many prediction ticks reach the horizon (7200 at 60 s and 120 ticks per second). Rounded, not truncated:
// in float, 60 / (1 / 120) is 7199.9994, which a plain cast would cut down to 7199
static const int TRAJECTORY_HORIZON_TICKS = static_cast<int>(std::lround(TRAJECTORY_HORIZON_SECONDS * INV_PHYSICS_TIME_STEP));

// Exact comparison on purpose: while the game is paused for planning, nothing changes these numbers except an edit
static bool sameVector(const DynamoEngine::Vector2D& a, const DynamoEngine::Vector2D& b)
{
    return a.x_val == b.x_val && a.y_val == b.y_val;
}

void PhysicsSystem::updateTrajectories(float frame_seconds, const UniverseState& universe_state,
                                       Trajectories& trajectories, UIState& ui_state)
{
    DEPRECATED_InputState& input_state = ui_state.getMutableDEPRECATED_InputState();

    // Leaving planning mode makes the lines vanish: once time runs, they no longer describe the future
    if (!ui_state.getIsPlanning())
    {
        if (trajectories.status != TraceStatus::Idle)
        {
            clearTrajectories(trajectories);
        }
        input_state.toggleTrajectoriesRequested = false;
        return;
    }

    // I: start a trace from the universe as it is now, or clear the lines that are there
    if (input_state.toggleTrajectoriesRequested)
    {
        input_state.toggleTrajectoriesRequested = false;
        if (trajectories.status == TraceStatus::Idle)
        {
            startTrajectories(universe_state, trajectories);
        }
        else
        {
            clearTrajectories(trajectories);
            return;
        }
    }

    if (trajectories.status == TraceStatus::Idle)
    {
        return; // no lines: waiting for I
    }

    // Any edit (a spawn, a delete, ...) makes the lines wrong: wipe them and wait for I again
    if (universeChangedSinceTraceStart(universe_state))
    {
        clearTrajectories(trajectories);
        return;
    }

    trajectories.clock += frame_seconds; // keeps running after Finished, so the newest points still fade in
    advanceIntegratedTrajectories(frame_seconds, trajectories);
}

void PhysicsSystem::startTrajectories(const UniverseState& universe_state, Trajectories& trajectories)
{
    m_prediction_universe_state = universe_state;
    m_prediction_universe_state.ship.cutThrust();
    m_trace_start_universe_state = universe_state;
    m_trace_time_accumulator = 0.0f;
    m_trace_ticks_done = 0;

    trajectories.status = TraceStatus::Tracing;
    trajectories.clock = 0.0f;
    trajectories.macro_body_paths.assign(universe_state.macro_bodies.size(), TracedPath());
    trajectories.ship_path.clear();
    recordTrajectoryPoints(trajectories); // every line starts where its body is now
}

void PhysicsSystem::clearTrajectories(Trajectories& trajectories)
{
    trajectories.status = TraceStatus::Idle;
    trajectories.clock = 0.0f;
    trajectories.macro_body_paths.clear();
    trajectories.ship_path.clear();

    // Hand back the copies' memory (up to ~3 MB each with a full particle budget)
    m_prediction_universe_state = UniverseState();
    m_trace_start_universe_state = UniverseState();
}

bool PhysicsSystem::universeChangedSinceTraceStart(const UniverseState& universe_state) const
{
    const UniverseState& start = m_trace_start_universe_state;

    if (universe_state.macro_bodies.size() != start.macro_bodies.size() ||
        universe_state.particles.size() != start.particles.size())
    {
        return true;
    }
    for (size_t i = 0; i < universe_state.macro_bodies.size(); ++i)
    {
        const GravitationalBody& now = universe_state.macro_bodies[i];
        const GravitationalBody& then = start.macro_bodies[i];
        if (!sameVector(now.position, then.position) || !sameVector(now.velocity, then.velocity) ||
            now.mass != then.mass || now.radius != then.radius || now.isForceStatic != then.isForceStatic)
        {
            return true;
        }
    }
    return !sameVector(universe_state.ship.center(), start.ship.center()) ||
           !sameVector(universe_state.ship.velocity(), start.ship.velocity());
}

void PhysicsSystem::advanceIntegratedTrajectories(float frame_seconds, Trajectories& trajectories)
{
    if (trajectories.status != TraceStatus::Tracing)
    {
        return;
    }

    // Same accumulator idea as Game::Run, at TRAJECTORY_TRACE_SPEED, capped so one slow frame can't snowball
    m_trace_time_accumulator += frame_seconds * TRAJECTORY_TRACE_SPEED;
    const float max_backlog = MAX_TRACE_TICKS_PER_FRAME * PHYSICS_TIME_STEP;
    if (m_trace_time_accumulator > max_backlog)
    {
        m_trace_time_accumulator = max_backlog;
    }

    // The time budget keeps the frame rate smooth in heavy scenes (thousands of debris particles make every tick
    // expensive): when this frame's share is used up, the rest waits for the next frame and the lines grow slower
    const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
    while (m_trace_time_accumulator >= PHYSICS_TIME_STEP && m_trace_ticks_done < TRAJECTORY_HORIZON_TICKS)
    {
        integrateUniversePhysics(m_prediction_universe_state); // the integrator only: no collisions, no controls
        m_trace_time_accumulator -= PHYSICS_TIME_STEP;
        m_trace_ticks_done += 1;
        if (m_trace_ticks_done % TRAJECTORY_TICKS_PER_POINT == 0)
        {
            recordTrajectoryPoints(trajectories);
        }

        std::chrono::duration<double, std::milli> spent = std::chrono::steady_clock::now() - started;
        if (spent.count() >= MAX_TRACE_MILLISECONDS_PER_FRAME)
        {
            break;
        }
    }

    if (m_trace_ticks_done >= TRAJECTORY_HORIZON_TICKS)
    {
        trajectories.status = TraceStatus::Finished;
    }
}

void PhysicsSystem::recordTrajectoryPoints(Trajectories& trajectories)
{
    // Path i belongs to macro body i: the prediction never collides or cleans up, so no body moves in the list
    const std::vector<GravitationalBody>& macro_bodies = m_prediction_universe_state.macro_bodies;
    for (size_t i = 0; i < macro_bodies.size(); ++i)
    {
        trajectories.macro_body_paths[i].push_back({macro_bodies[i].position, trajectories.clock});
    }
    trajectories.ship_path.push_back({m_prediction_universe_state.ship.center(), trajectories.clock});
}
