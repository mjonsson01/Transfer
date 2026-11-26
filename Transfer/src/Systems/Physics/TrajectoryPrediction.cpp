// File: Transfer/src/Systems/Physics/TrajectoryPrediction.cpp

// PhysicsSystem's trajectory prediction (planning mode). These are PhysicsSystem member functions, declared in
// PhysicsSystem.hpp; they live in their own file to keep PhysicsSystem.cpp readable.
//
// The prediction is a COPY of the universe, stepped with integrateUniversePhysics: exactly the forces and motion of
// the real tick, without collisions. So a line is exactly what the game will do, as long as nothing collides.

#include "Systems/Physics/PhysicsSystem.hpp"

// Custom Imports
#include "Systems/Physics/PhysicsConstants.hpp"

// Standard Library Imports
#include <chrono>
#include <cmath>

// How many prediction ticks reach the horizon (21600 at 180 s and 120 ticks per second). Rounded, not truncated:
// in float, seconds / (1 / 120) can land just under the whole number (60 s gave 7199.9994), and a cast cuts it down
static const int TRAJECTORY_HORIZON_TICKS = static_cast<int>(std::lround(TRAJECTORY_HORIZON_SECONDS * INV_PHYSICS_TIME_STEP));

// Exact comparison on purpose: while the game is paused for planning, nothing changes these numbers except an edit
static bool sameVector(const DynamoEngine::Vector2D& a, const DynamoEngine::Vector2D& b)
{
    return a.x_val == b.x_val && a.y_val == b.y_val;
}

void PhysicsSystem::updateTrajectories(float frame_seconds, const UniverseState& universe_state,
                                       const PlayerInput& player_input, Trajectories& trajectories,
                                       const UIState& ui_state)
{
    // Leaving planning mode makes the lines vanish: once time runs, they no longer describe the future
    if (!ui_state.getIsPlanning())
    {
        if (trajectories.status != TraceStatus::Idle)
        {
            clearTrajectories(trajectories);
        }
        clearSpawnPreviewTrajectory(trajectories);
        return;
    }

    // The placement preview runs on its own, with or without the traced lines
    updateSpawnPreviewTrajectory(universe_state, player_input, ui_state.spawnSettings(), trajectories);

    // I: start a trace from the universe as it is now, or clear the lines that are there
    if (player_input.commands.toggle_trajectories)
    {
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
    // The ship keeps its current throttle: the line shows where it goes if you hold this burn
    m_prediction_universe_state = universe_state;
    m_trace_start_universe_state = universe_state;
    m_trace_time_accumulator = 0.0f;
    m_trace_ticks_done = 0;

    trajectories.status = TraceStatus::Tracing;
    trajectories.clock = 0.0f;
    trajectories.predicted_seconds = 0.0f;
    trajectories.trace_id += 1;
    trajectories.macro_body_paths.assign(universe_state.macro_bodies.size(), TracedPath());
    trajectories.ship_path.clear();
    recordTrajectoryPoints(trajectories); // every line starts where its body is now
}

void PhysicsSystem::clearTrajectories(Trajectories& trajectories)
{
    trajectories.status = TraceStatus::Idle;
    trajectories.clock = 0.0f;
    trajectories.predicted_seconds = 0.0f;
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
    // The throttle counts too: it can be changed while planning, and the prediction keeps it
    return !sameVector(universe_state.ship.center(), start.ship.center()) ||
           !sameVector(universe_state.ship.velocity(), start.ship.velocity()) ||
           universe_state.ship.throttle() != start.ship.throttle();
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

    trajectories.predicted_seconds = m_trace_ticks_done * PHYSICS_TIME_STEP;
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
    if (m_prediction_universe_state.has_ship) // no ship, no ship line
    {
        trajectories.ship_path.push_back({m_prediction_universe_state.ship.center(), trajectories.clock});
    }
}

void PhysicsSystem::updateSpawnPreviewTrajectory(const UniverseState& universe_state, const PlayerInput& player_input,
                                                 const SpawnSettings& spawn_settings, Trajectories& trajectories)
{
    GravitationalBody body;
    const bool is_aiming_planet = player_input.is_aiming_spawn && !player_input.spawn_aim.is_cluster;
    if (!is_aiming_planet || !makeSpawnPreviewBody(player_input.spawn_aim, spawn_settings, body))
    {
        clearSpawnPreviewTrajectory(trajectories); // not placing a planet (or none would spawn)
        return;
    }

    // A different planet than last frame (moved, new velocity, other size or mass...): start again from the beginning
    const GravitationalBody& last = m_preview_spawn_body;
    // (or the throttle was changed while planning: the preview's ship keeps the throttle it started with)
    bool changed = !m_has_spawn_preview || !sameVector(body.position, last.position) ||
                   !sameVector(body.velocity, last.velocity) || body.mass != last.mass || body.radius != last.radius ||
                   body.isForceStatic != last.isForceStatic ||
                   universe_state.ship.throttle() != m_preview_universe_state.ship.throttle();
    if (changed)
    {
        // The planets and the ship, plus the planet being placed. NO debris, on purpose: the preview restarts every
        // time the mouse moves, and debris would make every restart slow (thousands of particles to copy and pull).
        // So it ignores particle clouds; the traced lines (I) include them.
        m_preview_universe_state.macro_bodies = universe_state.macro_bodies;
        m_preview_universe_state.macro_bodies.push_back(body);
        m_preview_universe_state.particles.clear();
        m_preview_universe_state.ship = universe_state.ship;
        m_preview_universe_state.has_ship = universe_state.has_ship;
        m_preview_universe_state.world_min_corner = universe_state.world_min_corner; // so the preview bounces too
        m_preview_universe_state.world_max_corner = universe_state.world_max_corner;
        m_preview_universe_state.reflective_walls = universe_state.reflective_walls;
        m_preview_universe_state.has_collision_box = universe_state.has_collision_box;
        m_preview_universe_state.collision_box_min_corner = universe_state.collision_box_min_corner;
        m_preview_universe_state.collision_box_max_corner = universe_state.collision_box_max_corner;
        m_preview_spawn_body = body;
        m_preview_ticks_done = 0;
        m_has_spawn_preview = true;
        trajectories.spawn_preview_path.clear();
        trajectories.spawn_preview_path.push_back({body.position, 0.0f});
    }

    // No speed limit, only the time budget: the line appears as fast as the frame allows
    const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
    while (m_preview_ticks_done < TRAJECTORY_HORIZON_TICKS)
    {
        integrateUniversePhysics(m_preview_universe_state);
        m_preview_ticks_done += 1;
        if (m_preview_ticks_done % TRAJECTORY_TICKS_PER_POINT == 0)
        {
            trajectories.spawn_preview_path.push_back({m_preview_universe_state.macro_bodies.back().position, 0.0f});
        }

        std::chrono::duration<double, std::milli> spent = std::chrono::steady_clock::now() - started;
        if (spent.count() >= MAX_PREVIEW_MILLISECONDS_PER_FRAME)
        {
            break;
        }
    }
    trajectories.preview_predicted_seconds = m_preview_ticks_done * PHYSICS_TIME_STEP;
}

void PhysicsSystem::clearSpawnPreviewTrajectory(Trajectories& trajectories)
{
    if (!m_has_spawn_preview)
    {
        return;
    }
    trajectories.spawn_preview_path.clear();
    trajectories.preview_predicted_seconds = 0.0f;
    m_preview_universe_state = UniverseState();
    m_has_spawn_preview = false;
}
