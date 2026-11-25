// File: Transfer/src/Systems/Physics/Trajectories.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/Math/Vector2.hpp"

// Standard Library Imports
#include <cstdint>
#include <vector>

// One recorded point of a predicted path
struct TracedPoint
{
    DynamoEngine::Vector2D position;
    float drawn_at = 0.0f; // the trace clock when this point was recorded: its fade-in starts here
};

// Everything predicted for one body so far, oldest point first
using TracedPath = std::vector<TracedPoint>;

enum class TraceStatus
{
    Idle,     // no lines: waiting for I
    Tracing,  // the lines are growing
    Finished, // the lines reached the horizon
};

// The predicted paths planning mode draws. PhysicsSystem writes them; RenderSystem and the readout only read them.
struct Trajectories
{
    TraceStatus status = TraceStatus::Idle;
    std::vector<TracedPath> macro_body_paths; // path i = macro body i as it was when the trace started
    TracedPath ship_path;
    // While placing a planet in planning mode: where it would go (empty otherwise). Live, so it doesn't fade in.
    TracedPath spawn_preview_path;
    float clock = 0.0f; // real seconds since the trace started (drives the fade-in and the bright tip)
    float predicted_seconds = 0.0f;         // how far ahead the traced lines reach so far, in GAME seconds
    float preview_predicted_seconds = 0.0f; // the same for the placement preview (0 when not placing)
    uint32_t trace_id = 0; // +1 every time a trace starts: tells the renderer its uploaded lines are from an old trace
};