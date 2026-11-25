// File: Transfer/src/Utilities/Constants/GameSystemConstants.hpp

#pragma once

// Standard Library Imports
#include <cstddef>
// #ifndef GAME_SYSTEM_CONSTANTS_H
// #define GAME_SYSTEM_CONSTANTS_H

// FPS target for system
// constexpr int TARGET_FPS = 165;
constexpr int TARGET_FPS = 60;

// Target frame delay in ms
const double FRAME_DELAY_MS = 1000.0 / TARGET_FPS;

// Delay to update the text of the rolling average frame counter UI Element
constexpr int FPS_UPDATE_DELTA_MS = 250;

// Time scaling factors to control rendering speed without changing the Physics
// System behaviors
const double MIN_TIME_SCALE_FACTOR = 0.1;     // negative is allowed and runs only due to alpha interpolation. Prevent 0
                                              // speed on slider since we now have Planning mode
const double REGULAR_TIME_SCALE_FACTOR = 1.0; // 1x speed
const double MAX_TIME_SCALE_FACTOR = 3.0;     // 3x speed

// Vertical and horizontal resolutions (starting resolution)

// can scale up as needed

// constexpr int SCREEN_HEIGHT = 1080;
// constexpr int SCREEN_WIDTH = 1920;
constexpr int SCREEN_HEIGHT = 720;
constexpr int SCREEN_WIDTH = 1280;

// Background Star count
constexpr int STAR_NUM = 100000;

constexpr bool VIEW_DEBUG = false;
// The F4 profiler reports averages and 1% lows over this many frames (5 s at 60 FPS)
constexpr size_t PROFILER_WINDOW_FRAMES = 300;
// Most lines the F4 profiler overlay shows (a header plus one line per section)
constexpr size_t MAX_PROFILER_LINES = 12;

const double MIN_ZOOM = 0.01;

const double MAX_ZOOM = 10.0;

const double STARTUP_ZOOM_VALUE = 0.35;

const double STAR_PARALLAX_FACTOR = 0.5;

constexpr int MAX_VELOCITY_VECTOR_VERTICES = 9; // vertices on a previewbody vector
