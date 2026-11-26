// File: Transfer/src/Constants/GameSystemConstants.hpp

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
// Whether the player's ship is in the universe when the game starts (off for now, while other systems are ironed out)
constexpr bool RUN_GAME_WITH_PLAYER_SHIP = false;

// --- Test scene tools (TestVisual scene only: entering it applies these, leaving it removes them) ---
// A collision box around the view when the test scene starts: bodies inside stay inside, bodies outside bounce off it
constexpr double TEST_COLLISION_BOX_WIDTH = 2000.0; // world pixels
constexpr double TEST_COLLISION_BOX_HEIGHT = 2000.0;
// Planets, debris and the ship bounce off the WORLD's edges (the star field's extent, red outline) instead of leaving
// the world (WALL_RESTITUTION)
constexpr bool REFLECTIVE_WORLD_WALLS = false;

// The F4 profiler reports averages and 1% lows over this many frames, and refreshes after each window (1 s at 60
// FPS). Under 100 frames, the "1% low" is simply the slowest frame of the window.
constexpr size_t PROFILER_WINDOW_FRAMES = 60;
// Most lines the profiler overlay shows (the body counts, the ticks/frame line, a header, then one line per section:
// F5 has 15 sections). 18 is also the most that fits above the sliders at 1280 x 720 (see the static_assert in
// TransferScenes.cpp): a 16th F5 section needs the overlay made shorter, e.g. merged rows or a smaller line height.
constexpr size_t MAX_PROFILER_LINES = 18;

const double MIN_ZOOM = 0.01;

const double MAX_ZOOM = 10.0;

const double STARTUP_ZOOM_VALUE = 0.35;

const double STAR_PARALLAX_FACTOR = 0.5;
// How fast the arrow keys pan the camera, in UI points per second (screen points at 1280 x 720; it grows with the
// window like the UI), so it feels the same at any zoom and any window size
const double CAMERA_KEY_PAN_SPEED = 900.0;

constexpr int MAX_VELOCITY_VECTOR_VERTICES = 9; // vertices on a previewbody vector
