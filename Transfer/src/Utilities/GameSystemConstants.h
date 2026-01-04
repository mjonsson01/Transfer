// File: Transfer/src/Utilities/GameSystemConstants.h

#pragma once

// FPS target for system
constexpr int TARGET_FPS = 165;

// Target frame delay in ms
const double FRAME_DELAY_MS = 1000.0 / TARGET_FPS;

// Delay to update the text of the rolling average frame counter UI Element
constexpr int FPS_UPDATE_DELTA_MS = 250;

// Time scaling factors to control rendering speed without changing the Physics System behaviors
const double SLOW_TIME_SCALE_FACTOR = 0.0625; // (1/16)x speed
const double REGULAR_TIME_SCALE_FACTOR = 1.0f; // 1x speed 
const double FAST_TIME_SCALE_FACTOR = 2.0f; // 2x speed

// Delay between particle formation when holding to spawn
constexpr int SPAWN_DELAY_MS = 10;

// Vertical and horizontal resolutions
constexpr int SCREEN_HEIGHT = 1080;
constexpr int SCREEN_WIDTH = 1920;

// Background Star count
constexpr int STAR_NUM = 200;