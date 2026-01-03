// File: Transfer/src/Utilities/GameSystemConstants.h

#pragma once

constexpr int TARGET_FPS = 165;

const double FRAME_DELAY_MS = 1000.0 / TARGET_FPS;

const double SLOW_TIME_SCALE_FACTOR = 0.0625; // 1/16 speed

const double REGULAR_TIME_SCALE_FACTOR = 1.0f;

// add speedup time scale factor?

constexpr int FPS_UPDATE_DELTA_MS = 250;

constexpr int SPAWN_DELAY_MS = 20;

constexpr int SCREEN_HEIGHT = 1080;
constexpr int SCREEN_WIDTH = 1920;

// update based on user prefs
constexpr int STAR_NUM = 200;