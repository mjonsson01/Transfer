// File: Transfer/src/Utilities/Constants/GameSystemConstants.h

#pragma once
// #ifndef GAME_SYSTEM_CONSTANTS_H
// #define GAME_SYSTEM_CONSTANTS_H

// FPS target for system
// constexpr int TARGET_FPS = 165;
constexpr int TARGET_FPS = 120;

// Target frame delay in ms
const double FRAME_DELAY_MS = 1000.0 / TARGET_FPS;

// Delay to update the text of the rolling average frame counter UI Element
constexpr int FPS_UPDATE_DELTA_MS = 250;

// Time scaling factors to control rendering speed without changing the Physics
// System behaviors
const double MIN_TIME_SCALE_FACTOR = 0.0;     // negative is allowed and runs only due to alpha interpolation
const double REGULAR_TIME_SCALE_FACTOR = 1.0; // 1x speed
const double MAX_TIME_SCALE_FACTOR = 2.0;     // 2x speed

// Delay between particle formation when holding to spawn
constexpr int SPAWN_DELAY_MS = 10;

// Vertical and horizontal resolutions (starting resolution)

// can scale up as needed

constexpr int SCREEN_HEIGHT = 1080;
constexpr int SCREEN_WIDTH = 1920;
// constexpr int SCREEN_HEIGHT = 720;
// constexpr int SCREEN_WIDTH = 1280;

// Background Star count
constexpr int STAR_NUM = 100000;

constexpr bool VIEW_DEBUG = false;

const double MIN_ZOOM = 0.01;

const double MAX_ZOOM = 10.0;

const double STARTUP_ZOOM_VALUE = 0.35;

const double STAR_PARALLAX_FACTOR = 0.5;

constexpr int MAX_VELOCITY_VECTOR_VERTICES = 9; // vertices on a previewbody vector

constexpr int NUM_SLIDER_TICKS = 30;

static const double GOLDEN_ANGLE = 2.399963229728653; // pi(3-sqrt(5))