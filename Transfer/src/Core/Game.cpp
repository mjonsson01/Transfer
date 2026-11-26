// File: Transfer/src/Core/Game.cpp

// Custom Imports
#include "Core/Game.hpp"

// Standard Library Imports
#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

Game::Game()
    : game_state(), ui_state(), m_player_input(), inputSystem(), physicsSystem(), renderSystem(game_state),
      audioSystem(), scenes(), m_profiler(PROFILER_WINDOW_FRAMES)
{
    addTransferScenes(scenes, ui_state, game_state);

    m_profile.input = m_profiler.addSection("Input");
    m_profile.spawns = m_profiler.addSection("Spawns");
    m_profile.audio = m_profiler.addSection("Audio");
    m_profile.physics = m_profiler.addSection("Physics");
    m_profile.trajectories = m_profiler.addSection("Trajectories");
    m_profile.render = m_profiler.addSection("Render (CPU)");
    m_profile.frame_work = m_profiler.addSection("Frame work");
    m_profile.frame_total = m_profiler.addSection("Frame total");
}

// Handle destruction of any new allocations. None for now, just default.
Game::~Game()
{
    TTF_Quit();
    SDL_Quit();
}

// Initializes SDL windows and starts the main game loop
void Game::StartGame()
{
    // Initialize any other useful game_state variables here.
    game_state.SetPlaying(true);
    game_state.getUniverseStateMutable().has_ship = RUN_GAME_WITH_PLAYER_SHIP;
    setUpWorldEdges();

    // Default to starting in the level scene since other scenes are not
    // implemented yet.

    scenes.requestSwitch(TransferScene::StartMenu);
    // scenes.requestSwitch(TransferScene::TestVisual);
    scenes.applyPendingSwitch(); // start there right away, before the first frame
    ui_state.setPlaySoundEffects(true);
    ui_state.setPlayMusic(true);
    // ui_state.setRequestedMusicMode(MusicMode::TITLE_THEME);
    ui_state.setRequestedMusicMode(MusicMode::MAIN_SHUFFLE);
    // Start the main game loop
    Game::Run();

    // End the game and clean up resources after exiting the loop
    Game::EndGame();
}
// The world's edges: the star field's extent, centered where the star field is (the screen's center at zoom 1).
// Runs after the RenderSystem measured the display (it sizes the star field). The camera can't pan past them; the
// test scene can make them reflective walls (REFLECTIVE_WORLD_WALLS) while it is active.
void Game::setUpWorldEdges()
{
    const CameraState& camera_state = game_state.getCameraState();
    const DynamoEngine::Vector2D center = {SCREEN_WIDTH / 2.0, SCREEN_HEIGHT / 2.0};
    const DynamoEngine::Vector2D half_size = {camera_state.max_display_width / (2.0 * MIN_ZOOM),
                                              camera_state.max_display_height / (2.0 * MIN_ZOOM)};
    UniverseState& universe_state = game_state.getUniverseStateMutable();
    universe_state.world_min_corner = center - half_size;
    universe_state.world_max_corner = center + half_size;
    universe_state.reflective_walls = false;
}

// Tears down the 'systems' and cleans up allocated resources.
void Game::EndGame()
{
    inputSystem.cleanUp();
    audioSystem.cleanUp();
    physicsSystem.cleanUp();
    renderSystem.cleanUp();
}
void Game::Run()
{
    // High-resolution frequency (ticks per second)
    const Uint64 perf_freq = SDL_GetPerformanceFrequency();

    // Timing variables
    Uint64 last_frame_start_tick = SDL_GetPerformanceCounter();
    Uint64 last_physics_update_tick = last_frame_start_tick;

    // Accumulators
    float physics_time_accumulator = 0.0f;
    float fps_time_accumulator = 0.0f;
    float current_fps = (float)TARGET_FPS;
    while (game_state.IsPlaying())
    {
        Uint64 frame_start = SDL_GetPerformanceCounter();
        float frame_seconds = (float)(frame_start - last_frame_start_tick) / perf_freq; // how long the last frame took

        Game::updateFPS(frame_start, last_frame_start_tick, fps_time_accumulator, current_fps);
        last_frame_start_tick = frame_start;

        // 1. Input
        Uint64 input_start = SDL_GetPerformanceCounter();
        Game::ProcessInput(frame_seconds);
        if (!game_state.IsPlaying())
            break;
        Uint64 input_end = SDL_GetPerformanceCounter();
        m_profiler.record(m_profile.input, DynamoEngine::FrameProfiler::millisecondsBetween(input_start, input_end));

        // 2. Spawns (instantiations)
        Uint64 inst_start = SDL_GetPerformanceCounter();
        Game::UpdateInstantiations();
        Uint64 inst_end = SDL_GetPerformanceCounter();
        m_profiler.record(m_profile.spawns, DynamoEngine::FrameProfiler::millisecondsBetween(inst_start, inst_end));

        // 3. Audio
        Uint64 audio_start = SDL_GetPerformanceCounter();
        Game::PlayAudio();
        Uint64 audio_end = SDL_GetPerformanceCounter();
        m_profiler.record(m_profile.audio, DynamoEngine::FrameProfiler::millisecondsBetween(audio_start, audio_end));

        // Timekeeping for Physics Logic
        Uint64 now_tick = SDL_GetPerformanceCounter();
        float frame_delta = (float)(now_tick - last_physics_update_tick) / perf_freq;
        last_physics_update_tick = now_tick;

        // Physics Scaling Logic
        const bool runs_simulation = scenes.currentScene().settings().runs_simulation;
        if (runs_simulation)
        {
            // Planning mode stops time; the slider keeps its value, so leaving planning resumes at that speed
            const float time_scale = ui_state.getIsPlanning() ? 0.0f : ui_state.getTimeScaleFactor();
            physics_time_accumulator += (frame_delta * time_scale);
            // Cap catchup to 8 steps maximum to prevent spiralling frozen frames
            const float max_physics_backlog = MAX_PHYSICS_STEPS_PER_FRAME * PHYSICS_TIME_STEP;
            if (physics_time_accumulator > max_physics_backlog)
            {
                physics_time_accumulator = max_physics_backlog;
            }
        }
        else
        {
            physics_time_accumulator = 0.0f;
        }

        // 4. Physics: every tick this frame
        Uint64 physics_start = SDL_GetPerformanceCounter();
        int ticks_this_frame = 0;
        while (physics_time_accumulator >= PHYSICS_TIME_STEP && runs_simulation)
        {
            Game::IntegratePhysicsFrame();
            physics_time_accumulator -= PHYSICS_TIME_STEP;
            ++ticks_this_frame;
        }
        Uint64 physics_end = SDL_GetPerformanceCounter();
        m_ticks_in_window += ticks_this_frame;
        m_most_ticks_in_a_frame = std::max(m_most_ticks_in_a_frame, ticks_this_frame);
        m_profiler.record(m_profile.physics,
                          DynamoEngine::FrameProfiler::millisecondsBetween(physics_start, physics_end));

        // 5. Planning mode's trajectory lines: once per frame in REAL time (they grow while the simulation is stopped)
        Uint64 trajectories_start = SDL_GetPerformanceCounter();
        if (runs_simulation && ui_state.getIsPlanning())
        {
            // Time is stopped, but W/S still set up the burn (before the prediction, so it sees the new throttle)
            physicsSystem.adjustThrottleWhilePlanning(game_state.getUniverseStateMutable(), m_player_input.ship,
                                                      frame_seconds);
        }
        if (runs_simulation)
        {
            physicsSystem.updateTrajectories(frame_seconds, game_state.getUniverseState(), m_player_input,
                                             game_state.getTrajectoriesMutable(), ui_state);
        }
        Uint64 trajectories_end = SDL_GetPerformanceCounter();
        m_profiler.record(m_profile.trajectories,
                          DynamoEngine::FrameProfiler::millisecondsBetween(trajectories_start, trajectories_end));

        // 6. Rendering
        game_state.setAlpha(physics_time_accumulator / PHYSICS_TIME_STEP);

        Uint64 render_start = SDL_GetPerformanceCounter();
        Game::RenderFrame();
        Uint64 render_end = SDL_GetPerformanceCounter();
        m_profiler.record(m_profile.render, DynamoEngine::FrameProfiler::millisecondsBetween(render_start, render_end));

        // Scene switches requested this frame (buttons, Esc) happen here, once nothing is using the old scene
        scenes.applyPendingSwitch();
        Uint64 work_end = SDL_GetPerformanceCounter();
        m_profiler.record(m_profile.frame_work,
                          DynamoEngine::FrameProfiler::millisecondsBetween(frame_start, work_end));

        // --- FRAME LIMITING ---
        // We limit based on how much work we did since frame_start
        Game::limitFrameRate(frame_start, render_end, perf_freq);

        Uint64 frame_end = SDL_GetPerformanceCounter();
        m_profiler.record(m_profile.frame_total,
                          DynamoEngine::FrameProfiler::millisecondsBetween(frame_start, frame_end));
        m_profiler.endFrame();
        if (m_profiler.takeNewReport())
        {
            Game::publishProfilerReport();
        }
        // (the physics profiler closes its own entries, one per tick)
        if (physicsSystem.profiler().takeNewReport())
        {
            Game::publishPhysicsProfilerReport();
        }
    }
}

// One profiler report as overlay lines: the body counts, a header, then one row per section. The columns line up in
// the monospaced UI font (SpaceMono). `title` names the report in its header.
static std::vector<std::string>
profilerReportLines(const char* title, const std::vector<DynamoEngine::FrameProfiler::SectionReport>& report,
                    const UniverseState& universe_state, const std::string& ticks_line)
{
    std::vector<std::string> lines;
    char line[128];
    // How busy the universe is right now (sampled when the report is made, so it refreshes with the timings)
    std::snprintf(line, sizeof(line), "Particles: %zu   Planets: %zu", universe_state.particles.size(),
                  universe_state.macro_bodies.size());
    lines.push_back(line);
    lines.push_back(ticks_line);
    std::snprintf(line, sizeof(line), "%-19s %8s %8s %8s", title, "avg", "1% low", "worst");
    lines.push_back(line);
    for (const DynamoEngine::FrameProfiler::SectionReport& section : report)
    {
        std::snprintf(line, sizeof(line), "%-19s %8.3f %8.3f %8.3f", section.name.c_str(), section.average_ms,
                      section.one_percent_low_ms, section.worst_ms);
        lines.push_back(line);
    }
    return lines;
}

// Also to the console while that overlay is on, so the numbers can be copied out
// `window_unit`: what the window counts ("frames" for F4, "ticks" for F5)
static void printProfilerLines(const std::vector<std::string>& lines, size_t window_size, const char* window_unit)
{
    std::printf("--- %zu %s ---\n", window_size, window_unit);
    for (const std::string& text : lines)
    {
        std::printf("%s\n", text.c_str());
    }
}

void Game::publishProfilerReport()
{
    // Physics ticks per frame over the same window. The physics runs 120 ticks per game second, so a 60 FPS frame
    // runs 2 at 1x speed; a slow frame owes more game time, so the next frame runs extra "catch-up" ticks (at most
    // MAX_PHYSICS_STEPS_PER_FRAME). If one tick takes longer than the game time it covers, the catch-up never
    // ends and every frame runs the maximum: "AT CAP".
    const double average_ticks = static_cast<double>(m_ticks_in_window) / m_profiler.windowFrames();
    char ticks_line[128];
    std::snprintf(ticks_line, sizeof(ticks_line), "Ticks/frame: %.2f avg, %d max (cap %d)%s", average_ticks,
                  m_most_ticks_in_a_frame, MAX_PHYSICS_STEPS_PER_FRAME,
                  m_most_ticks_in_a_frame >= MAX_PHYSICS_STEPS_PER_FRAME ? "  AT CAP" : "");
    m_ticks_line = ticks_line;
    m_ticks_in_window = 0;
    m_most_ticks_in_a_frame = 0;

    std::vector<std::string> lines =
        profilerReportLines("Frame (ms)", m_profiler.report(), game_state.getUniverseState(), m_ticks_line);
    ui_state.setFrameProfilerLines(lines);
    if (ui_state.getProfilerOverlay() == ProfilerOverlay::Frame)
    {
        printProfilerLines(lines, m_profiler.windowFrames(), "frames");
    }
}

void Game::publishPhysicsProfilerReport()
{
    DynamoEngine::FrameProfiler& physics_profiler = physicsSystem.profiler();
    std::vector<std::string> lines =
        profilerReportLines("Physics/tick (ms)", physics_profiler.report(), game_state.getUniverseState(),
                            m_ticks_line);
    ui_state.setPhysicsProfilerLines(lines);
    if (ui_state.getProfilerOverlay() == ProfilerOverlay::Physics)
    {
        printProfilerLines(lines, physics_profiler.windowFrames(), "ticks");
    }
}

// --------- DISPATCH TO SYSTEM METHODS --------- //

void Game::ProcessInput(float frame_seconds)
{
    // Dispatch to Input System (which gives the current scene's UI the first look)
    inputSystem.processSystemInputFrame(game_state, ui_state, scenes, m_player_input, frame_seconds);
}

void Game::IntegratePhysicsFrame()
{
    // Dispatch to Physics System
    physicsSystem.advanceUniverseState(game_state.getUniverseStateMutable(), m_player_input.ship);
}

void Game::UpdateInstantiations()
{
    physicsSystem.updateGravBodyInstantiations(game_state.getUniverseStateMutable(), m_player_input,
                                               ui_state.spawnSettings());
}
void Game::RenderFrame()
{
    // Dispatch to Render System -- renders the current scene's UI as well.
    renderSystem.renderFullFrame(game_state, ui_state, m_player_input, scenes.currentScene());
}

void Game::PlayAudio() { audioSystem.ProcessSystemAudioFrame(game_state, ui_state); }
// --------- UTILITY METHODS FOR FPS --------- //
void Game::updateFPS(Uint64 render_end, Uint64 last_render, float& fps_accumulator, float& current_fps)
{
    static const Uint64 perf_freq = SDL_GetPerformanceFrequency();

    // Calculate the duration of this specific frame in seconds
    float frameTimeSeconds = (float)(render_end - last_render) / perf_freq;

    // Accumulate the time in milliseconds for the update interval logic
    fps_accumulator += (frameTimeSeconds * 1000.0f);

    // Update the average FPS every FPS_UPDATE_DELTA_MS (e.g., 500ms)
    if (fps_accumulator >= FPS_UPDATE_DELTA_MS)
    {
        // Avoid division by zero; cap minimum frame time to 1 microsecond
        float safe_frame_time = std::max(frameTimeSeconds, 0.000001f);
        float instantFPS = 1.0f / safe_frame_time;

        // Exponential moving average for smoothing
        current_fps = (0.9f * current_fps) + (0.1f * instantFPS);

        // Clamp for stability (prevents massive spikes from affecting the UI)
        float target_fps_max = static_cast<float>(TARGET_FPS) * 1.1f;
        current_fps = std::min(current_fps, target_fps_max);

        ui_state.setFPS(current_fps);
        fps_accumulator = 0.0f;
    }
}
void Game::limitFrameRate(Uint64 renderStart, Uint64 render_end, Uint64 perfFreq)
{
    // 1. Calculate how many ticks our target frame duration is
    // (1.0 / TARGET_FPS) * perfFreq
    const Uint64 targetTicksPerFrame = perfFreq / TARGET_FPS;

    Uint64 frameTicks = render_end - renderStart;

    if (frameTicks < targetTicksPerFrame)
    {
        Uint64 ticksToWait = targetTicksPerFrame - frameTicks;

        // 2. Convert ticks to milliseconds for SDL_Delay
        // We subtract 1ms to avoid oversleeping (SDL_Delay is imprecise)
        uint32_t msToWait = (uint32_t)((ticksToWait * 1000) / perfFreq);

        if (msToWait > 1)
        {
            SDL_Delay(msToWait - 1);
        }

        // 3. Busy-wait for the remaining sub-millisecond precision
        // This ensures we hit the exact target tick
        while (SDL_GetPerformanceCounter() - renderStart < targetTicksPerFrame)
        {
            // Do nothing, just wait out the remaining microseconds
        }
    }
}