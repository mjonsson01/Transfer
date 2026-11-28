// File: Transfer/src/Core/Game.cpp

// Custom Imports
#include "Core/Game.hpp"

// Standard Library Imports
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

Game::Game()
    : game_state(), ui_state(), m_player_input(), inputSystem(), predictionPhysics(), renderSystem(game_state),
      audioSystem(), scenes(), m_profiler(PROFILER_WINDOW_FRAMES)
{
    addTransferScenes(scenes, ui_state, game_state);

    m_profile.input = m_profiler.addSection("Input");
    m_profile.to_physics = m_profiler.addSection("To physics");
    m_profile.audio = m_profiler.addSection("Audio");
    m_profile.from_physics = m_profiler.addSection("From physics");
    m_profile.physics_thread = m_profiler.addSection("Physics thread");
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
    // The simulation runs on its own thread from here on, starting from the universe set up above
    m_physics_thread.start(game_state.getUniverseState());
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
    m_physics_thread.stop(); // first: it finishes its tick and ends before anything else is torn down
    inputSystem.cleanUp();
    audioSystem.cleanUp();
    predictionPhysics.cleanUp();
    renderSystem.cleanUp();
}
void Game::Run()
{
    // High-resolution frequency (ticks per second)
    const Uint64 perf_freq = SDL_GetPerformanceFrequency();

    // Timing variables
    Uint64 last_frame_start_tick = SDL_GetPerformanceCounter();

    // Accumulators
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

        // 2. Requests to the physics thread: controls, speed, planning, and the queued changes (spawns, scene edits)
        Uint64 to_physics_start = SDL_GetPerformanceCounter();
        Game::sendRequestsToPhysics();
        Uint64 to_physics_end = SDL_GetPerformanceCounter();
        m_profiler.record(m_profile.to_physics,
                          DynamoEngine::FrameProfiler::millisecondsBetween(to_physics_start, to_physics_end));

        // 3. Audio
        Uint64 audio_start = SDL_GetPerformanceCounter();
        Game::PlayAudio();
        Uint64 audio_end = SDL_GetPerformanceCounter();
        m_profiler.record(m_profile.audio, DynamoEngine::FrameProfiler::millisecondsBetween(audio_start, audio_end));

        // 4. The physics thread's newest universe (it runs the ticks on its own: the frame never waits for them)
        Uint64 from_physics_start = SDL_GetPerformanceCounter();
        Game::takeUniverseFromPhysics();
        Uint64 from_physics_end = SDL_GetPerformanceCounter();
        m_profiler.record(m_profile.from_physics,
                          DynamoEngine::FrameProfiler::millisecondsBetween(from_physics_start, from_physics_end));
        const bool runs_simulation = scenes.currentScene().settings().runs_simulation;

        // 5. Planning mode's trajectory lines: once per frame in REAL time (they grow while the simulation is stopped)
        Uint64 trajectories_start = SDL_GetPerformanceCounter();
        // (While planning, W/S still move the throttle: the physics thread does that, and the copy shows it)
        if (runs_simulation)
        {
            predictionPhysics.updateTrajectories(frame_seconds, game_state.getUniverseState(), m_player_input,
                                                 game_state.getTrajectoriesMutable(), ui_state);
        }
        Uint64 trajectories_end = SDL_GetPerformanceCounter();
        m_profiler.record(m_profile.trajectories,
                          DynamoEngine::FrameProfiler::millisecondsBetween(trajectories_start, trajectories_end));

        // 6. Rendering. Bodies are drawn between where they started and ended the copy's last tick (previousPosition
        //    and position), `alpha` of the way: how far into the next tick real time has got since the copy was
        //    published, at the speed the simulation is really running. Time stopped: the end of the tick.
        float alpha = 1.0f;
        const double effective_speed = ui_state.getEffectiveSimulationSpeed();
        if (runs_simulation && !ui_state.getIsPlanning() && effective_speed > 0.0)
        {
            double seconds_since_publish = static_cast<double>(SDL_GetPerformanceCounter() - m_universe_publish_time) /
                                           static_cast<double>(perf_freq);
            double ticks_since_publish = seconds_since_publish * effective_speed / PHYSICS_TIME_STEP;
            alpha = static_cast<float>(std::clamp(ticks_since_publish, 0.0, 1.0));
        }
        game_state.setAlpha(alpha);

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
        // (the F5 report arrives from the physics thread with a universe: takeUniverseFromPhysics publishes it)
    }
}

// One profiler report as overlay lines: the body counts, a header, then one row per section. The columns line up in
// the monospaced UI font (SpaceMono). `title` names the report in its header.
// A section whose code didn't run at all during the window (switched off, e.g. debris-debris gravity) is left out:
// its 0 ms isn't a measurement. Every row shown is a real timing.
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
        if (section.frames_recorded == 0)
        {
            continue;
        }
        std::snprintf(line, sizeof(line), "%-19s %8.3f %8.3f %8.3f", section.name.c_str(), section.average_ms,
                      section.one_percent_low_ms, section.worst_ms);
        lines.push_back(line);
    }
    return lines;
}

// The overlay has MAX_PROFILER_LINES labels. If a report ever has more lines, the last label says how many rows didn't
// fit instead of silently dropping them (the console print always has every line).
static std::vector<std::string> fitToOverlay(const std::vector<std::string>& lines)
{
    if (lines.size() <= MAX_PROFILER_LINES)
    {
        return lines;
    }
    std::vector<std::string> fitted(lines.begin(), lines.begin() + (MAX_PROFILER_LINES - 1));
    fitted.push_back("(+" + std::to_string(lines.size() - fitted.size()) + " rows: see console)");
    return fitted;
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
    // Physics ticks per frame over the same window: the physics thread runs 120 ticks per game second, so about 2 per
    // 60 FPS frame at 1x speed and 6 at 3x. "behind" says in how many frames it couldn't keep up (it dropped game
    // time beyond MAX_PHYSICS_BACKLOG_TICKS, so the simulation ran slower than the slider asks).
    const double average_ticks = static_cast<double>(m_ticks_in_window) / m_profiler.windowFrames();
    const int behind_percent =
        static_cast<int>(100.0 * m_frames_behind_in_window / static_cast<double>(m_profiler.windowFrames()));
    char ticks_line[128];
    if (m_frames_behind_in_window > 0)
    {
        std::snprintf(ticks_line, sizeof(ticks_line), "Ticks/frame: %.2f avg, %d max, behind %d%%", average_ticks,
                      m_most_ticks_in_a_frame, behind_percent);
    }
    else
    {
        std::snprintf(ticks_line, sizeof(ticks_line), "Ticks/frame: %.2f avg, %d max", average_ticks,
                      m_most_ticks_in_a_frame);
    }
    m_ticks_line = ticks_line;
    m_ticks_in_window = 0;
    m_most_ticks_in_a_frame = 0;
    m_frames_behind_in_window = 0;

    std::vector<std::string> lines =
        profilerReportLines("Frame (ms)", m_profiler.report(), game_state.getUniverseState(), m_ticks_line);
    ui_state.setFrameProfilerLines(fitToOverlay(lines));
    if (ui_state.getProfilerOverlay() == ProfilerOverlay::Frame)
    {
        printProfilerLines(lines, m_profiler.windowFrames(), "frames");
    }
}

void Game::publishPhysicsProfilerReport()
{
    std::vector<std::string> lines =
        profilerReportLines("Physics/tick (ms)", m_physics_report, game_state.getUniverseState(), m_ticks_line);
    ui_state.setPhysicsProfilerLines(fitToOverlay(lines));
    if (ui_state.getProfilerOverlay() == ProfilerOverlay::Physics)
    {
        printProfilerLines(lines, PROFILER_WINDOW_FRAMES, "ticks");
    }
}

// --------- DISPATCH TO SYSTEM METHODS --------- //

void Game::ProcessInput(float frame_seconds)
{
    // Dispatch to Input System (which gives the current scene's UI the first look)
    inputSystem.processSystemInputFrame(game_state, ui_state, scenes, m_player_input, frame_seconds);
}

void Game::sendRequestsToPhysics()
{
    // Edits the scenes queued (e.g. the test scene's collision box), in the order they were queued
    std::vector<std::function<void(UniverseState&)>>& edits = game_state.pendingUniverseEdits();
    for (std::function<void(UniverseState&)>& edit : edits)
    {
        m_physics_thread.queueUniverseEdit(std::move(edit));
    }
    edits.clear();

    // The player's spawn / undo / clear this frame, with the spawn panel's settings as they are right now
    const PlayerCommands& commands = m_player_input.commands;
    if (commands.spawn || commands.undo_newest_spawn || commands.clear_all_bodies)
    {
        m_physics_thread.queueSpawnCommands(m_player_input, ui_state.spawnSettings());
    }

    PhysicsThread::FrameRequest request;
    request.ship_controls = m_player_input.ship;
    request.speed = ui_state.getTimeScaleFactor();
    request.is_planning = ui_state.getIsPlanning();
    request.runs_simulation = scenes.currentScene().settings().runs_simulation;
    m_physics_thread.sendFrameRequest(request);
}

void Game::takeUniverseFromPhysics()
{
    PhysicsThread::Handoff handoff;
    if (!m_physics_thread.takeNewestUniverse(game_state.getUniverseStateMutable(), handoff))
    {
        return; // nothing newer yet: keep drawing the copy we have
    }
    m_universe_publish_time = handoff.publish_time;
    m_profiler.record(m_profile.physics_thread, handoff.busy_ms);
    m_ticks_in_window += handoff.ticks;
    m_most_ticks_in_a_frame = std::max(m_most_ticks_in_a_frame, handoff.ticks);
    if (handoff.dropped_time)
    {
        m_frames_behind_in_window += 1;
    }
    ui_state.setEffectiveSimulationSpeed(handoff.effective_speed);
    if (handoff.has_new_profiler_report)
    {
        m_physics_report = std::move(handoff.profiler_report);
        Game::publishPhysicsProfilerReport();
    }
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