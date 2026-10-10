// File: Transfer/src/Core/Game.cpp

// Custom Imports
#include "Core/Game.hpp"

Game::Game()
    : game_state(), ui_state(), inputSystem(), physicsSystem(), renderSystem(game_state), audioSystem(), scenes()
{
    addTransferScenes(scenes, ui_state, game_state);
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

    // Default to starting in the level scene since other scenes are not
    // implemented yet.

    // scenes.requestSwitch(TransferScene::StartMenu);
    scenes.requestSwitch(TransferScene::TestVisual);
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
// Tears down the 'systems' and cleans up allocated resources.
void Game::EndGame()
{
    inputSystem.cleanUp();
    audioSystem.CleanUp();
    physicsSystem.CleanUp();
    renderSystem.CleanUp();
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

        // 1. Profile Input
        Uint64 input_start = SDL_GetPerformanceCounter();
        Game::ProcessInput(frame_seconds);
        if (!game_state.IsPlaying())
            break;
        Uint64 input_end = SDL_GetPerformanceCounter();

        // 2. Profile Instantiations
        Uint64 inst_start = SDL_GetPerformanceCounter();
        Game::UpdateInstantiations();
        Uint64 inst_end = SDL_GetPerformanceCounter();

        // 3. Profile Audio
        Uint64 audio_start = SDL_GetPerformanceCounter();
        Game::PlayAudio();
        Uint64 audio_end = SDL_GetPerformanceCounter();

        // Timekeeping for Physics Logic
        Uint64 now_tick = SDL_GetPerformanceCounter();
        float frame_delta = (float)(now_tick - last_physics_update_tick) / perf_freq;
        last_physics_update_tick = now_tick;

        // Physics Scaling Logic
        const bool runs_simulation = scenes.currentScene().settings().runs_simulation;
        if (runs_simulation)
        {
            physics_time_accumulator += (frame_delta * ui_state.getTimeScaleFactor());
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

        // 4. Profile Physics Integration
        Uint64 phys_total_start = SDL_GetPerformanceCounter();
        while (physics_time_accumulator >= PHYSICS_TIME_STEP && runs_simulation)
        {
            Game::IntegratePhysicsFrame();
            physics_time_accumulator -= PHYSICS_TIME_STEP;
        }
        Uint64 phys_total_end = SDL_GetPerformanceCounter();

        // 5. Profile Rendering
        game_state.setAlpha(physics_time_accumulator / PHYSICS_TIME_STEP);

        Uint64 render_start = SDL_GetPerformanceCounter();
        Game::RenderFrame();
        Uint64 render_end = SDL_GetPerformanceCounter();

        // Scene switches requested this frame (buttons, Esc) happen here, once nothing is using the old scene
        scenes.applyPendingSwitch();

        // --- FRAME LIMITING ---
        // We limit based on how much work we did since frame_start
        Game::limitFrameRate(frame_start, render_end, perf_freq);

        // Calculate metrics in milliseconds
        float input_time = (float)((input_end - input_start) * 1000) / perf_freq;
        float instantiation_time = (float)((inst_end - inst_start) * 1000) / perf_freq;
        float audio_playback_time = (float)((audio_end - audio_start) * 1000) / perf_freq;
        float physics_time = (float)((phys_total_end - phys_total_start) * 1000) / perf_freq;
        float rendering_time = (float)((render_end - render_start) * 1000) / perf_freq;

        // printf("N=%zu | Rend: %.3f | Phys: %.3f\n", game_state.getParticles().size() +
        // game_state.getMacroBodies().size(),
        //        rendering_time, physics_time);
        // printf("Profile Time [ms] | Input: %.3f | Inst: %.3f | Audio: %.3f | Phys: %.3f | Rend: %.3f\n", input_time,
        //        instantiation_time, audio_playback_time, physics_time, rendering_time);
    }
}

// --------- DISPATCH TO SYSTEM METHODS --------- //

void Game::ProcessInput(float frame_seconds)
{
    // Dispatch to Input System (which gives the current scene's UI the first look)
    inputSystem.processSystemInputFrame(game_state, ui_state, scenes, frame_seconds);
}

void Game::IntegratePhysicsFrame()
{
    // Dispatch to Physics System
    physicsSystem.UpdateSystemFrame(game_state, ui_state);
}

void Game::UpdateInstantiations() { physicsSystem.UpdateGravBodyInstantiations(game_state, ui_state); }
void Game::RenderFrame()
{
    // Dispatch to Render System -- renders the current scene's UI as well.
    renderSystem.RenderFullFrame(game_state, ui_state, scenes.currentScene());
}

void Game::PlayAudio() { audioSystem.ProcessSystemAudioFrame(game_state, ui_state); }
// --------- UTILITY METHODS FOR FPS --------- //
void Game::updateFPS(Uint64 renderEnd, Uint64 lastRender, float& fpsAccumulator, float& currentFPS)
{
    static const Uint64 perf_freq = SDL_GetPerformanceFrequency();

    // Calculate the duration of this specific frame in seconds
    float frameTimeSeconds = (float)(renderEnd - lastRender) / perf_freq;

    // Accumulate the time in milliseconds for the update interval logic
    fpsAccumulator += (frameTimeSeconds * 1000.0f);

    // Update the average FPS every FPS_UPDATE_DELTA_MS (e.g., 500ms)
    if (fpsAccumulator >= FPS_UPDATE_DELTA_MS)
    {
        // Avoid division by zero; cap minimum frame time to 1 microsecond
        float safeFrameTime = std::max(frameTimeSeconds, 0.000001f);
        float instantFPS = 1.0f / safeFrameTime;

        // Exponential moving average for smoothing
        currentFPS = (0.9f * currentFPS) + (0.1f * instantFPS);

        // Clamp for stability (prevents massive spikes from affecting the UI)
        float target_fps_max = static_cast<float>(TARGET_FPS) * 1.1f;
        currentFPS = std::min(currentFPS, target_fps_max);

        ui_state.setFPS(currentFPS);
        fpsAccumulator = 0.0f;
    }
}
void Game::limitFrameRate(Uint64 renderStart, Uint64 renderEnd, Uint64 perfFreq)
{
    // 1. Calculate how many ticks our target frame duration is
    // (1.0 / TARGET_FPS) * perfFreq
    const Uint64 targetTicksPerFrame = perfFreq / TARGET_FPS;

    Uint64 frameTicks = renderEnd - renderStart;

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