// File: Transfer/src/Core/Game.cpp

// Custom Imports
#include "Core/Game.hpp"

Game::Game()
    : gameState(), UIState(), inputSystem(), physicsSystem(), renderSystem(gameState), audioSystem(), UISystem()
{
    // fill in imp here
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
    // Initialize any other useful gameState variables here.
    gameState.SetPlaying(true);

    // Default to starting in the level scene since other scenes are not
    // implemented yet.

    UIState.setCurrentScene(SceneIdentifier::START_MENU_SCENE);
    // UIState.setCurrentScene(SceneIdentifier::GAME_SCENE);
    // UIState.setCurrentScene(SceneIdentifier::TEST_VISUAL_SCENE);
    UIState.setPlaySoundEffects(true);
    UIState.setPlayMusic(true);
    UIState.setRequestedMusicMode(MusicMode::TITLE_THEME);
    // Start the main game loop
    Game::Run();

    // End the game and clean up resources after exiting the loop
    Game::EndGame();
}
// Tears down the 'systems' and cleans up allocated resources.
void Game::EndGame()
{
    inputSystem.CleanUp();
    audioSystem.CleanUp();
    physicsSystem.CleanUp();
    UISystem.CleanUp();
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
    while (gameState.IsPlaying())
    {
        // printf("capacity: %.3zu", gameState.getParticles().capacity());
        Uint64 frame_start = SDL_GetPerformanceCounter();

        Game::updateFPS(frame_start, last_frame_start_tick, fps_time_accumulator, current_fps);
        last_frame_start_tick = frame_start;

        // 1. Profile Input
        Uint64 input_start = SDL_GetPerformanceCounter();
        Game::ProcessInput();
        if (!gameState.IsPlaying())
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
        if (UIState.getCurrentSceneID() == SceneIdentifier::GAME_SCENE)
        {
            physics_time_accumulator += (frame_delta * UIState.getTimeScaleFactor());
        }
        else
        {
            physics_time_accumulator = 0.0f;
        }

        // 4. Profile Physics Integration
        Uint64 phys_total_start = SDL_GetPerformanceCounter();
        while (physics_time_accumulator >= PHYSICS_TIME_STEP &&
               UIState.getCurrentSceneID() == SceneIdentifier::GAME_SCENE)
        {
            Game::IntegratePhysicsFrame();
            physics_time_accumulator -= PHYSICS_TIME_STEP;
        }
        Uint64 phys_total_end = SDL_GetPerformanceCounter();

        // 5. Profile Rendering
        gameState.setAlpha(physics_time_accumulator / PHYSICS_TIME_STEP);

        Uint64 render_start = SDL_GetPerformanceCounter();
        Game::RenderFrame();
        Uint64 render_end = SDL_GetPerformanceCounter();

        // --- FRAME LIMITING ---
        // We limit based on how much work we did since frame_start
        Game::limitFrameRate(frame_start, render_end, perf_freq);

        // Calculate metrics in milliseconds
        float input_time = (float)((input_end - input_start) * 1000) / perf_freq;
        float instantiation_time = (float)((inst_end - inst_start) * 1000) / perf_freq;
        float audio_playback_time = (float)((audio_end - audio_start) * 1000) / perf_freq;
        float physics_time = (float)((phys_total_end - phys_total_start) * 1000) / perf_freq;
        float rendering_time = (float)((render_end - render_start) * 1000) / perf_freq;

        // printf("N=%zu | Rend: %.3f | Phys: %.3f\n", gameState.getParticles().size() +
        // gameState.getMacroBodies().size(),
        //        rendering_time, physics_time);
        // printf("Profile Time [ms] | Input: %.3f | Inst: %.3f | Audio: %.3f | Phys: %.3f | Rend: %.3f\n", input_time,
        //        instantiation_time, audio_playback_time, physics_time, rendering_time);
    }
}

// --------- DISPATCH TO SYSTEM METHODS --------- //

void Game::ProcessInput()
{
    // Dispatch to Input System
    inputSystem.ProcessSystemInputFrame(gameState, UIState);
    UISystem.UpdateUIElements(gameState, UIState);
}

void Game::IntegratePhysicsFrame()
{
    // Dispatch to Physics System
    physicsSystem.UpdateSystemFrame(gameState, UIState);
}

void Game::UpdateInstantiations() { physicsSystem.UpdateGravBodyInstantiations(gameState, UIState); }
void Game::RenderFrame()
{
    // Dispatch to Render System -- renders UI as well.
    Scene* current_scene = UISystem.getScene(UIState.getCurrentSceneID());
    const std::unordered_map<UIElementIdentifier, UIElement*>& UI_elements_in_scene = current_scene->getSceneElements();
    renderSystem.RenderFullFrame(gameState, UIState, UI_elements_in_scene);
}

void Game::PlayAudio() { audioSystem.ProcessSystemAudioFrame(gameState, UIState); }
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

        UIState.setFPS(currentFPS);
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
        std::cout << "Camera Zoom" << gameState.getCameraState().zoom << std::endl;
    }
}