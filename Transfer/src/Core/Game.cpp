// File: Transfer/src/Core/Game.cpp

#include "Core/Game.h"


// Likely change the resolution to be scalable in the future? default 1920x1080 for now. Will be inside the Render system eventually.
Game::Game()
	:  gameState(), UIState(), inputSystem(), physicsSystem(), renderSystem()
{
	// fill in imp here
}

// Handle destruction of any new allocations. None for now, just default.
Game::~Game()
{
	TTF_Quit();
	SDL_Quit();
}

// Initializes SDL windows, renderer, and starts the main game loop
void Game::StartGame()
{
	// Initialize any other useful gameState variables here.
	gameState.SetPlaying(true);

	// Initialize UI elements (like FPS counter and Sliders) and other vars as we go.
	renderSystem.getUISystem()->InitializeUIElements(UIState);

	// Start the main game loop
	Game::Run();
    
    // End the game and clean up resources after exiting the loop
    Game::EndGame();
}
// Tears down the 'systems' and cleans up allocated resources.
void Game::EndGame()
{
    renderSystem.getUISystem()->DeleteUIElements(UIState);
    renderSystem.CleanUp();
    physicsSystem.CleanUp();
    inputSystem.CleanUp();
}
void Game::Run()
{	
	// Initialize the time management variables
	uint32_t last_physics_update_time = SDL_GetTicks();
	uint32_t last_render_time = 0;
    uint32_t now = 0;
    float frame_delta = 0.0f;
    float scaled_frame_delta = 0.0f;
    uint32_t render_start = 0.0f;
    uint32_t render_end = 0.0f;
	
	// Timing accumulators
	float physics_time_accumulator = 0.0f;
	float fps_time_accumulator = 0.0f;

	// Local FPS variable
	float current_fps = 0.0f;

	// Frame interpolation alpha (dynamic)
	float alpha = gameState.getAlpha();

	while (gameState.IsPlaying()){

		// Poll for SDL Events and Process Input
		Game::ProcessInput();

		// Timekeeping
        now = SDL_GetTicks();
        frame_delta = (now - last_physics_update_time) / 1000.0f;
        last_physics_update_time = now;
        
        // Deal with SLOWMO or SPEEDUP
        scaled_frame_delta = frame_delta * gameState.getTimeScaleFactor(); // Use the new time scale factor
        physics_time_accumulator += scaled_frame_delta;

        // Update Physics (remains untouched by time scaling of rendering, maintaining physics accuracy)
        while (physics_time_accumulator >= PHYSICS_TIME_STEP)
        {
            UpdatePhysicsFrame();
            physics_time_accumulator -= PHYSICS_TIME_STEP;
        }
        // Rendering
        alpha = physics_time_accumulator / PHYSICS_TIME_STEP;
        gameState.setAlpha(alpha);

        render_start = SDL_GetTicks();
        RenderFrame();
        render_end = SDL_GetTicks();

        // FPS Calculation
        UpdateFPS(render_end, last_render_time, fps_time_accumulator, current_fps);
        last_render_time = render_end;

        // Frame limiting (soft limiting)
        LimitFrameRate(render_start, render_end);
    }

}


// --------- DISPATCH TO SYSTEM METHODS --------- //

void Game::ProcessInput()
{
	// Dispatch to Input System
	// inputSystem.ProcessSystemFrame(gameState);
//	UISystem.ProcessUIFrame(gameState, UIState);
    // Updates the UI gameState and 
	inputSystem.ProcessSystemInputFrame(gameState, UIState);
}

void Game::UpdatePhysicsFrame()
{
	// Dispatch to Physics System
	physicsSystem.UpdateSystemFrame(gameState, UIState);
}

void Game::RenderFrame()
{
	// Dispatch to Renderer System -- fills in UI as well.
	renderSystem.RenderFullFrame(gameState, UIState);
}

// --------- UTILITY METHODS FOR FPS --------- //

void Game::UpdateFPS(uint32_t renderEnd, uint32_t lastRender, float& fpsAccumulator, float& currentFPS)
{
    float frameTime = (renderEnd - lastRender) / 1000.0f;
    fpsAccumulator += (renderEnd - lastRender);

    if (fpsAccumulator > FPS_UPDATE_DELTA_MS)
    {
        frameTime = std::max(frameTime, 0.001f);
        currentFPS = 0.9f * currentFPS + 0.1f * (1.0f / frameTime);
        float target_fps_max = static_cast<float>(TARGET_FPS)*1.1f;
        currentFPS = std::min(currentFPS, target_fps_max); // Clamp FPS for stability
        UIState.setFPS(currentFPS);
        fpsAccumulator = 0.0f;
    }
}

void Game::LimitFrameRate(uint32_t renderStart, uint32_t renderEnd)
{
    double frame_duration = static_cast<double>(renderEnd - renderStart);
    if (frame_duration < FRAME_DELAY_MS)
        SDL_Delay(static_cast<uint32_t>(FRAME_DELAY_MS - frame_duration));
}