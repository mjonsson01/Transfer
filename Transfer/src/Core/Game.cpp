// File: Transfer/src/Core/Game.cpp

#include "Core/Game.h"


// Likely change the resolution to be scalable in the future? default 1920x1080 for now. Will be inside the Render system eventually.
Game::Game()
	:  state(), UIState(), inputSystem(), physicsSystem(), renderSystem()
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
	// Initialize any other useful state variables here.
	state.SetPlaying(true);

	// Initialize UI elements (like FPS counter) and other vars as we go.
	renderSystem.getUISystem()->InitializeUIElements(UIState);

	// Start the main game loop
	Game::Run();
}
// Tears down the 'systems' and cleans up allocated resources.
void Game::EndGame()
{
    // renderSystem.getUISystem()->DeleteUIElements(UIState);
    renderSystem.getUISystem()->DeleteUIElements(UIState);
    renderSystem.CleanUp();
    physicsSystem.CleanUp();
    inputSystem.CleanUp();
}
void Game::Run()
{	
	// Initial time management variables
	Uint32 lastPhysicsUpdateTime = SDL_GetTicks();
	Uint32 lastRenderTime = lastPhysicsUpdateTime;
	
	// Timing accumulators
	float physicsTimeAccumulator = 0.0f;
	float fpsTimeAccumulator = 0.0f;

	// Local fps variable
	float currentFPS = 0.0f;

	// Frame interpolation alpha
	float alpha = 0.0f;

	while (state.IsPlaying()){

		// Poll for SDL Events and Process Input
		Game::ProcessInput();

		// Timekeeping
        Uint32 now = SDL_GetTicks();
        float frameDelta = (now - lastPhysicsUpdateTime) / 1000.0f;
        lastPhysicsUpdateTime = now;
        physicsTimeAccumulator += frameDelta;

        // Update Physics
        while (physicsTimeAccumulator >= PHYSICS_TIME_STEP)
        {
            UpdatePhysicsFrame();
            physicsTimeAccumulator -= PHYSICS_TIME_STEP;
        }

        // Rendering
        alpha = physicsTimeAccumulator / PHYSICS_TIME_STEP;
        state.setAlpha(alpha);

        Uint32 renderStart = SDL_GetTicks();
        RenderFrame();
        Uint32 renderEnd = SDL_GetTicks();

        // FPS Calculation
        UpdateFPS(renderEnd, lastRenderTime, fpsTimeAccumulator, currentFPS);
        lastRenderTime = renderEnd;

        // Frame limiting
        LimitFrameRate(renderStart, renderEnd);
    }

}


// Dispatch function methods.

// prob will dispatch to UI system as well. 

// Thoughts here that maybe inputSystem should still be the one dispatching to the UI? or maybe the UI system is a fully external overlay not tracked
// within the window, resulting in a separate system. Still want to have a shared pointer or maybe singleton instance of the game so that
// the state doesn't need to be passed around so much lol
void Game::ProcessInput()
{
	// Dispatch to Input System
	// inputSystem.ProcessSystemFrame(state);
//	UISystem.ProcessUIFrame(state, UIState);
    // Updates the UI state and 
	inputSystem.ProcessSystemInputFrame(state, UIState);
}

void Game::UpdatePhysicsFrame()
{
	// Dispatch to Physics System
	physicsSystem.UpdateSystemFrame(state);
}

void Game::RenderFrame()
{
	// Dispatch to Renderer System -- fills in UI as well.
	renderSystem.RenderFullFrame(state, UIState);
}

// Helpers to clean up the Run() method


// FPS Helpers
void Game::UpdateFPS(Uint32 renderEnd, Uint32 lastRender, float& fpsAccumulator, float& currentFPS)
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

void Game::LimitFrameRate(Uint32 renderStart, Uint32 renderEnd)
{
    double frameDuration = static_cast<double>(renderEnd - renderStart);
    if (frameDuration < FRAME_DELAY_MS)
        SDL_Delay(static_cast<Uint32>(FRAME_DELAY_MS - frameDuration));
}