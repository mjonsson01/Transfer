// File: Transfer/src/Core/Game.h

#pragma once

// Custom Imports
#include "Core/GameState.h"
#include "Core/UIState.h"
#include "Systems/RenderSystem.h"
#include "Systems/PhysicsSystem.h"
#include "Systems/InputSystem.h"
#include "Systems/UISystem.h"
#include "Utilities/EngineConstants.h"
#include "Utilities/GameSystemConstants.h"

// Standard Library Imports
#include <numeric>

class Game
{
	public: 
		// Constructor and Destructor
		Game();
		~Game();

		// Initializes SDL windows, renderer, and starts the main game loop by calling Run()
		void StartGame();
		// Tears down the 'systems' and cleans up allocated resources.
		void EndGame();
		// Game Loop Entry Point
		void Run();

	private:
		// Core Game Loop Methods
		void ProcessInput(); // Handles User Input from Keyboard and Mouse Events
		void UpdatePhysicsFrame();  // Updates Game State and Physics
		void RenderFrame();  // Renders the Current Frame to the Screen

		// Helpers for Run() method
		void UpdateFPS(Uint32 renderEnd, Uint32 lastRender, float& fpsAccumulator, float& currentFPS);
		void LimitFrameRate(Uint32 renderStart, Uint32 renderEnd);


	private:
		// Systems and State
		GameState state;		     // Contains all game entities and their states
		UIState UIState;             // Contains all UI related states
		InputSystem inputSystem;     // Manages all user input
		// UISystem UISystem;         	 // Manages UI logic and state - now inside RenderSystem
		PhysicsSystem physicsSystem; // Manages physics calculations and Frame Updates
		RenderSystem renderSystem;	 // Manages all rendering operations
};