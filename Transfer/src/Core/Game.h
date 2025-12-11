// File: Transfer/src/Core/Game.h

#pragma once

// Custom Imports
#include "Core/GameState.h"
#include "Core/UIState.h"
#include "Core/InputState.h"
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
	public:

	private:
		// Core Game Loop Methods
		void ProcessInput(); // Handles User Input from Keyboard and Mouse Events
		void UpdatePhysicsFrame();  // Updates Game State and Physics
		void RenderFrame();  // Renders the Current Frame to the Screen Including UI

		// Helpers for Run() method
		void UpdateFPS(uint32_t renderEnd, uint32_t lastRender, float& fpsAccumulator, float& currentFPS);
		void LimitFrameRate(uint32_t renderStart, uint32_t renderEnd);


	private:
		// Systems and State
		GameState state;		     // Contains all game entities and their states
		UIState UIState;             // Contains all UI related states
		// InputState inputState;	     // Contains all input related states // Moving to UIState
		InputSystem inputSystem;     // Manages all user input
		// UISystem UISystem;         	 // Manages UI logic and state - now inside RenderSystem
		PhysicsSystem physicsSystem; // Manages physics calculations and Frame Updates
		RenderSystem renderSystem;	 // Manages all rendering operations

	// private:
		// Likely move to the game state?
    	// float m_timeScale = 1.0f; // 1.0 is normal speed, 0.5 is half speed.
};