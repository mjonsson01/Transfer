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
		void ProcessInput(); 		// Handles User Input from keyboard and mouse events
		void UpdatePhysicsFrame();  // Updates Game State and Physics
		void RenderFrame();  		// Renders the current frame to the screen including UI

		// Helpers for Run() method
		void UpdateFPS(uint32_t renderEnd, uint32_t lastRender, float& fpsAccumulator, float& currentFPS); 	// Rolling average frame calculation and update
		void LimitFrameRate(uint32_t renderStart, uint32_t renderEnd); 									    // Soft Clamps frame rate based on static Game System Constant value


	private:
		// Systems and State
		GameState gameState;		     // Contains all game entities and their states
		UIState UIState;             // Contains all UI related states
		InputSystem inputSystem;     // Manages all user input
		PhysicsSystem physicsSystem; // Manages physics calculations and Frame Updates
		RenderSystem renderSystem;	 // Manages all rendering operations
};