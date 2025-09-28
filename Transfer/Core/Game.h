#pragma once
#include <SDL3/SDL.h>

#include "GameState.h"
#include "RendererSystem.h"
#include "PhysicsSystem.h"

class Game
{
public: 
	// Constructor and Destructor
	Game();
	~Game();

	// Game Loop Entry Point
	void Run();

private:
	// Core Game Loop Methods
	void ProcessInput(); // Handles User Input from Keyboard and Mouse Events
	void UpdateFrame();  // Updates Game State, Physics, and Logic
	void RenderFrame();  // Renders the Current Frame to the Screen

private:
	// SDL Components
	SDL_Window* window;
	SDL_Renderer* renderer;

private:
	// Systems and State
	GameState state;		     // Contains all game entities and their states
	InputSystem inputSystem;     // Manages all user input
	PhysicsSystem physicsSystem; // Manages physics calculations and Frame Updates
	Renderer rendererSystem;	 // Manages all rendering operations
};