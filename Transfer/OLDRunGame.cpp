//#include "RunGame.h"
//
//
//void StartGame()
//{
//	SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS | SDL_INIT_AUDIO);
//	SDL_Window* window = SDL_CreateWindow("title", 1920, 1080, 0);
//	SDL_Renderer* renderer = SDL_CreateRenderer(window, nullptr);
//	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // black background
//	SDL_RenderClear(renderer);
//	
//	GameState state; //instantiate game state
//	state.playing = true;
//	state.beginPhysicsSimulation = false;
//	state.isDragging = false; // in case it wasn't instantiated as false.
//	while (state.playing){
//		//std::cout << "Game state isDragging: " << state.isDragging << std::endl;
//		PlayGame(state); // updates game state by 1 frame
//		RenderFrame(renderer, state);
//	}
//	SDL_Delay(5000);
//	SDL_DestroyRenderer(renderer);
//	SDL_DestroyWindow(window);
//	SDL_Quit();
//}
//
//
//void PlayGame(GameState& state)
//{	
//	SDL_Event event;
//	while (SDL_PollEvent(&event)) {
//		if (event.type == SDL_EVENT_QUIT) {
//			state.playing = false;
//		}
//		HandleBodyPropertyKeyEvent(event, state.selectedRadius, state.selectedMass);
//		// check validity after handling key events
//		CheckBodyInstantiationValidity(state);
//
//
//		if (event.type == SDL_EVENT_MOUSE_BUTTON_DOWN && state.bodySelectionValidity) {
//			//slingShotGravitationalBody(renderer, body); // need to pass renderer and body here.
//			state.isDragging = true;
//			state.dragStartX = static_cast<int>(event.button.x);
//			state.dragStartY = static_cast<int>(event.button.y);
//			state.dragEndX = state.dragStartX; // initialize drag end to start
//			state.dragEndY = state.dragStartY; // initialize drag end to start
//		}
//
//		if (event.type == SDL_EVENT_MOUSE_MOTION && state.bodySelectionValidity && state.isDragging) {
//			state.dragEndX = static_cast<int>(event.motion.x);
//			state.dragEndY = static_cast<int>(event.motion.y);
//			// could render a line here from start to end for visual feedback.
//			
//		}
//		if (event.type == SDL_EVENT_MOUSE_BUTTON_UP && state.bodySelectionValidity && state.isDragging) {
//			
//			double velocityScale = 0.5;
//			double v_x = (state.dragEndX - state.dragStartX) * velocityScale;
//			double v_y = (state.dragEndY - state.dragStartY) * velocityScale;
//			int x_init = state.dragStartX;
//			int y_init = state.dragStartY;
//			int x_end = state.dragEndX;
//			int y_end = state.dragEndY;
//			state.initialVelocities.push_back({ x_init, y_init, x_end, y_end, {v_x, v_y} });
//			// instantiate the body at the start position
//			// refactor as helper function later.
//			// int xPos = static_cast<int>(event.button.x);
//			// int yPos = static_cast<int>(event.button.y);
//			int xPos = state.dragStartX;
//			int yPos = state.dragStartY;
//
//			GravitationalBody newBody;
//			newBody.x = xPos;
//			newBody.y = yPos;
//			newBody.radius = state.selectedRadius;
//			newBody.mass = state.selectedMass;
//			newBody.netVelocity = { v_x, v_y};
//			newBody.netForce = { 0.0, 0.0 };
//			state.bodies.push_back(newBody);
//			state.selectedRadius = 0;
//			state.selectedMass = 0;
//
//			// reset the state vars
//			state.isDragging = false;
//		}
//		// somehow call the slingShotGravitationalBody function here with mouse holds.
//
//
//		// toggle physics simulation with 'P' key
//		if (event.type == SDL_EVENT_KEY_DOWN && event.key.scancode == SDL_SCANCODE_P) {
//			state.beginPhysicsSimulation = !state.beginPhysicsSimulation; // play/pause toggle
//			std::cout << "toggling physics simulation" << std::endl;
//			
//			// empty initial velocities
//			state.initialVelocities.clear();
//		}
//		
//	}
//	if (state.beginPhysicsSimulation && state.bodies.size() >= 2) {
//		applyGravityToSystem(state.bodies);
//		updatePhysicsFrame(state.bodies);
//	}
//}
//
//void CheckBodyInstantiationValidity(GameState& state)
//{
//	if ((state.selectedRadius > 0) && (state.selectedMass > 0))
//	{
//		state.bodySelectionValidity = true;
//	}
//	else
//	{
//		state.bodySelectionValidity = false;
//	}
//}