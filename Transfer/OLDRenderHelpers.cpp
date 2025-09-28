//#include "RenderHelpers.h"
//void renderCircle(SDL_Renderer* renderer, int x, int y, int radius, Color color)
//{
//	SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);
//	for (int w = 0; w < radius * 2; w++)
//	{
//		for (int h = 0; h < radius * 2; h++)
//		{
//			int dx = radius - w; // horizontal offset
//			int dy = radius - h; // vertical offset
//			if ((dx * dx + dy * dy) <= (radius * radius))
//			{
//				SDL_RenderPoint(renderer, x + dx, y + dy);
//			}
//		}
//	}
//}
//
//
//void RenderFrame(SDL_Renderer* renderer, GameState& state)
//{
//	// render the frame
//	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // black background
//	SDL_RenderClear(renderer);
//	
//	// render the existing bodies
//	for (const auto& body : state.bodies) {
//		renderGravitationalBody(renderer, body);
//	}
//
//
//	// render it live
//	if (state.isDragging && state.bodySelectionValidity)
//	{
//	SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255); // white line
//	SDL_RenderLine(renderer, state.dragStartX, state.dragStartY, state.dragEndX, state.dragEndY);
//	}
//	// render old lines
//	for (const auto& initVel : state.initialVelocities) {
//		Color color = ColorLibrary::White; // green for initial velocity lines
//		SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);
//		SDL_RenderLine(renderer, initVel.x_init, initVel.y_init, initVel.x_end, initVel.y_end);
//	}
//	SDL_Delay(TIME_STEP_MS); // 62.5 FPS cap
//	SDL_RenderPresent(renderer);
//}
//
//void renderGravitationalBody(SDL_Renderer* renderer, const GravitationalBody& body)
//{
//	Color color = GetColorForMass(body.mass);
//	renderCircle(renderer, body.x, body.y, body.radius, color);
//}