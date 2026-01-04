// File: Transfer/src/Systems/RenderSystem.h

#pragma once

// Standard Library Imports
#include <string>
#include <iostream>
#include <unordered_map>
#include <random>
#include <numeric>
#include <algorithm>
#include <cmath>

// SDL3 Imports
#include "SDL3/SDL.h"
#include "SDL3_ttf/SDL_ttf.h"

// Custom Imports
#include "Core/GameState.h"
#include "Core/UIState.h"
#include "Utilities/Colors.h"
#include "Entities/TwinklingStars.h"
#include "Utilities/EngineConstants.h"
#include "Utilities/GameSystemConstants.h"
#include "Systems/UISystem.h"

// Circle (Gravitational body) texture cache.
struct CircleKey {
    int radius;
    SDL_Color color;

    bool operator==(const CircleKey& other) const {
        return radius == other.radius &&
               color.r == other.color.r &&
               color.g == other.color.g &&
               color.b == other.color.b &&
               color.a == other.color.a;
    }
};

// Hash function for unordered_map
namespace std {
    template <>
    struct hash<CircleKey> {
        size_t operator()(const CircleKey& k) const {
            return ((hash<int>()(k.radius)
                     ^ (hash<int>()(k.color.r) << 1))
                     ^ (hash<int>()(k.color.g) << 1))
                     ^ (hash<int>()(k.color.b) << 1)
                     ^ (hash<int>()(k.color.a) << 1);
        }
    };
}

class RenderSystem
{
	public:
		//Constructor and Destructor
		// No arguments for now, but will need to pass through resolution and other info later
		RenderSystem();
		~RenderSystem(); // make sure to teardown destructor and window
		
		// Main Loop Rendering Function, renders engine state and UI state
		void RenderFullFrame(GameState& state, UIState& UIState);

		// Main Cleanup method (tears down all the SDL components)
		void CleanUp();
		// Getters for SDL Components
		SDL_Renderer* getRenderer() const { return renderer; }
		TTF_Font* getUIFont() const { return UIFont; }

		// Getter for UI System
		UISystem* getUISystem() { return &uiSystem; }

	private:
		//SDL Components
		SDL_Window* window = nullptr;
		SDL_Renderer* renderer = nullptr;
		// Font for UI Elements that require text
		TTF_Font* UIFont = nullptr;

		// Container for background twinkling stars
		std::vector<TwinklingStar> twinklingStars;
		// Container for textures of all background twinkling stars
		std::vector<SDL_Texture*> twinklingStarTextures; // pre-created tiny textures (1-3 px)
			
	private:
		// Subordinate Rendering Functions
		void renderBodies(GameState& state); // Renders all the gravitational bodies (both Macro and Particle)

		// Renders Input Artifacts
		void renderInputArtifacts(GameState& state);

		// Utility Rendering Helper Functions
		SDL_Color getColorForMass(double mass);
		
		// Store for circle textures
		std::unordered_map<CircleKey, SDL_Texture*> circleTextureCache;
		// Circle textures helpers
		SDL_Texture* getCircleTexture(int radius, SDL_Color color);
		SDL_Texture* createCircleTexture(int radius, SDL_Color color);

		// Texture cleanup helper
		void clearCachedCircleTextures();

		// Single Call GenerateStar pattern
		void createStarField(int numStars);
		void createStarTextures();
		void updateStars();
		void renderStars();
	private:
		// Managing system for UI overlay
		UISystem uiSystem;
};