// File: Transfer/src/Systems/RenderSystem.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3/SDL_render.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Core/GameState.h"
#include "Core/UIState.h"
#include "Entities/UIElements/UIElement.h"
#include "Entities/VisualElements/TwinklingStars.h"
#include "Utilities/Constants/EngineConstants.h"
#include "Utilities/Constants/GameSystemConstants.h"
#include "Utilities/Rendering/Colors.h"
#include "Utilities/System/SystemPathUtility.h"

// Standard Library Imports
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <string>

class RenderSystem
{
  public:
    // Constructor and Destructor
    //  No arguments for now, but will need to pass through resolution and other
    //  info later
    RenderSystem();
    ~RenderSystem(); // make sure to teardown destructor and window

    // Main Loop Rendering Function, renders engine state and UI state
    void RenderFullFrame(GameState& gameState, UIState& UIState, const std::vector<UIElement*>& allUIElementsInScope);

    // Main Cleanup method (tears down all the SDL components)
    void CleanUp();
    // Getters for SDL Components
    SDL_Renderer* getRenderer() const { return renderer; }
    TTF_Font* getUIFontRegular() const { return UIFontRegular; }
    TTF_Font* getUIFontTitle() const { return UIFontTitle; }

  private:
    // SDL Components
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;
    // Font for UI Elements that require text
    TTF_Font* UIFontRegular = nullptr;
    TTF_Font* UIFontTitle = nullptr;

  private:
    // Subordinate Rendering Functions
    void renderPreviewBodies(UIState& UIState);

    void renderBodies(GameState& gameState); // Renders all the gravitational
                                             // bodies (both Macro and Particle)

    // Renders Drag Lines on Preview Bodies
    void renderDragLine(Vector2D lineStart, Vector2D lineEnd);

    // Renders UI Elements
    void renderUIElements(UIState& UIState, std::vector<UIElement*> allUIElementsInScope);

    // Utility Rendering Helper Functions
    SDL_Color getColorForProperty(const GravitationalBody& body);

    // Container for background twinkling stars
    std::vector<TwinklingStar> twinklingStars;
    // Container for textures of all background twinkling stars
    std::vector<SDL_Texture*> twinklingStarTextures; // pre-created tiny textures (1-3 px)
    std::vector<SDL_Texture*> circleTextureCache;
    void buildCircleTextureCache();
    // Texture cleanup helper
    void clearCachedCircleTextures();

    // Single Call GenerateStar pattern
    void createStarField(int numStars);
    void createStarTextures();
    void updateStars();
    void renderStars();
};