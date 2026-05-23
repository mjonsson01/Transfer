// File: Transfer/src/Systems/RenderSystem.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3/SDL_pixels.h>
#include <SDL3/SDL_render.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "Entities/UIElements/UIElement.hpp"
#include "Entities/VisualElements/TwinklingStars.hpp"
#include "Utilities/Constants/EngineConstants.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"
#include "Utilities/Rendering/CameraData.hpp"
#include "Utilities/Rendering/Colors.hpp"
#include "Utilities/Rendering/UnifiedBodyVertex.hpp"
#include "Utilities/System/SystemPathUtility.hpp"

// Standard Library Imports
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>

class RenderSystem
{
  public:
    // Constructor and Destructor
    //  No arguments for now, but will need to pass through resolution and other
    //  info later
    RenderSystem();
    ~RenderSystem(); // make sure to teardown destructor and window

    // Main Loop Rendering Function, renders engine state and UI state
    void RenderFullFrame(GameState& gameState, UIState& UIState,
                         const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope);

    // Main Cleanup method (tears down all the SDL components)
    void CleanUp();
    // Getters for SDL Components
    TTF_Font* getUIFontRegular() const { return UIFontRegular; }
    TTF_Font* getUIFontTitle() const { return UIFontTitle; }

  private:
    // SDL Components
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;

    SDL_GPUDevice* gpu = nullptr;
    SDL_GPUBuffer* unified_body_vertex_buffer;
    SDL_GPUBuffer* camera_uniform_buffer;
    SDL_GPUGraphicsPipeline* unified_body_pipeline;
    SDL_GPUTransferBuffer* body_transfer_buffer = nullptr;
    SDL_GPUTransferBuffer* camera_transfer_buffer = nullptr;
    // Font for UI Elements that require text
    TTF_Font* UIFontRegular = nullptr;
    TTF_Font* UIFontTitle = nullptr;

  private:
    // Subordinate Rendering Functions
    void renderGameFrame(GameState& gameState, UIState& UIState,
                         const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope);
    void renderNonGameFrame(GameState& gameState, UIState& UIState,
                            const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope);

    void appendPreviewBodies(std::vector<UnifiedBodyVertex>& vertexData, UIState& UIState);

    void renderBodies(GameState& gameState, UIState& UIState); // Renders all the gravitational
                                                               // bodies (both Macro and Particle)
    SDL_GPUShader* LoadShader(SDL_GPUDevice* device, const char* fileName);

    // Renders Drag Lines on Preview Bodies
    void renderDragLine(Vector2D lineStart, Vector2D lineEnd);

    // Renders UI Elements
    void renderUIElements(UIState& UIState,
                          const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope);

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