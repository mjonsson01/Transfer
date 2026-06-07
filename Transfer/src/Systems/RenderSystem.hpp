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
#include "Utilities/Rendering/GPUTypes.hpp"
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
    SDL_GPUDevice* gpu = nullptr;

    SDL_GPUBuffer* unifiedBodyVertexBuffer;
    SDL_GPUBuffer* cameraUniformBuffer;
    SDL_GPUBuffer* twinklingStarVertexBuffer;
    SDL_GPUGraphicsPipeline* unifiedBodyPipeline;
    SDL_GPUGraphicsPipeline* twinklingStarPipeline;
    SDL_GPUTransferBuffer* unifiedBodyTransferBuffer = nullptr;
    SDL_GPUTransferBuffer* twinklingStarTransferBuffer = nullptr;
    SDL_GPUTransferBuffer* cameraTransferBuffer = nullptr;
    // Font for UI Elements that require text
    TTF_Font* UIFontRegular = nullptr;
    TTF_Font* UIFontTitle = nullptr;
    std::vector<TwinklingStarVertex> twinklingStars;

  private:
    // Subordinate Rendering Functions
    void renderGameFrame(GameState& gameState, UIState& UIState,
                         const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope,
                         SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf);
    void renderNonGameFrame(GameState& gameState, UIState& UIState,
                            const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope,
                            SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf);

    void appendPreviewBodies(std::vector<UnifiedBodyVertex>& vertexData, UIState& UIState);

    void renderBodies(GameState& gameState, UIState& UIState, SDL_GPURenderPass* pass,
                      SDL_GPUCommandBuffer* cmdbuf); // Renders all the gravitational
                                                     // bodies (both Macro and Particle)

    void uploadBodies(GameState& gameState, UIState& UIState, SDL_GPUCommandBuffer* cmdbuf);
    SDL_GPUShader* LoadShader(SDL_GPUDevice* device, const char* fileName);

    void createGravBodyGPUBuffer();
    void createTwinklingStarGPUBuffer();

    void createTwinklingStarField();
    void uploadTwinklingStarField(SDL_GPUCommandBuffer* cmdbuf);
    void renderTwinklingStarField(SDL_GPURenderPass* pass);
    // Utility Rendering Helper Functions
    SDL_Color getColorForProperty(const GravitationalBody& body);
};