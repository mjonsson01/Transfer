// File: Transfer/src/Systems/RenderSystem.hpp

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3/SDL_pixels.h>
#include <SDL3/SDL_render.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Core/CameraState.hpp"
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "Entities/UIElements/UIElement.hpp"
#include "Entities/VisualElements/TwinklingStars.hpp"
#include "Utilities/Constants/EngineConstants.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"
#include "Utilities/Rendering/CameraData.hpp"
#include "Utilities/Rendering/CameraTransform.hpp"
#include "Utilities/Rendering/Colors.hpp"
#include "Utilities/Rendering/FontAtlasUtility.hpp"
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
#include <vector>

class RenderSystem
{
  public:
    // Constructor and Destructor
    //  No arguments for now, but will need to pass through resolution and other
    //  info later
    RenderSystem(GameState& gameState);
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

    // Unified Body Rendering Components
    std::vector<UnifiedBodyVertex> unifiedBodyVertices;
    SDL_GPUBuffer* unifiedBodyVertexBuffer = nullptr;
    SDL_GPUGraphicsPipeline* unifiedBodyPipeline = nullptr;
    SDL_GPUTransferBuffer* unifiedBodyTransferBuffer = nullptr;

    // Twinkling Star Rendering Components
    std::vector<TwinklingStarVertex> twinklingStarVertices;
    SDL_GPUBuffer* twinklingStarVertexBuffer = nullptr;
    SDL_GPUGraphicsPipeline* twinklingStarPipeline = nullptr;
    SDL_GPUTransferBuffer* twinklingStarTransferBuffer = nullptr;

    // UI Element Rendering Components
    std::vector<UIElementVertex> uiVertices;
    SDL_GPUBuffer* uiVertexBuffer = nullptr;
    SDL_GPUGraphicsPipeline* uiPipeline = nullptr;
    SDL_GPUTransferBuffer* uiTransferBuffer = nullptr;

    // Velocity Vector Rendering Components
    std::vector<VelocityVectorVertex> velocityVectorVertices;
    SDL_GPUBuffer* velocityVectorVertexBuffer = nullptr;
    SDL_GPUGraphicsPipeline* velocityVectorPipeline = nullptr;
    SDL_GPUTransferBuffer* velocityVectorTransferBuffer = nullptr;

    // Player Starship Rendering Components
    std::vector<StarshipVertex> starshipVertices;
    SDL_GPUBuffer* starshipVertexBuffer = nullptr;
    SDL_GPUGraphicsPipeline* starshipPipeline = nullptr;
    SDL_GPUTransferBuffer* starshipTransferBuffer = nullptr;

    // Text Rendering Components
    SDL_GPUTexture* fontAtlasTexture = nullptr;
    SDL_GPUSampler* fontAtlasSampler = nullptr;
    FontAtlasUtility fontAtlas;

    // Font for UI Elements that require text
    TTF_Font* UIFontRegular = nullptr;
    TTF_Font* UIFontTitle = nullptr;

  private:
    // Subordinate Rendering Functions
    void renderGameFrame(GameState& gameState, UIState& UIState,
                         const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope,
                         SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf);
    void renderNonGameFrame(GameState& gameState, UIState& UIState,
                            const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope,
                            SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf);

    void appendPreviewBodies(std::vector<UnifiedBodyVertex>& vertexData, UIState& UIState,
                             const CameraState& cameraState);

    void renderBodies(GameState& gameState, UIState& UIState, SDL_GPURenderPass* pass,
                      SDL_GPUCommandBuffer* cmdbuf); // Renders all the gravitational
                                                     // bodies (both Macro and Particle)

    void uploadBodies(GameState& gameState, UIState& UIState, SDL_GPUCommandBuffer* cmdbuf);
    SDL_GPUShader* LoadShader(SDL_GPUDevice* device, const char* baseFileName, uint32_t numSamplers = 0,
                              uint32_t numUniformBuffers = 0);

    void createGravBodyGPUBuffer();
    void createTwinklingStarGPUBuffer();

    void createUIPipeline();
    void createVelocityVectorPipeline(); // creates pipeline for Velocity vectors in preview body.
    void createFontAtlasTexture();       // bakes fontAtlas from UIFontRegular and uploads it to the GPU
    void uploadUIVertices(const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope,
                          SDL_GPUCommandBuffer* cmdbuf);
    void renderUIElements(SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf, const CameraState& cameraState);

    void createTwinklingStarField(float fieldMaxWidth, float fieldMaxHeight);
    void uploadTwinklingStarField(SDL_GPUCommandBuffer* cmdbuf);
    void renderTwinklingStarField(SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf,
                                  const CameraState& cameraState);

    CameraConstants buildCameraConstants(const CameraState& cameraState, const Vector2D& offset);
    // Utility Rendering Helper Functions
    void buildVelocityVectorGeometry(Vector2D lineStart, Vector2D lineEnd);
    void uploadVelocityVectorVertices(SDL_GPUCommandBuffer* cmdbuf);
    void renderVelocityVector(SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf, const CameraState& cameraState);
    SDL_Color getColorForProperty(const GravitationalBody& body);
};