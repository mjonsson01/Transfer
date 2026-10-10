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
#include "DynamoEngine/Rendering/FontAtlas.hpp"
#include "DynamoEngine/Rendering/UIVertex.hpp"
#include "DynamoEngine/Scenes/Scene.hpp"
#include "DynamoEngine/UI/UIGeometryBuilder.hpp"
#include "DynamoEngine/UI/UIRoot.hpp"
#include "Entities/VisualElements/TwinklingStars.hpp"
#include "Utilities/Constants/EngineConstants.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"
#include "Utilities/Rendering/CameraData.hpp"
#include "Utilities/Rendering/CameraTransform.hpp"
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
    RenderSystem(GameState& game_state);
    ~RenderSystem(); // make sure to teardown destructor and window

    // Main Loop Rendering Function, renders engine state and the current scene's UI
    void RenderFullFrame(GameState& game_state, UIState& ui_state, const DynamoEngine::Scene& scene);

    // Main Cleanup method (tears down all the SDL components)
    void CleanUp();
    // Getters for SDL Components
    TTF_Font* getUIFontRegular() const { return UIFontRegular; }

  private:
    // SDL Components
    SDL_Window* window = nullptr;
    SDL_GPUDevice* gpu = nullptr;

    // Unified Body Rendering Components
    std::vector<UnifiedBodyVertex> unifiedBodyVertices;
    SDL_GPUBuffer* unifiedBodyVertexBuffer = nullptr;
    SDL_GPUTransferBuffer* unifiedBodyTransferBuffer = nullptr;
    SDL_GPUGraphicsPipeline* unifiedBodyPipeline = nullptr;
    uint32_t m_unified_body_capacity = 0; // how many bodies the two buffers above can hold right now

    // Twinkling Star Rendering Components
    std::vector<TwinklingStarVertex> twinklingStarVertices;
    SDL_GPUBuffer* twinklingStarVertexBuffer = nullptr;
    SDL_GPUTransferBuffer* twinklingStarTransferBuffer = nullptr;
    SDL_GPUGraphicsPipeline* twinklingStarPipeline = nullptr;

    // UI Element Rendering Components
    std::vector<DynamoEngine::UIVertex> m_ui_vertices;
    SDL_GPUBuffer* uiVertexBuffer = nullptr;
    SDL_GPUTransferBuffer* uiTransferBuffer = nullptr;
    SDL_GPUGraphicsPipeline* uiPipeline = nullptr;

    // Velocity Vector Rendering Components
    std::vector<VelocityVectorVertex> velocityVectorVertices;
    SDL_GPUBuffer* velocityVectorVertexBuffer = nullptr;
    SDL_GPUTransferBuffer* velocityVectorTransferBuffer = nullptr;
    SDL_GPUGraphicsPipeline* velocityVectorPipeline = nullptr;

    // Player Starship Rendering Components
    std::vector<StarshipVertex> starshipVertices;
    SDL_GPUBuffer* starshipVertexBuffer = nullptr;
    SDL_GPUTransferBuffer* starshipTransferBuffer = nullptr;
    SDL_GPUGraphicsPipeline* starshipPipeline = nullptr;
    SDL_GPUTexture* m_starship_texture = nullptr; // the ship sprite, alpha premultiplied, with a full mipmap chain
    SDL_GPUSampler* m_sprite_sampler = nullptr;   // smooth filtering, between pixels AND between mipmap levels

    // Text Rendering Components
    SDL_GPUTexture* fontAtlasTexture = nullptr;
    SDL_GPUSampler* fontAtlasSampler = nullptr;
    DynamoEngine::FontAtlas fontAtlas;

    // Font for UI Elements that require text
    TTF_Font* UIFontRegular = nullptr;

  private:
    // Subordinate Rendering Functions
    void renderGameFrame(GameState& game_state, UIState& ui_state, const DynamoEngine::UIRoot& ui,
                         SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf);
    void renderNonGameFrame(const DynamoEngine::UIRoot& ui, SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf);

    void appendPreviewBodies(std::vector<UnifiedBodyVertex>& vertexData, UIState& ui_state,
                             const CameraState& camera_state);

    void renderBodies(GameState& game_state, UIState& ui_state, SDL_GPURenderPass* pass,
                      SDL_GPUCommandBuffer* cmdbuf); // Renders all the gravitational
                                                     // bodies (both Macro and Particle)

    void uploadUnifiedBodies(GameState& game_state, UIState& ui_state, SDL_GPUCommandBuffer* cmdbuf);
    // (Re)creates both body buffers to hold `capacity` bodies, releasing the old ones.
    // Returns false if the GPU couldn't provide the memory; nothing may be uploaded then.
    bool createUnifiedBodyBuffers(uint32_t capacity);

    SDL_GPUShader* LoadShader(SDL_GPUDevice* device, const char* base_file_name, uint32_t numSamplers = 0,
                              uint32_t numUniformBuffers = 0);

    void createUnifiedBodyGPUBufferAndPipeline();
    void createUIGPUBufferAndPipeline();
    void createVelocityVectorGPUBufferAndPipeline();
    void createTwinklingStarGPUBufferAndPipeline();
    void createStarshipGPUBufferAndPipeline();
    void createFontAtlasSampler();
    // Loads a PNG (path relative to Assets/) into a GPU texture: alpha premultiplied, full mipmap chain.
    // Returns nullptr and logs why if it fails. The caller owns the texture and must SDL_ReleaseGPUTexture it.
    SDL_GPUTexture* loadSpriteTexture(const std::string& asset_path);
    void createSpriteSampler();
    // Bakes fontAtlas from UIFontRegular for `pixel_scale` (screen pixels per UI point) and uploads it to the GPU,
    // replacing the previous atlas texture
    void rebuildFontAtlas(float pixel_scale);
    void uploadUIVertices(const DynamoEngine::UIRoot& ui, SDL_GPUCommandBuffer* cmdbuf);
    void renderUIElements(SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf, const DynamoEngine::UIRoot& ui);

    void createTwinklingStarField(float fieldMaxWidth, float fieldMaxHeight);
    void uploadTwinklingStarField(SDL_GPUCommandBuffer* cmdbuf);
    void uploadStarship(GameState& game_state, UIState& ui_state, SDL_GPUCommandBuffer* cmdbuf);
    void renderTwinklingStarField(SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf,
                                  const CameraState& camera_state);
    void renderStarship(GameState& game_state, SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf,
                        const CameraState& camera_state);

    CameraConstants buildCameraConstants(const CameraState& camera_state, const DynamoEngine::Vector2D& offset);
    // Utility Rendering Helper Functions
    void buildVelocityVectorGeometry(DynamoEngine::Vector2D lineStart, DynamoEngine::Vector2D lineEnd);
    void uploadVelocityVectorVertices(SDL_GPUCommandBuffer* cmdbuf);
    void renderVelocityVectors(SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf, const CameraState& camera_state);
    SDL_Color getColorForProperty(const GravitationalBody& body);
};