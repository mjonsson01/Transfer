// File: Transfer/src/Systems/Rendering/RenderSystem.hpp

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3/SDL_pixels.h>
#include <SDL3/SDL_render.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Camera/CameraState.hpp"
#include "Core/GameState.hpp"
#include "Systems/Input/PlayerInput.hpp"
#include "Core/UIState.hpp"
#include "DynamoEngine/Rendering/FontAtlas.hpp"
#include "DynamoEngine/Rendering/UIVertex.hpp"
#include "DynamoEngine/Scenes/Scene.hpp"
#include "DynamoEngine/UI/UIGeometryBuilder.hpp"
#include "DynamoEngine/UI/UIRoot.hpp"
#include "Systems/Rendering/TwinklingStars.hpp"
#include "Constants/EngineConstants.hpp"
#include "Constants/GameSystemConstants.hpp"
#include "Systems/Physics/FragmentLayout.hpp"
#include "Systems/Rendering/CameraData.hpp"
#include "Camera/CameraTransform.hpp"
#include "Systems/Rendering/GPUTypes.hpp"
#include "Utilities/SystemPathUtility.hpp"

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
    ~RenderSystem(); // make sure to teardown destructor and m_window

    // Main Loop Rendering Function, renders engine state and the current scene's UI
    // `player_input`: for the see-through preview of a spawn being aimed
    void renderFullFrame(GameState& game_state, UIState& ui_state, const PlayerInput& player_input,
                         const DynamoEngine::Scene& scene);

    // Main Cleanup method (tears down all the SDL components)
    void cleanUp();
    // Getters for SDL Components
    TTF_Font* getUIFontRegular() const { return UIFontRegular; }

  private:
    // SDL Components
    SDL_Window* m_window = nullptr;
    SDL_GPUDevice* m_gpu = nullptr;

    // Unified Body Rendering Components
    std::vector<UnifiedBodyVertex> m_unified_body_vertices;
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

    // Debug Overlay Rendering Components (F3: hitboxes)
    std::vector<DebugLineVertex> m_debug_line_vertices;
    SDL_GPUBuffer* m_debug_line_vertex_buffer = nullptr;
    SDL_GPUTransferBuffer* m_debug_line_transfer_buffer = nullptr;
    SDL_GPUGraphicsPipeline* m_debug_line_pipeline = nullptr;

    // Planning mode's trajectory lines (TrajectoryLine shaders; the fade-in happens in the shader).
    // Traced lines are append-only: each frame uploads only the segments added since the last frame.
    SDL_GPUGraphicsPipeline* m_trajectory_line_pipeline = nullptr;
    SDL_GPUBuffer* m_trajectory_line_vertex_buffer = nullptr;
    SDL_GPUTransferBuffer* m_trajectory_line_transfer_buffer = nullptr;
    std::vector<TrajectoryLineVertex> m_trajectory_new_vertices; // this frame's new traced segments, to upload
    uint32_t m_trajectory_uploaded_vertex_count = 0;             // traced-line vertices already in the GPU buffer
    bool m_trajectory_restart_upload = false; // the next upload starts the buffer over (a new trace began)
    uint32_t m_trajectory_trace_id_uploaded = 0; // which trace (Trajectories::trace_id) the buffer holds
    // Per traced path (the macro bodies in order, then the ship): how many points already have their segment uploaded
    std::vector<size_t> m_trajectory_points_uploaded;
    // Everything that changes every frame (the placement preview, the time marks, the bright tips): rebuilt and
    // uploaded whole each frame. Small next to the traced lines.
    SDL_GPUBuffer* m_trajectory_frame_vertex_buffer = nullptr;
    SDL_GPUTransferBuffer* m_trajectory_frame_transfer_buffer = nullptr;
    std::vector<TrajectoryLineVertex> m_trajectory_frame_vertices;

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

    // The see-through body (or cluster) that would spawn on release, while a spawn is being aimed
    void appendPreviewBodies(std::vector<UnifiedBodyVertex>& vertex_data, const PlayerInput& player_input,
                             const SpawnSettings& spawn_settings, double zoom);

    void renderBodies(GameState& game_state, UIState& ui_state, SDL_GPURenderPass* pass,
                      SDL_GPUCommandBuffer* cmdbuf); // Renders all the gravitational
                                                     // bodies (both Macro and Particle)

    void uploadUnifiedBodies(GameState& game_state, UIState& ui_state, const PlayerInput& player_input,
                             SDL_GPUCommandBuffer* cmdbuf);
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
    void createDebugLineGPUBufferAndPipeline();

    // Debug overlay: fills m_debug_line_vertices with the hitboxes when ui_state's debug flag is on (empty otherwise)
    void buildDebugLines(const GameState& game_state, UIState& ui_state);
    void addDebugLine(const DynamoEngine::Vector2D& start, const DynamoEngine::Vector2D& prev_start,
                      const DynamoEngine::Vector2D& end, const DynamoEngine::Vector2D& prev_end, SDL_FColor color);
    void addDebugCircle(const DynamoEngine::Vector2D& center, const DynamoEngine::Vector2D& prev_center, double radius,
                        SDL_FColor color);
    void uploadDebugLines(SDL_GPUCommandBuffer* cmdbuf);
    void renderDebugLines(SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf, const CameraState& camera_state);
    // Trajectory lines (planning mode)
    void createTrajectoryLineBuffersAndPipeline();
    // Collects this frame's NEW traced segments, and rebuilds the per-frame lines (preview, time marks, tips)
    void buildTrajectoryLines(const GameState& game_state, const UIState& ui_state);
    // How a path is drawn: its colour, and whether new pieces fade in (traced lines) or not (the live preview)
    struct TrajectoryStyle
    {
        SDL_FColor color;
        bool fades_in;
    };
    // Appends the segments ending at points first_point .. the last point. `start_previous` is where the first point
    // was one tick earlier, so the line starts exactly at the body as drawn (bodies are drawn interpolated between
    // ticks, even while time is stopped).
    void appendTrajectorySegments(const TracedPath& path, size_t first_point, const DynamoEngine::Vector2D& start_previous,
                                  const TrajectoryStyle& style, std::vector<TrajectoryLineVertex>& vertices_out);
    // A short line across the path every TRAJECTORY_TIME_MARK_SECONDS; `zoom` keeps them the same size on screen
    void appendTrajectoryTimeMarks(const TracedPath& path, const TrajectoryStyle& style, float zoom,
                                   std::vector<TrajectoryLineVertex>& vertices_out);
    // The bright front of a line that is still being traced: its newest segments, brightest at the tip
    void appendTrajectoryTip(const TracedPath& path, std::vector<TrajectoryLineVertex>& vertices_out);
    void uploadTrajectoryLines(SDL_GPUCommandBuffer* cmdbuf);
    void renderTrajectoryLines(SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf, const CameraState& camera_state,
                               float clock);
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
    // The spawn velocity arrow from lineStart to lineEnd (world space). Its thickness and head are sized in SCREEN
    // pixels (divided by `zoom`), so the arrow looks the same at any zoom; only its length follows the drag.
    void buildVelocityVectorGeometry(DynamoEngine::Vector2D lineStart, DynamoEngine::Vector2D lineEnd, double zoom);
    void uploadVelocityVectorVertices(SDL_GPUCommandBuffer* cmdbuf);
    void renderVelocityVectors(SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf, const CameraState& camera_state);
    SDL_Color getColorForProperty(const GravitationalBody& body);
};