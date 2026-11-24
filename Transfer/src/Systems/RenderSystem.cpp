// File: Transfer/src/Systems/RenderSystem.cpp

#include "Systems/RenderSystem.hpp"

#include "DynamoEngine/Constants/GlobalConstants.hpp"

namespace
{
constexpr float UI_FONT_SIZE = 18.0f; // UI points
} // namespace

// Constructor: Initializes SDL Window and GPU
RenderSystem::RenderSystem(GameState& game_state)
{
    SDL_InitSubSystem(SDL_INIT_VIDEO);
    TTF_Init();

    // 1. Window & Device Setup
    int window_flags = SDL_WINDOW_RESIZABLE | SDL_WINDOW_HIGH_PIXEL_DENSITY;
    m_window = SDL_CreateWindow("Transfer", SCREEN_WIDTH, SCREEN_HEIGHT, window_flags);

    SDL_DisplayID display_id = SDL_GetDisplayForWindow(m_window);
    const SDL_DisplayMode* desktop_mode = SDL_GetDesktopDisplayMode(display_id);
    if (desktop_mode)
    {
        game_state.getCameraStateMutable().max_display_width = (float)desktop_mode->w;
        game_state.getCameraStateMutable().max_display_height = (float)desktop_mode->h;
    }

#ifdef __APPLE__
    SDL_GPUShaderFormat formats = SDL_GPU_SHADERFORMAT_MSL;
#else
    SDL_GPUShaderFormat formats = SDL_GPU_SHADERFORMAT_SPIRV;
#endif
    m_gpu = SDL_CreateGPUDevice(formats, true, nullptr);
    if (m_gpu)
        SDL_ClaimWindowForGPUDevice(m_gpu, m_window);
    // 2. Resource/Font Setup
    UIFontRegular = TTF_OpenFont(Utilities::GetResourcePath("Fonts/SpaceMono-Regular.ttf").c_str(), UI_FONT_SIZE);

    createUnifiedBodyGPUBufferAndPipeline();
    createTwinklingStarGPUBufferAndPipeline();
    createVelocityVectorGPUBufferAndPipeline();
    createUIGPUBufferAndPipeline();
    createFontAtlasSampler(); // the atlas itself is baked on the first frame, once the UI scale is known
    createSpriteSampler();
    m_starship_texture = loadSpriteTexture("Visual/Ships/FutureShip.png");
    createStarshipGPUBufferAndPipeline();
    createDebugLineGPUBufferAndPipeline();
    createTrajectoryLineBuffersAndPipeline();
    createTwinklingStarField(game_state.getCameraState().max_display_width,
                             game_state.getCameraState().max_display_height);
    if (m_gpu)
    {
        SDL_GPUCommandBuffer* initCmdBuf = SDL_AcquireGPUCommandBuffer(m_gpu);
        // The star field never changes (twinkling and parallax happen in the shader), so upload it once, here
        uploadTwinklingStarField(initCmdBuf);
        SDL_SubmitGPUCommandBuffer(initCmdBuf);
    }
}
// Destructor: Cleans up SDL Window
RenderSystem::~RenderSystem()
{
    // TTF_Quit() and SDL_QUIT() handled at the Game level
}

// --------- CLEANUP METHOD --------- //
void RenderSystem::cleanUp()
{
    // Release GPU-specific resources

    // Release Unified Body Resources
    if (unifiedBodyVertexBuffer != nullptr)
        SDL_ReleaseGPUBuffer(m_gpu, unifiedBodyVertexBuffer);
    if (unifiedBodyPipeline != nullptr)
        SDL_ReleaseGPUGraphicsPipeline(m_gpu, unifiedBodyPipeline);
    if (unifiedBodyTransferBuffer != nullptr)
        SDL_ReleaseGPUTransferBuffer(m_gpu, unifiedBodyTransferBuffer);

    // Release Twinkling Star Resources
    if (twinklingStarVertexBuffer != nullptr)
        SDL_ReleaseGPUBuffer(m_gpu, twinklingStarVertexBuffer);
    if (twinklingStarPipeline != nullptr)
        SDL_ReleaseGPUGraphicsPipeline(m_gpu, twinklingStarPipeline);
    if (twinklingStarTransferBuffer != nullptr)
        SDL_ReleaseGPUTransferBuffer(m_gpu, twinklingStarTransferBuffer);

    // Release UI Resources
    if (uiVertexBuffer != nullptr)
        SDL_ReleaseGPUBuffer(m_gpu, uiVertexBuffer);
    if (uiPipeline != nullptr)
        SDL_ReleaseGPUGraphicsPipeline(m_gpu, uiPipeline);
    if (uiTransferBuffer != nullptr)
        SDL_ReleaseGPUTransferBuffer(m_gpu, uiTransferBuffer);

    // Release Velocity Vector Resources
    if (velocityVectorVertexBuffer != nullptr)
        SDL_ReleaseGPUBuffer(m_gpu, velocityVectorVertexBuffer);
    if (velocityVectorPipeline != nullptr)
        SDL_ReleaseGPUGraphicsPipeline(m_gpu, velocityVectorPipeline);
    if (velocityVectorTransferBuffer != nullptr)
        SDL_ReleaseGPUTransferBuffer(m_gpu, velocityVectorTransferBuffer);

    // Release Starship Pipeline
    if (starshipVertexBuffer != nullptr)
        SDL_ReleaseGPUBuffer(m_gpu, starshipVertexBuffer);
    if (starshipPipeline != nullptr)
        SDL_ReleaseGPUGraphicsPipeline(m_gpu, starshipPipeline);
    if (starshipTransferBuffer != nullptr)
        SDL_ReleaseGPUTransferBuffer(m_gpu, starshipTransferBuffer);
    if (m_starship_texture != nullptr)
        SDL_ReleaseGPUTexture(m_gpu, m_starship_texture);
    if (m_sprite_sampler != nullptr)
        SDL_ReleaseGPUSampler(m_gpu, m_sprite_sampler);

    // Release Debug Overlay Pipeline
    if (m_debug_line_vertex_buffer != nullptr)
        SDL_ReleaseGPUBuffer(m_gpu, m_debug_line_vertex_buffer);
    if (m_debug_line_transfer_buffer != nullptr)
        SDL_ReleaseGPUTransferBuffer(m_gpu, m_debug_line_transfer_buffer);
    if (m_debug_line_pipeline != nullptr)
        SDL_ReleaseGPUGraphicsPipeline(m_gpu, m_debug_line_pipeline);

    // Release Trajectory Line Pipeline and Buffers
    if (m_trajectory_line_pipeline != nullptr)
        SDL_ReleaseGPUGraphicsPipeline(m_gpu, m_trajectory_line_pipeline);
    if (m_trajectory_line_vertex_buffer != nullptr)
        SDL_ReleaseGPUBuffer(m_gpu, m_trajectory_line_vertex_buffer);
    if (m_trajectory_line_transfer_buffer != nullptr)
        SDL_ReleaseGPUTransferBuffer(m_gpu, m_trajectory_line_transfer_buffer);
    if (m_trajectory_frame_vertex_buffer != nullptr)
        SDL_ReleaseGPUBuffer(m_gpu, m_trajectory_frame_vertex_buffer);
    if (m_trajectory_frame_transfer_buffer != nullptr)
        SDL_ReleaseGPUTransferBuffer(m_gpu, m_trajectory_frame_transfer_buffer);

    // Release Font Resources
    if (fontAtlasTexture != nullptr)
        SDL_ReleaseGPUTexture(m_gpu, fontAtlasTexture);
    if (fontAtlasSampler != nullptr)
        SDL_ReleaseGPUSampler(m_gpu, fontAtlasSampler);

    // release the window from the gpu
    if (m_gpu != nullptr && m_window != nullptr)
        SDL_ReleaseWindowFromGPUDevice(m_gpu, m_window);

    // destroy the gpu
    if (m_gpu != nullptr)
        SDL_DestroyGPUDevice(m_gpu);

    // destroy the window
    if (m_window)
    {
        SDL_DestroyWindow(m_window);
        m_window = nullptr;
    }
}

SDL_GPUShader* RenderSystem::LoadShader(SDL_GPUDevice* device, const char* base_file_name, uint32_t numSamplers,
                                        uint32_t numUniformBuffers)
{
    size_t size;

#ifdef __APPLE__
    std::string fileName = std::string(base_file_name) + ".msl";
    const char* entrypoint = "main0"; // spirv-cross renames the MSL entry point away from "main"
    SDL_GPUShaderFormat format = SDL_GPU_SHADERFORMAT_MSL;
#else
    std::string fileName = std::string(base_file_name) + ".spv";
    const char* entrypoint = "main";
    SDL_GPUShaderFormat format = SDL_GPU_SHADERFORMAT_SPIRV;
#endif

    // Uses SystemPathUtility to find the shader directory
    std::string fullPath = Utilities::GetResourcePath(fileName);
    void* code = SDL_LoadFile(fullPath.c_str(), &size);

    if (!code)
    {
        std::cerr << "Failed to load shader file: " << fileName << std::endl;
        return nullptr;
    }

    SDL_GPUShaderCreateInfo info = {.code_size = size,
                                    .code = (const uint8_t*)code,
                                    .entrypoint = entrypoint,
                                    .format = format,
                                    .stage = (fileName.find("vert") != std::string::npos)
                                                 ? SDL_GPU_SHADERSTAGE_VERTEX
                                                 : SDL_GPU_SHADERSTAGE_FRAGMENT,
                                    .num_samplers = numSamplers,
                                    .num_uniform_buffers = numUniformBuffers};

    SDL_GPUShader* shader = SDL_CreateGPUShader(device, &info);
    SDL_free(code);
    return shader;
}
// --------- RENDER FULL FRAME METHOD --------- //

void RenderSystem::renderFullFrame(GameState& game_state, UIState& ui_state, const DynamoEngine::Scene& scene)
{

    // Re-bake the text whenever the UI scale or the display's pixel density changes (m_window resized, or moved to a
    // screen with a different density), so text is always drawn at exactly one atlas pixel per screen pixel
    const float pixel_scale = scene.ui().uiScale() * SDL_GetWindowPixelDensity(m_window);
    if (fontAtlasTexture == nullptr || pixel_scale != fontAtlas.pixelScale())
    {
        rebuildFontAtlas(pixel_scale);
    }

    SDL_GPUCommandBuffer* cmdbuf = SDL_AcquireGPUCommandBuffer(m_gpu);

    const bool draws_world = scene.settings().draws_world;

    if (draws_world)
    {
        buildDebugLines(game_state, ui_state);
        uploadDebugLines(cmdbuf);
        buildTrajectoryLines(game_state, ui_state);
        uploadTrajectoryLines(cmdbuf);
        uploadUnifiedBodies(game_state, ui_state, cmdbuf);
        uploadStarship(game_state, ui_state, cmdbuf);
    }
    uploadUIVertices(scene.ui(), cmdbuf);
    // Acquire the display target
    SDL_GPUTexture* swapchainTexture = nullptr;
    Uint32 w = 0, h = 0;
    if (!SDL_AcquireGPUSwapchainTexture(cmdbuf, m_window, &swapchainTexture, &w, &h))
    {
        SDL_SubmitGPUCommandBuffer(cmdbuf);
        return;
    }

    if (swapchainTexture)
    {
        SDL_GPUColorTargetInfo color_info = {};
        color_info.texture = swapchainTexture;
        color_info.clear_color = {0.0f, 0.0f, 0.0f, 1.0f};
        color_info.load_op = SDL_GPU_LOADOP_CLEAR;
        color_info.store_op = SDL_GPU_STOREOP_STORE;

        SDL_GPURenderPass* pass = SDL_BeginGPURenderPass(cmdbuf, &color_info, 1, nullptr);

        if (draws_world)
        {
            renderGameFrame(game_state, ui_state, scene.ui(), pass, cmdbuf);
        }
        else
        {
            renderNonGameFrame(scene.ui(), pass, cmdbuf);
        }

        SDL_EndGPURenderPass(pass);
    }

    SDL_SubmitGPUCommandBuffer(cmdbuf);
}

void RenderSystem::renderGameFrame(GameState& game_state, UIState& ui_state, const DynamoEngine::UIRoot& ui,
                                   SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf)
{
    game_state.getCameraStateMutable().render_alpha = game_state.getAlpha();
    renderTwinklingStarField(pass, cmdbuf, game_state.getCameraState());
    renderTrajectoryLines(pass, cmdbuf, game_state.getCameraState(), game_state.getTrajectories().clock); // under bodies
    renderBodies(game_state, ui_state, pass, cmdbuf);
    renderStarship(game_state, pass, cmdbuf, game_state.getCameraState());
    renderVelocityVectors(pass, cmdbuf, game_state.getCameraState());
    renderDebugLines(pass, cmdbuf, game_state.getCameraState()); // over the whole world, under the UI
    renderUIElements(pass, cmdbuf, ui);
}

void RenderSystem::renderNonGameFrame(const DynamoEngine::UIRoot& ui, SDL_GPURenderPass* pass,
                                      SDL_GPUCommandBuffer* cmdbuf)
{
    renderUIElements(pass, cmdbuf, ui);
}
void RenderSystem::uploadUnifiedBodies(GameState& game_state, UIState& ui_state, SDL_GPUCommandBuffer* cmdbuf)
{
    m_unified_body_vertices.clear();

    auto& particles = game_state.getUniverseState().particles;
    auto& bodies = game_state.getUniverseState().macro_bodies;

    m_unified_body_vertices.reserve(particles.size() + bodies.size());

    for (auto& p : particles)
    {
        if (p.visible)
        {
            m_unified_body_vertices.push_back(p.toUnifiedVertex());
        }
    }

    for (auto& b : bodies)
    {
        if (b.visible)
        {
            m_unified_body_vertices.push_back(b.toUnifiedVertex());
        }
    }

    appendPreviewBodies(m_unified_body_vertices, ui_state, game_state.getCameraState());
    uploadVelocityVectorVertices(cmdbuf);

    // Emergency only: the particle budget should make this impossible. If a rule change ever breaks it, grow the
    // buffers (doubling, so it happens rarely) and say so, instead of hiding bodies or writing past the buffer.
    if (m_unified_body_vertices.size() > m_unified_body_capacity)
    {
        uint32_t new_capacity = std::max(m_unified_body_capacity, INITIAL_UNIFIED_BODY_CAPACITY);
        while (new_capacity < m_unified_body_vertices.size())
        {
            new_capacity *= 2;
        }
        printf("WARNING: %zu bodies don't fit the body buffers (%u), growing them to %u. Is the particle budget "
               "still being enforced?\n",
               m_unified_body_vertices.size(), m_unified_body_capacity, new_capacity);
        if (!createUnifiedBodyBuffers(new_capacity))
        {
            m_unified_body_vertices
                .clear(); // no GPU memory: draw no bodies this frame rather than write past the buffer
            return;
        }
    }
    // Copy pass
    if (!m_unified_body_vertices.empty())
    {
        void* map = SDL_MapGPUTransferBuffer(m_gpu, unifiedBodyTransferBuffer, true);
        SDL_memcpy(map, m_unified_body_vertices.data(), m_unified_body_vertices.size() * sizeof(UnifiedBodyVertex));
        SDL_UnmapGPUTransferBuffer(m_gpu, unifiedBodyTransferBuffer);

        SDL_GPUCopyPass* copyPass = SDL_BeginGPUCopyPass(cmdbuf);

        SDL_GPUTransferBufferLocation src = {.transfer_buffer = unifiedBodyTransferBuffer, .offset = 0};
        SDL_GPUBufferRegion dst = {.buffer = unifiedBodyVertexBuffer,
                                   .offset = 0,
                                   .size = (uint32_t)(m_unified_body_vertices.size() * sizeof(UnifiedBodyVertex))};

        SDL_UploadToGPUBuffer(copyPass, &src, &dst, true);
        SDL_EndGPUCopyPass(copyPass);
    }
}

void RenderSystem::uploadStarship(GameState& game_state, UIState& ui_state, SDL_GPUCommandBuffer* cmdbuf)
{
    starshipVertices.clear();
    starshipVertices.reserve(MAX_STARSHIP_VERTICES);

    game_state.getUniverseStateMutable().ship.buildGeometry(starshipVertices);

    // Never copy more than the buffer holds; renderStarship draws starshipVertices.size(), so it stays in sync
    if (starshipVertices.size() > MAX_STARSHIP_VERTICES)
    {
        printf("Too many starship vertices %zu > %u\n", starshipVertices.size(), MAX_STARSHIP_VERTICES);
        starshipVertices.resize(MAX_STARSHIP_VERTICES);
    }
    if (starshipVertices.empty())
    {
        printf("No starship vertices received");
    }
    if (!starshipVertices.empty())
    {
        void* map = SDL_MapGPUTransferBuffer(m_gpu, starshipTransferBuffer, true);
        SDL_memcpy(map, starshipVertices.data(), starshipVertices.size() * sizeof(StarshipVertex));
        SDL_UnmapGPUTransferBuffer(m_gpu, starshipTransferBuffer);
        SDL_GPUCopyPass* copyPass = SDL_BeginGPUCopyPass(cmdbuf);
        SDL_GPUTransferBufferLocation src = {.transfer_buffer = starshipTransferBuffer, .offset = 0};
        SDL_GPUBufferRegion dst = {.buffer = starshipVertexBuffer,
                                   .offset = 0,
                                   .size = (uint32_t)(starshipVertices.size() * sizeof(StarshipVertex))};

        SDL_UploadToGPUBuffer(copyPass, &src, &dst, true);
        SDL_EndGPUCopyPass(copyPass);
    }
}

void RenderSystem::renderBodies(GameState& game_state, UIState& ui_state, SDL_GPURenderPass* pass,
                                SDL_GPUCommandBuffer* cmdbuf)
{
    // Draw exactly what uploadUnifiedBodies put in the buffer (one instance per body)
    uint32_t instance_count = (uint32_t)m_unified_body_vertices.size();
    if (instance_count > 0)
    {
        SDL_BindGPUGraphicsPipeline(pass, unifiedBodyPipeline);
        CameraConstants camera_constants =
            buildCameraConstants(game_state.getCameraState(), game_state.getCameraState().offset);
        SDL_PushGPUVertexUniformData(cmdbuf, 0, &camera_constants, sizeof(camera_constants));

        SDL_GPUBufferBinding vbo = {.buffer = unifiedBodyVertexBuffer, .offset = 0};
        SDL_BindGPUVertexBuffers(pass, 0, &vbo, 1);

        SDL_DrawGPUPrimitives(pass,
                              6,              // vertices per quad
                              instance_count, // instances
                              0, 0);
    }
}

void RenderSystem::renderStarship(GameState& game_state, SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf,
                                  const CameraState& camera_state)
{
    // No sprite (loadSpriteTexture already logged why at startup) or nothing to draw: skip the ship
    if (m_starship_texture == nullptr || starshipVertices.empty())
    {
        return;
    }
    SDL_BindGPUGraphicsPipeline(pass, starshipPipeline);
    CameraConstants camera_constants =
        buildCameraConstants(game_state.getCameraState(), game_state.getCameraState().offset);
    SDL_PushGPUVertexUniformData(cmdbuf, 0, &camera_constants, sizeof(camera_constants));

    SDL_GPUBufferBinding vbo = {.buffer = starshipVertexBuffer, .offset = 0};
    SDL_BindGPUVertexBuffers(pass, 0, &vbo, 1);

    SDL_GPUTextureSamplerBinding sprite_binding = {.texture = m_starship_texture, .sampler = m_sprite_sampler};
    SDL_BindGPUFragmentSamplers(pass, 0, &sprite_binding, 1);

    SDL_DrawGPUPrimitives(pass,
                          (uint32_t)starshipVertices.size(), // vertices per quad
                          1,                                 // instances
                          0, 0);
}

void RenderSystem::createUIGPUBufferAndPipeline()
{
    SDL_GPUBufferCreateInfo vb_info = {.usage = SDL_GPU_BUFFERUSAGE_VERTEX,
                                       .size = MAX_UI_VERTICES * sizeof(DynamoEngine::UIVertex)};

    uiVertexBuffer = SDL_CreateGPUBuffer(m_gpu, &vb_info);
    SDL_GPUTransferBufferCreateInfo tb_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD,
                                               .size = MAX_UI_VERTICES * sizeof(DynamoEngine::UIVertex)};
    uiTransferBuffer = SDL_CreateGPUTransferBuffer(m_gpu, &tb_info);

    SDL_GPUShader* vert_shader = LoadShader(m_gpu, "Shaders/UIElement.vert", 0, 1);
    SDL_GPUShader* frag_shader = LoadShader(m_gpu, "Shaders/UIElement.frag", 1, 0);

    SDL_GPUVertexAttribute attrs[6];
    attrs[0] = {.location = 0,
                .buffer_slot = 0,
                .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2,
                .offset = offsetof(DynamoEngine::UIVertex, x)};
    attrs[1] = {.location = 1,
                .buffer_slot = 0,
                .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2,
                .offset = offsetof(DynamoEngine::UIVertex, u)};
    attrs[2] = {.location = 2,
                .buffer_slot = 0,
                .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT4,
                .offset = offsetof(DynamoEngine::UIVertex, r)};
    attrs[3] = {.location = 3,
                .buffer_slot = 0,
                .format = SDL_GPU_VERTEXELEMENTFORMAT_UINT,
                .offset = offsetof(DynamoEngine::UIVertex, z_index)};
    attrs[4] = {.location = 4,
                .buffer_slot = 0,
                .format = SDL_GPU_VERTEXELEMENTFORMAT_UINT,
                .offset = offsetof(DynamoEngine::UIVertex, mode)};
    SDL_GPUGraphicsPipelineCreateInfo pipeline_info = {};
    pipeline_info.target_info.num_color_targets = 1;

    SDL_GPUColorTargetDescription color_target = {};
    color_target.format = SDL_GetGPUSwapchainTextureFormat(m_gpu, m_window);
    color_target.blend_state.enable_blend = true;
    color_target.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
    color_target.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
    color_target.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
    color_target.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
    color_target.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
    color_target.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;

    pipeline_info.target_info.color_target_descriptions = &color_target;
    pipeline_info.vertex_shader = vert_shader;
    pipeline_info.fragment_shader = frag_shader;
    pipeline_info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
    pipeline_info.vertex_input_state.vertex_attributes = attrs;
    pipeline_info.vertex_input_state.num_vertex_attributes = 5;

    SDL_GPUVertexBufferDescription vbo_desc = {
        .slot = 0, .pitch = sizeof(DynamoEngine::UIVertex), .input_rate = SDL_GPU_VERTEXINPUTRATE_VERTEX};
    pipeline_info.vertex_input_state.vertex_buffer_descriptions = &vbo_desc;
    pipeline_info.vertex_input_state.num_vertex_buffers = 1;

    uiPipeline = SDL_CreateGPUGraphicsPipeline(m_gpu, &pipeline_info);
    SDL_ReleaseGPUShader(m_gpu, vert_shader);
    SDL_ReleaseGPUShader(m_gpu, frag_shader);
}

void RenderSystem::createFontAtlasSampler()
{
    SDL_GPUSamplerCreateInfo sampler_info = {.min_filter = SDL_GPU_FILTER_LINEAR,
                                             .mag_filter = SDL_GPU_FILTER_LINEAR,
                                             .address_mode_u = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE,
                                             .address_mode_v = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE};
    fontAtlasSampler = SDL_CreateGPUSampler(m_gpu, &sampler_info);
}

void RenderSystem::createSpriteSampler()
{
    // max_lod must be raised from its default of 0, or the GPU is never allowed past mip level 0 (no mipmaps at all)
    SDL_GPUSamplerCreateInfo sampler_info = {.min_filter = SDL_GPU_FILTER_LINEAR,
                                             .mag_filter = SDL_GPU_FILTER_LINEAR,
                                             .mipmap_mode = SDL_GPU_SAMPLERMIPMAPMODE_LINEAR,
                                             .address_mode_u = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE,
                                             .address_mode_v = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE,
                                             .address_mode_w = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE,
                                             .max_lod = 1000.0f};
    m_sprite_sampler = SDL_CreateGPUSampler(m_gpu, &sampler_info);
}

SDL_GPUTexture* RenderSystem::loadSpriteTexture(const std::string& asset_path)
{
    // 1. Load the PNG into a CPU-side image
    SDL_Surface* loaded = SDL_LoadPNG(Utilities::GetResourcePath(asset_path).c_str());
    if (loaded == nullptr)
    {
        std::cerr << "Couldn't load sprite '" << asset_path << "': " << SDL_GetError() << std::endl;
        return nullptr;
    }

    // 2. Make the bytes R, G, B, A in that order, which is what the R8G8B8A8 texture below expects
    SDL_Surface* image = SDL_ConvertSurface(loaded, SDL_PIXELFORMAT_RGBA32);
    SDL_DestroySurface(loaded);
    if (image == nullptr)
    {
        std::cerr << "Couldn't convert sprite '" << asset_path << "': " << SDL_GetError() << std::endl;
        return nullptr;
    }

    // 3. Premultiply (colour *= alpha), so filtering never blends the black of transparent pixels into the edges
    SDL_PremultiplySurfaceAlpha(image, false);

    // 4. One mip level per halving of the larger side, plus the full-size image itself (2048 -> 12 levels)
    uint32_t mip_levels = 1;
    uint32_t side = (uint32_t)std::max(image->w, image->h);
    while (side > 1)
    {
        side /= 2;
        mip_levels++;
    }

    // SAMPLER: shaders read it. COLOR_TARGET: SDL renders into the smaller levels when generating the mipmaps.
    SDL_GPUTextureCreateInfo tex_info = {.type = SDL_GPU_TEXTURETYPE_2D,
                                         .format = SDL_GPU_TEXTUREFORMAT_R8G8B8A8_UNORM,
                                         .usage = SDL_GPU_TEXTUREUSAGE_SAMPLER | SDL_GPU_TEXTUREUSAGE_COLOR_TARGET,
                                         .width = (Uint32)image->w,
                                         .height = (Uint32)image->h,
                                         .layer_count_or_depth = 1,
                                         .num_levels = mip_levels};
    SDL_GPUTexture* texture = SDL_CreateGPUTexture(m_gpu, &tex_info);
    if (texture == nullptr)
    {
        std::cerr << "Couldn't create a texture for sprite '" << asset_path << "': " << SDL_GetError() << std::endl;
        SDL_DestroySurface(image);
        return nullptr;
    }

    // 5. Copy the full-size image into a transfer buffer, row by row (a surface's rows can be padded: pitch)
    Uint32 row_bytes = (Uint32)image->w * 4;
    SDL_GPUTransferBufferCreateInfo tb_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD,
                                               .size = row_bytes * (Uint32)image->h};
    SDL_GPUTransferBuffer* transfer_buffer = SDL_CreateGPUTransferBuffer(m_gpu, &tb_info);
    Uint8* dst = (Uint8*)SDL_MapGPUTransferBuffer(m_gpu, transfer_buffer, false); // brand-new buffer: nothing to cycle
    Uint8* src = (Uint8*)image->pixels;
    for (int row = 0; row < image->h; row++)
    {
        SDL_memcpy(dst + row * row_bytes, src + row * image->pitch, row_bytes);
    }
    SDL_UnmapGPUTransferBuffer(m_gpu, transfer_buffer);

    // 6. Upload it into mip level 0, then let the GPU shrink level 0 into all the smaller levels
    SDL_GPUCommandBuffer* cmdbuf = SDL_AcquireGPUCommandBuffer(m_gpu);
    SDL_GPUCopyPass* copy_pass = SDL_BeginGPUCopyPass(cmdbuf);
    SDL_GPUTextureTransferInfo src_info = {.transfer_buffer = transfer_buffer,
                                           .offset = 0,
                                           .pixels_per_row = (Uint32)image->w,
                                           .rows_per_layer = (Uint32)image->h};
    SDL_GPUTextureRegion dst_region = {
        .texture = texture, .mip_level = 0, .w = (Uint32)image->w, .h = (Uint32)image->h, .d = 1};
    SDL_UploadToGPUTexture(copy_pass, &src_info, &dst_region, false);
    SDL_EndGPUCopyPass(copy_pass);
    SDL_GenerateMipmapsForGPUTexture(cmdbuf, texture); // must be outside any pass, so after EndGPUCopyPass
    SDL_SubmitGPUCommandBuffer(cmdbuf);

    SDL_ReleaseGPUTransferBuffer(m_gpu, transfer_buffer); // SDL frees it once the upload has finished
    SDL_DestroySurface(image);
    return texture;
}

void RenderSystem::rebuildFontAtlas(float pixel_scale)
{
    SDL_Surface* atlas_surface = fontAtlas.buildAtlas(UIFontRegular, UI_FONT_SIZE, pixel_scale);
    if (!atlas_surface)
    {
        std::cerr << "Failed to bake font atlas" << std::endl;
        return;
    }

    // A frame still on the GPU may be using the old texture; SDL_GPU only frees it once that frame is done
    if (fontAtlasTexture != nullptr)
    {
        SDL_ReleaseGPUTexture(m_gpu, fontAtlasTexture);
    }

    SDL_GPUTextureCreateInfo tex_info = {.type = SDL_GPU_TEXTURETYPE_2D,
                                         .format = SDL_GPU_TEXTUREFORMAT_R8G8B8A8_UNORM,
                                         .usage = SDL_GPU_TEXTUREUSAGE_SAMPLER,
                                         .width = (Uint32)atlas_surface->w,
                                         .height = (Uint32)atlas_surface->h,
                                         .layer_count_or_depth = 1,
                                         .num_levels = 1};
    fontAtlasTexture = SDL_CreateGPUTexture(m_gpu, &tex_info);

    Uint32 pixelDataSize = (Uint32)(atlas_surface->w * atlas_surface->h * 4);
    SDL_GPUTransferBufferCreateInfo tb_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD, .size = pixelDataSize};
    SDL_GPUTransferBuffer* atlasTransferBuffer = SDL_CreateGPUTransferBuffer(m_gpu, &tb_info);

    // Copy row by row in case the surface pitch isn't tightly packed.
    Uint8* dst = (Uint8*)SDL_MapGPUTransferBuffer(m_gpu, atlasTransferBuffer, false);
    Uint8* src = (Uint8*)atlas_surface->pixels;
    Uint32 rowBytes = (Uint32)atlas_surface->w * 4;
    for (int row = 0; row < atlas_surface->h; row++)
    {
        SDL_memcpy(dst + row * rowBytes, src + row * atlas_surface->pitch, rowBytes);
    }
    SDL_UnmapGPUTransferBuffer(m_gpu, atlasTransferBuffer);

    SDL_GPUCommandBuffer* cmdbuf = SDL_AcquireGPUCommandBuffer(m_gpu);
    SDL_GPUCopyPass* copyPass = SDL_BeginGPUCopyPass(cmdbuf);

    SDL_GPUTextureTransferInfo src_info = {.transfer_buffer = atlasTransferBuffer,
                                           .offset = 0,
                                           .pixels_per_row = (Uint32)atlas_surface->w,
                                           .rows_per_layer = (Uint32)atlas_surface->h};
    SDL_GPUTextureRegion dst_region = {.texture = fontAtlasTexture,
                                       .mip_level = 0,
                                       .layer = 0,
                                       .x = 0,
                                       .y = 0,
                                       .z = 0,
                                       .w = (Uint32)atlas_surface->w,
                                       .h = (Uint32)atlas_surface->h,
                                       .d = 1};
    SDL_UploadToGPUTexture(copyPass, &src_info, &dst_region, false);

    SDL_EndGPUCopyPass(copyPass);
    SDL_SubmitGPUCommandBuffer(cmdbuf);

    SDL_ReleaseGPUTransferBuffer(m_gpu, atlasTransferBuffer);
    SDL_DestroySurface(atlas_surface);
}

void RenderSystem::uploadUIVertices(const DynamoEngine::UIRoot& ui, SDL_GPUCommandBuffer* cmdbuf)
{
    m_ui_vertices.clear();
    DynamoEngine::UIGeometryBuilder builder(m_ui_vertices, fontAtlas); // one builder for this frame

    // Every visible element of the scene's UI, back to front
    ui.drawElements(builder);

    // Never copy more than the buffer holds; renderUIElements draws m_ui_vertices.size(), so it stays in sync
    if (m_ui_vertices.size() > MAX_UI_VERTICES)
    {
        printf("Too many UI vertices %zu > %u\n", m_ui_vertices.size(), MAX_UI_VERTICES);
        m_ui_vertices.resize(MAX_UI_VERTICES);
    }

    if (m_ui_vertices.empty())
        return;

    void* map = SDL_MapGPUTransferBuffer(m_gpu, uiTransferBuffer, true);
    SDL_memcpy(map, m_ui_vertices.data(), m_ui_vertices.size() * sizeof(DynamoEngine::UIVertex));
    SDL_UnmapGPUTransferBuffer(m_gpu, uiTransferBuffer);

    SDL_GPUCopyPass* copyPass = SDL_BeginGPUCopyPass(cmdbuf);
    SDL_GPUTransferBufferLocation src = {.transfer_buffer = uiTransferBuffer, .offset = 0};
    SDL_GPUBufferRegion dst = {.buffer = uiVertexBuffer,
                               .offset = 0,
                               .size = (uint32_t)(m_ui_vertices.size() * sizeof(DynamoEngine::UIVertex))};
    SDL_UploadToGPUBuffer(copyPass, &src, &dst, true);
    SDL_EndGPUCopyPass(copyPass);
}

void RenderSystem::renderUIElements(SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf,
                                    const DynamoEngine::UIRoot& ui)
{
    if (m_ui_vertices.empty())
        return;
    SDL_BindGPUGraphicsPipeline(pass, uiPipeline);

    // The UI is built in UI space (m_window size / UI scale); the shader stretches that space over the whole m_window
    const DynamoEngine::Vector2F ui_space_size = ui.uiSpaceSize();
    float ui_space[2] = {ui_space_size.x_val, ui_space_size.y_val};
    SDL_PushGPUVertexUniformData(cmdbuf, 0, ui_space, sizeof(ui_space));

    SDL_GPUBufferBinding vbo = {.buffer = uiVertexBuffer, .offset = 0};
    SDL_BindGPUVertexBuffers(pass, 0, &vbo, 1);
    SDL_GPUTextureSamplerBinding texBinding = {.texture = fontAtlasTexture, .sampler = fontAtlasSampler};
    SDL_BindGPUFragmentSamplers(pass, 0, &texBinding, 1);
    SDL_DrawGPUPrimitives(pass, (uint32_t)m_ui_vertices.size(), 1, 0, 0);
}
void RenderSystem::appendPreviewBodies(std::vector<UnifiedBodyVertex>& vertex_data, UIState& ui_state,
                                       const CameraState& camera_state)
{
    DEPRECATED_InputState& input_state = ui_state.getMutableDEPRECATED_InputState();
    velocityVectorVertices.clear();
    if (!input_state.isPreviewingMacro && !input_state.isPreviewingParticleCluster)
    {
        return;
    }

    // The body that will spawn: the chosen size and mass, at the cursor (or at the drag start, with an arrow to the
    // cursor, when it's being given a velocity)
    GravitationalBody preview_body = {};
    preview_body.mass = input_state.selectedMass;
    preview_body.radius = input_state.selectedRadius;
    if (input_state.isPreviewingWithInitialVelocity)
    {
        preview_body.position = ScreenToWorldCoordinates(input_state.mouseDragStartPosition, camera_state);

        DynamoEngine::Vector2D arrow_end = ScreenToWorldCoordinates(input_state.mouseCurrPosition, camera_state);
        buildVelocityVectorGeometry(preview_body.position, arrow_end);
    }
    else
    {
        preview_body.position = ScreenToWorldCoordinates(input_state.mouseCurrPosition, camera_state);
    }
    preview_body.previousPosition = preview_body.position; // to prevent alpha interpolation artifacts

    if (input_state.isPreviewingMacro)
    {
        preview_body.isPreview = true;
        vertex_data.push_back(preview_body.toUnifiedVertex());
        return;
    }

    // A cluster: show the particles it will be made of, laid out exactly as createParticleCluster will place them
    std::vector<GravitationalBody> preview_particles;
    layOutFragments(preview_body, DEFAULT_FRAGMENT_COUNT, preview_particles);
    for (GravitationalBody& particle : preview_particles)
    {
        particle.isPreview = true;
        vertex_data.push_back(particle.toUnifiedVertex());
    }
}

CameraConstants RenderSystem::buildCameraConstants(const CameraState& camera_state,
                                                   const DynamoEngine::Vector2D& offset)
{
    CameraConstants camera_constants = {};
    camera_constants.screenWidth = camera_state.window_width;
    camera_constants.screenHeight = camera_state.window_height;
    camera_constants.zoom = (float)camera_state.zoom;
    camera_constants.offsetX = (float)offset.x_val;
    camera_constants.offsetY = (float)offset.y_val;
    camera_constants.viewMode = static_cast<uint32_t>(camera_state.visor_view);
    camera_constants.rendering_alpha = camera_state.render_alpha;
    camera_constants._padding1 = 0.0f;
    return camera_constants;
}
// Smooth interpolation color lookup function

// --------- RENDER UTILITY HELPERS --------- //
// Need to rework to fix constant calls (store color in grav body) and also fix beyond max mass opacity.
SDL_Color RenderSystem::getColorForProperty(const GravitationalBody& body)
// NEED TO OPTIMIZE OUT CONSTANT ACCESS
{
    // 1. The "Event Horizon" (Beyond Max Mass)
    double absMass = std::abs(body.mass);
    // std::cout << "body.mass: " << body.mass << std::endl;

    if (absMass > MAX_MASS)
        return SDL_Color{0, 0, 0, 255};

    // 2. The Scaling Factor
    // Power of 0.1 stretches the scale so 10^3 and 10^12 actually look different.
    static const double exponent = 0.1;
    static const double maxScaled = std::pow(MAX_MASS, exponent);
    double t = std::pow(absMass, exponent) / maxScaled;
    t = std::clamp(t, 0.0, 1.0);

    Uint8 r = 0, g = 0, b = 0;

    if (body.mass < 0)
    {
        r = static_cast<Uint8>(100 + (155 * t));
        g = static_cast<Uint8>(200 * std::pow(t, 3)); // Adds "heat" (yellowish) at high mass
        b = static_cast<Uint8>(200 * std::pow(t, 5)); // Becomes white at the very limit
    }
    else
    {
        r = static_cast<Uint8>(200 * std::pow(t, 5));
        g = static_cast<Uint8>(200 * std::pow(t, 3));
        b = static_cast<Uint8>(100 + (155 * t));
    }

    Uint8 opacity = (body.isMacroGhost || !body.isCollidable) ? 175 : 255;
    return SDL_Color{r, g, b, opacity};
}
bool RenderSystem::createUnifiedBodyBuffers(uint32_t capacity)
{
    // SDL frees these once the GPU has finished with them, so releasing while last frame's draw may still be
    // running is safe. The contents don't need copying: every frame uploads all bodies from scratch anyway.
    if (unifiedBodyVertexBuffer != nullptr)
    {
        SDL_ReleaseGPUBuffer(m_gpu, unifiedBodyVertexBuffer);
    }
    if (unifiedBodyTransferBuffer != nullptr)
    {
        SDL_ReleaseGPUTransferBuffer(m_gpu, unifiedBodyTransferBuffer);
    }

    const uint32_t size_in_bytes = capacity * (uint32_t)sizeof(UnifiedBodyVertex);

    SDL_GPUBufferCreateInfo vb_info = {.usage = SDL_GPU_BUFFERUSAGE_VERTEX, .size = size_in_bytes};
    unifiedBodyVertexBuffer = SDL_CreateGPUBuffer(m_gpu, &vb_info);

    SDL_GPUTransferBufferCreateInfo tb_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD, .size = size_in_bytes};
    unifiedBodyTransferBuffer = SDL_CreateGPUTransferBuffer(m_gpu, &tb_info);

    if (unifiedBodyVertexBuffer == nullptr || unifiedBodyTransferBuffer == nullptr)
    {
        printf("Couldn't create body buffers for %u bodies: %s\n", capacity, SDL_GetError());
        m_unified_body_capacity = 0; // so the next frame tries again
        return false;
    }

    m_unified_body_capacity = capacity;
    return true;
}

void RenderSystem::createUnifiedBodyGPUBufferAndPipeline()
{
    // Reserved at full size up front, so normal play never reallocates
    createUnifiedBodyBuffers(INITIAL_UNIFIED_BODY_CAPACITY);

    SDL_GPUShader* vert_shader = LoadShader(m_gpu, "Shaders/UnifiedGravBody.vert", 0, 1);
    SDL_GPUShader* frag_shader = LoadShader(m_gpu, "Shaders/UnifiedGravBody.frag", 0, 0);

    // Define Pipeline State
    SDL_GPUVertexAttribute vertex_attributes[8];

    // [0] Position, float2
    vertex_attributes[0].location = 0;
    vertex_attributes[0].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2;
    vertex_attributes[0].offset = offsetof(UnifiedBodyVertex, x);
    vertex_attributes[0].buffer_slot = 0;

    // [1] Prev Position, float2
    vertex_attributes[1].location = 1;
    vertex_attributes[1].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2;
    vertex_attributes[1].offset = offsetof(UnifiedBodyVertex, prevX);
    vertex_attributes[1].buffer_slot = 0;

    // [2] Radius, float1
    vertex_attributes[2].location = 2;
    vertex_attributes[2].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT;
    vertex_attributes[2].offset = offsetof(UnifiedBodyVertex, radius);
    vertex_attributes[2].buffer_slot = 0;

    // [3] log_mass, float1
    vertex_attributes[3].location = 3;
    vertex_attributes[3].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT;
    vertex_attributes[3].offset = offsetof(UnifiedBodyVertex, logMass);
    vertex_attributes[3].buffer_slot = 0;

    // [4] temperature, float1
    vertex_attributes[4].location = 4;
    vertex_attributes[4].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT;
    vertex_attributes[4].offset = offsetof(UnifiedBodyVertex, temperature);
    vertex_attributes[4].buffer_slot = 0;

    // [5] charge, float1
    vertex_attributes[5].location = 5;
    vertex_attributes[5].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT;
    vertex_attributes[5].offset = offsetof(UnifiedBodyVertex, charge);
    vertex_attributes[5].buffer_slot = 0;

    // [6] flags, uint1
    vertex_attributes[6].location = 6;
    vertex_attributes[6].format = SDL_GPU_VERTEXELEMENTFORMAT_UINT;
    vertex_attributes[6].offset = offsetof(UnifiedBodyVertex, flags);
    vertex_attributes[6].buffer_slot = 0;

    // [7] seed, uint1
    vertex_attributes[7].location = 7;
    vertex_attributes[7].format = SDL_GPU_VERTEXELEMENTFORMAT_UINT;
    vertex_attributes[7].offset = offsetof(UnifiedBodyVertex, seed);
    vertex_attributes[7].buffer_slot = 0;

    SDL_GPUGraphicsPipelineCreateInfo pipeline_info = {};
    pipeline_info.target_info.num_color_targets = 1;

    SDL_GPUColorTargetDescription color_target = {};
    color_target.format = SDL_GetGPUSwapchainTextureFormat(m_gpu, m_window);
    color_target.blend_state.enable_blend = true;

    color_target.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
    color_target.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
    color_target.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
    color_target.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;

    // MUST SET THESE EXPLICITLY:
    color_target.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
    color_target.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;

    pipeline_info.target_info.color_target_descriptions = &color_target;
    pipeline_info.vertex_shader = vert_shader;
    pipeline_info.fragment_shader = frag_shader;

    // Change to TRIANGLESTRIP ``TODO: CHECK?`` for your Quads!
    pipeline_info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;

    pipeline_info.vertex_input_state.vertex_attributes = vertex_attributes;
    pipeline_info.vertex_input_state.num_vertex_attributes = 8;

    SDL_GPUVertexBufferDescription vbo_desc = {
        .slot = 0, .pitch = sizeof(UnifiedBodyVertex), .input_rate = SDL_GPU_VERTEXINPUTRATE_INSTANCE};
    pipeline_info.vertex_input_state.vertex_buffer_descriptions = &vbo_desc;
    pipeline_info.vertex_input_state.num_vertex_buffers = 1;

    unifiedBodyPipeline = SDL_CreateGPUGraphicsPipeline(m_gpu, &pipeline_info);

    // 7. Cleanup temp shader handles
    SDL_ReleaseGPUShader(m_gpu, vert_shader);
    SDL_ReleaseGPUShader(m_gpu, frag_shader);
}

void RenderSystem::createTwinklingStarGPUBufferAndPipeline()
{

    SDL_GPUBufferCreateInfo vb_info = {.usage = SDL_GPU_BUFFERUSAGE_VERTEX,
                                       .size = STAR_NUM * sizeof(TwinklingStarVertex)};
    twinklingStarVertexBuffer = SDL_CreateGPUBuffer(m_gpu, &vb_info);

    SDL_GPUTransferBufferCreateInfo tb_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD,
                                               .size = STAR_NUM * sizeof(TwinklingStarVertex)};

    twinklingStarTransferBuffer = SDL_CreateGPUTransferBuffer(m_gpu, &tb_info);

    SDL_GPUShader* vert_shader = LoadShader(m_gpu, "Shaders/TwinklingStar.vert", 0, 1);
    SDL_GPUShader* frag_shader = LoadShader(m_gpu, "Shaders/TwinklingStar.frag", 0, 1);

    SDL_GPUVertexAttribute vertex_attributes[5];

    // Position
    vertex_attributes[0].location = 0;
    vertex_attributes[0].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2;
    vertex_attributes[0].offset = offsetof(TwinklingStarVertex, x);
    vertex_attributes[0].buffer_slot = 0; // Need to check?

    vertex_attributes[1].location = 1;
    vertex_attributes[1].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT;
    vertex_attributes[1].offset = offsetof(TwinklingStarVertex, radius);
    vertex_attributes[1].buffer_slot = 0;

    vertex_attributes[2].location = 2;
    vertex_attributes[2].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT;
    vertex_attributes[2].offset = offsetof(TwinklingStarVertex, alpha);
    vertex_attributes[2].buffer_slot = 0;

    vertex_attributes[3].location = 3;
    vertex_attributes[3].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT;
    vertex_attributes[3].offset = offsetof(TwinklingStarVertex, twinkleSpeed);
    vertex_attributes[3].buffer_slot = 0;

    vertex_attributes[4].location = 4;
    vertex_attributes[4].format = SDL_GPU_VERTEXELEMENTFORMAT_UINT;
    vertex_attributes[4].offset = offsetof(TwinklingStarVertex, seed);
    vertex_attributes[4].buffer_slot = 0;

    SDL_GPUGraphicsPipelineCreateInfo pipeline_info = {};
    pipeline_info.target_info.num_color_targets = 1;

    SDL_GPUColorTargetDescription color_target = {};
    color_target.format = SDL_GetGPUSwapchainTextureFormat(m_gpu, m_window);
    color_target.blend_state.enable_blend = true;

    color_target.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
    color_target.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
    color_target.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
    color_target.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;

    // MUST SET THESE EXPLICITLY:
    color_target.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
    color_target.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;

    pipeline_info.target_info.color_target_descriptions = &color_target;
    pipeline_info.vertex_shader = vert_shader;
    pipeline_info.fragment_shader = frag_shader;

    // Change to TRIANGLESTRIP ``TODO: CHECK?`` for your Quads!
    pipeline_info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;

    pipeline_info.vertex_input_state.vertex_attributes = vertex_attributes;
    pipeline_info.vertex_input_state.num_vertex_attributes = 5;

    SDL_GPUVertexBufferDescription vbo_desc = {
        .slot = 0, .pitch = sizeof(TwinklingStarVertex), .input_rate = SDL_GPU_VERTEXINPUTRATE_INSTANCE};
    pipeline_info.vertex_input_state.vertex_buffer_descriptions = &vbo_desc;
    pipeline_info.vertex_input_state.num_vertex_buffers = 1;

    twinklingStarPipeline = SDL_CreateGPUGraphicsPipeline(m_gpu, &pipeline_info);

    SDL_ReleaseGPUShader(m_gpu, vert_shader);
    SDL_ReleaseGPUShader(m_gpu, frag_shader);
}

void RenderSystem::createTwinklingStarField(float fieldMaxWidth, float fieldMaxHeight)
{
    twinklingStarVertices.clear();
    twinklingStarVertices.reserve(STAR_NUM);

    std::mt19937 rng(12345);

    double field_width = fieldMaxWidth / MIN_ZOOM;
    double field_height = fieldMaxHeight / MIN_ZOOM;
    double marginX = (field_width - SCREEN_WIDTH) / 2.0;
    double marginY = (field_height - SCREEN_HEIGHT) / 2.0;
    std::uniform_real_distribution<float> posX((float)(-marginX), (float)(field_width + marginX));
    std::uniform_real_distribution<float> posY((float)(-marginY), (float)(field_height + marginY));

    std::uniform_real_distribution<float> radiusDist(8.0f, 24.0f);

    std::uniform_real_distribution<float> alphaDist(0.3f, 1.0f);

    std::uniform_real_distribution<float> twinkleDist(0.25f, 2.5f);

    for (uint32_t i = 0; i < STAR_NUM; i++)
    {
        TwinklingStarVertex star;
        star.x = posX(rng);
        star.y = posY(rng);
        star.radius = radiusDist(rng);
        star.alpha = alphaDist(rng);
        star.twinkleSpeed = twinkleDist(rng);
        star.seed = rng();
        twinklingStarVertices.push_back(star);
    }
}

void RenderSystem::uploadTwinklingStarField(SDL_GPUCommandBuffer* cmdbuf)
{

    void* map = SDL_MapGPUTransferBuffer(m_gpu, twinklingStarTransferBuffer, true);

    SDL_memcpy(map, twinklingStarVertices.data(), twinklingStarVertices.size() * sizeof(TwinklingStarVertex));

    SDL_UnmapGPUTransferBuffer(m_gpu, twinklingStarTransferBuffer);

    SDL_GPUCopyPass* copyPass = SDL_BeginGPUCopyPass(cmdbuf);

    SDL_GPUTransferBufferLocation src = {.transfer_buffer = twinklingStarTransferBuffer, .offset = 0};

    SDL_GPUBufferRegion dst = {.buffer = twinklingStarVertexBuffer,
                               .offset = 0,
                               .size = (uint32_t)(twinklingStarVertices.size() * sizeof(TwinklingStarVertex))};

    SDL_UploadToGPUBuffer(copyPass, &src, &dst, true);

    SDL_EndGPUCopyPass(copyPass);
}

void RenderSystem::renderTwinklingStarField(SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf,
                                            const CameraState& camera_state)
{
    SDL_BindGPUGraphicsPipeline(pass, twinklingStarPipeline);

    CameraConstants camera_constants = buildCameraConstants(camera_state, camera_state.twinkling_star_offset);

    SDL_PushGPUVertexUniformData(cmdbuf, 0, &camera_constants, sizeof(camera_constants));

    float elapsedSeconds = (float)SDL_GetTicks() / 1000.0f;
    SDL_PushGPUFragmentUniformData(cmdbuf, 0, &elapsedSeconds, sizeof(elapsedSeconds));
    SDL_GPUBufferBinding vbo = {.buffer = twinklingStarVertexBuffer, .offset = 0};
    SDL_BindGPUVertexBuffers(pass, 0, &vbo, 1);
    SDL_DrawGPUPrimitives(pass, 6, STAR_NUM, 0, 0);
}

void RenderSystem::createVelocityVectorGPUBufferAndPipeline()
{
    SDL_GPUBufferCreateInfo vb_info = {.usage = SDL_GPU_BUFFERUSAGE_VERTEX,
                                       .size = MAX_VELOCITY_VECTOR_VERTICES * sizeof(VelocityVectorVertex)};
    velocityVectorVertexBuffer = SDL_CreateGPUBuffer(m_gpu, &vb_info);

    SDL_GPUTransferBufferCreateInfo tb_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD,
                                               .size = MAX_VELOCITY_VECTOR_VERTICES * sizeof(VelocityVectorVertex)};
    velocityVectorTransferBuffer = SDL_CreateGPUTransferBuffer(m_gpu, &tb_info);

    SDL_GPUShader* vert_shader = LoadShader(m_gpu, "Shaders/VelocityVector.vert", 0, 1);
    SDL_GPUShader* frag_shader = LoadShader(m_gpu, "Shaders/VelocityVector.frag", 0, 0);

    SDL_GPUVertexAttribute vertex_attributes[2];
    vertex_attributes[0].location = 0;
    vertex_attributes[0].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2;
    vertex_attributes[0].offset = offsetof(VelocityVectorVertex, x);
    vertex_attributes[0].buffer_slot = 0;

    vertex_attributes[1].location = 1;
    vertex_attributes[1].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT4;
    vertex_attributes[1].offset = offsetof(VelocityVectorVertex, r);
    vertex_attributes[1].buffer_slot = 0;

    SDL_GPUGraphicsPipelineCreateInfo pipeline_info = {};
    pipeline_info.target_info.num_color_targets = 1;

    SDL_GPUColorTargetDescription color_target = {};
    color_target.format = SDL_GetGPUSwapchainTextureFormat(m_gpu, m_window);
    color_target.blend_state.enable_blend = true;
    color_target.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
    color_target.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
    color_target.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
    color_target.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
    color_target.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
    color_target.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;

    pipeline_info.target_info.color_target_descriptions = &color_target;
    pipeline_info.vertex_shader = vert_shader;
    pipeline_info.fragment_shader = frag_shader;
    pipeline_info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;

    pipeline_info.vertex_input_state.vertex_attributes = vertex_attributes;
    pipeline_info.vertex_input_state.num_vertex_attributes = 2;

    SDL_GPUVertexBufferDescription vbo_desc = {
        .slot = 0, .pitch = sizeof(VelocityVectorVertex), .input_rate = SDL_GPU_VERTEXINPUTRATE_VERTEX};
    pipeline_info.vertex_input_state.vertex_buffer_descriptions = &vbo_desc;
    pipeline_info.vertex_input_state.num_vertex_buffers = 1;

    velocityVectorPipeline = SDL_CreateGPUGraphicsPipeline(m_gpu, &pipeline_info);

    SDL_ReleaseGPUShader(m_gpu, vert_shader);
    SDL_ReleaseGPUShader(m_gpu, frag_shader);
}

void RenderSystem::createStarshipGPUBufferAndPipeline()
{
    SDL_GPUBufferCreateInfo vb_info = {.usage = SDL_GPU_BUFFERUSAGE_VERTEX,
                                       .size = MAX_STARSHIP_VERTICES * sizeof(StarshipVertex)};

    starshipVertexBuffer = SDL_CreateGPUBuffer(m_gpu, &vb_info);

    SDL_GPUTransferBufferCreateInfo tb_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD,
                                               .size = MAX_STARSHIP_VERTICES * sizeof(StarshipVertex)};
    starshipTransferBuffer = SDL_CreateGPUTransferBuffer(m_gpu, &tb_info);

    SDL_GPUShader* vert_shader = LoadShader(m_gpu, "Shaders/Starship.vert", 0, 1);
    SDL_GPUShader* frag_shader = LoadShader(m_gpu, "Shaders/Starship.frag", 1, 0); // 1 sampler

    SDL_GPUVertexAttribute vertex_attributes[5];

    // Position float2
    vertex_attributes[0].location = 0;
    vertex_attributes[0].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2;
    vertex_attributes[0].offset = offsetof(StarshipVertex, x);
    vertex_attributes[0].buffer_slot = 0;

    // Previous position float2
    vertex_attributes[1].location = 1;
    vertex_attributes[1].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2;
    vertex_attributes[1].offset = offsetof(StarshipVertex, prevX);
    vertex_attributes[1].buffer_slot = 0;

    // square size, float
    vertex_attributes[2].location = 2;
    vertex_attributes[2].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT;
    vertex_attributes[2].offset = offsetof(StarshipVertex, square_size);
    vertex_attributes[2].buffer_slot = 0;

    // texture coordinates (u, v) float2
    vertex_attributes[3].location = 3;
    vertex_attributes[3].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2;
    vertex_attributes[3].offset = offsetof(StarshipVertex, u);
    vertex_attributes[3].buffer_slot = 0;

    // Color float 4
    vertex_attributes[4].location = 4;
    vertex_attributes[4].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT4;
    vertex_attributes[4].offset = offsetof(StarshipVertex, r);
    vertex_attributes[4].buffer_slot = 0;

    SDL_GPUGraphicsPipelineCreateInfo pipeline_info = {};
    pipeline_info.target_info.num_color_targets = 1;

    SDL_GPUColorTargetDescription color_target = {};
    color_target.format = SDL_GetGPUSwapchainTextureFormat(m_gpu, m_window);
    color_target.blend_state.enable_blend = true;

    // Premultiplied alpha (the sprite's colours were already multiplied by alpha when it was loaded), so the
    // source is added as-is (ONE) instead of being multiplied by its alpha a second time
    color_target.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color_target.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
    color_target.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color_target.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;

    // MUST SET THESE EXPLICITLY:
    color_target.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
    color_target.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;

    pipeline_info.target_info.color_target_descriptions = &color_target;
    pipeline_info.vertex_shader = vert_shader;
    pipeline_info.fragment_shader = frag_shader;

    // // Change to TRIANGLESTRIP ``TODO: CHECK?`` for quads
    pipeline_info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;

    pipeline_info.vertex_input_state.vertex_attributes = vertex_attributes;
    pipeline_info.vertex_input_state.num_vertex_attributes = 5;

    SDL_GPUVertexBufferDescription vbo_desc = {
        .slot = 0, .pitch = sizeof(StarshipVertex), .input_rate = SDL_GPU_VERTEXINPUTRATE_VERTEX};
    pipeline_info.vertex_input_state.vertex_buffer_descriptions = &vbo_desc;
    pipeline_info.vertex_input_state.num_vertex_buffers = 1;

    starshipPipeline = SDL_CreateGPUGraphicsPipeline(m_gpu, &pipeline_info);
    if (!starshipPipeline)
        std::cerr << "Starship pipeline creation failed: " << SDL_GetError() << std::endl;
    SDL_ReleaseGPUShader(m_gpu, vert_shader);
    SDL_ReleaseGPUShader(m_gpu, frag_shader);
}

void RenderSystem::buildVelocityVectorGeometry(DynamoEngine::Vector2D lineStart, DynamoEngine::Vector2D lineEnd)
{
    velocityVectorVertices.clear();

    float dx = static_cast<float>(lineEnd.x_val - lineStart.x_val);
    float dy = static_cast<float>(lineEnd.y_val - lineStart.y_val);
    float length = static_cast<float>((lineEnd - lineStart).magnitude());
    if (length <= DynamoEngine::EPSILON)
    {
        return;
    }
    dx /= length;
    dy /= length;
    float px = -dy;
    float py = dx;

    const float thickness = 12.0f;
    const float arrow_length = 30.0f;
    const float arrow_width = 36.0f;
    float half_T = thickness / 2.0f;

    // Shorten the shaft so the arrowhead has room at the tip
    DynamoEngine::Vector2D line_end = {lineEnd.x_val - dx * arrow_length, lineEnd.y_val - dy * arrow_length};

    float r = 1.0f, g = 1.0f, b = 1.0f, a = 1.0f; // white, matching the original

    // Shaft quad, expanded into 2 raw triangles
    float v0x = float(lineStart.x_val) + px * half_T, v0y = float(lineStart.y_val) + py * half_T;
    float v1x = float(lineStart.x_val) - px * half_T, v1y = float(lineStart.y_val) - py * half_T;
    float v2x = float(line_end.x_val) - px * half_T, v2y = float(line_end.y_val) - py * half_T;
    float v3x = float(line_end.x_val) + px * half_T, v3y = float(line_end.y_val) + py * half_T;

    velocityVectorVertices.push_back({v0x, v0y, r, g, b, a});
    velocityVectorVertices.push_back({v1x, v1y, r, g, b, a});
    velocityVectorVertices.push_back({v2x, v2y, r, g, b, a});
    velocityVectorVertices.push_back({v0x, v0y, r, g, b, a});
    velocityVectorVertices.push_back({v2x, v2y, r, g, b, a});
    velocityVectorVertices.push_back({v3x, v3y, r, g, b, a});

    // Arrowhead triangle
    velocityVectorVertices.push_back({float(lineEnd.x_val), float(lineEnd.y_val), r, g, b, a});
    velocityVectorVertices.push_back({float(line_end.x_val) + px * (arrow_width / 2.0f),
                                      float(line_end.y_val) + py * (arrow_width / 2.0f), r, g, b, a});
    velocityVectorVertices.push_back({float(line_end.x_val) - px * (arrow_width / 2.0f),
                                      float(line_end.y_val) - py * (arrow_width / 2.0f), r, g, b, a});
}

void RenderSystem::uploadVelocityVectorVertices(SDL_GPUCommandBuffer* cmdbuf)
{
    if (velocityVectorVertices.empty())
        return;

    // Never copy more than the buffer holds; renderVelocityVectors draws velocityVectorVertices.size()
    if (velocityVectorVertices.size() > MAX_VELOCITY_VECTOR_VERTICES)
    {
        printf("Too many velocity vector vertices %zu > %d\n", velocityVectorVertices.size(),
               MAX_VELOCITY_VECTOR_VERTICES);
        velocityVectorVertices.resize(MAX_VELOCITY_VECTOR_VERTICES);
    }

    void* map = SDL_MapGPUTransferBuffer(m_gpu, velocityVectorTransferBuffer, true);
    SDL_memcpy(map, velocityVectorVertices.data(), velocityVectorVertices.size() * sizeof(VelocityVectorVertex));
    SDL_UnmapGPUTransferBuffer(m_gpu, velocityVectorTransferBuffer);

    SDL_GPUCopyPass* copyPass = SDL_BeginGPUCopyPass(cmdbuf);
    SDL_GPUTransferBufferLocation src = {.transfer_buffer = velocityVectorTransferBuffer, .offset = 0};
    SDL_GPUBufferRegion dst = {.buffer = velocityVectorVertexBuffer,
                               .offset = 0,
                               .size = (uint32_t)(velocityVectorVertices.size() * sizeof(VelocityVectorVertex))};
    SDL_UploadToGPUBuffer(copyPass, &src, &dst, true);
    SDL_EndGPUCopyPass(copyPass);
}

void RenderSystem::renderVelocityVectors(SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf,
                                         const CameraState& camera_state)
{
    if (velocityVectorVertices.empty())
        return;

    SDL_BindGPUGraphicsPipeline(pass, velocityVectorPipeline);

    CameraConstants camera_constants = buildCameraConstants(camera_state, camera_state.offset);
    SDL_PushGPUVertexUniformData(cmdbuf, 0, &camera_constants, sizeof(camera_constants));

    SDL_GPUBufferBinding vbo = {.buffer = velocityVectorVertexBuffer, .offset = 0};
    SDL_BindGPUVertexBuffers(pass, 0, &vbo, 1);
    SDL_DrawGPUPrimitives(pass, (uint32_t)velocityVectorVertices.size(), 1, 0, 0);
}

void RenderSystem::createDebugLineGPUBufferAndPipeline()
{
    SDL_GPUBufferCreateInfo vb_info = {.usage = SDL_GPU_BUFFERUSAGE_VERTEX,
                                       .size = MAX_DEBUG_LINE_VERTICES * sizeof(DebugLineVertex)};
    m_debug_line_vertex_buffer = SDL_CreateGPUBuffer(m_gpu, &vb_info);

    SDL_GPUTransferBufferCreateInfo tb_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD,
                                               .size = MAX_DEBUG_LINE_VERTICES * sizeof(DebugLineVertex)};
    m_debug_line_transfer_buffer = SDL_CreateGPUTransferBuffer(m_gpu, &tb_info);

    SDL_GPUShader* vert_shader = LoadShader(m_gpu, "Shaders/DebugLine.vert", 0, 1);
    SDL_GPUShader* frag_shader = LoadShader(m_gpu, "Shaders/DebugLine.frag", 0, 0);

    // The three fields of DebugLineVertex, in the order the vertex shader lists its inputs
    SDL_GPUVertexAttribute vertex_attributes[3];

    // Position float2
    vertex_attributes[0].location = 0;
    vertex_attributes[0].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2;
    vertex_attributes[0].offset = offsetof(DebugLineVertex, x);
    vertex_attributes[0].buffer_slot = 0;

    // Previous-tick position float2
    vertex_attributes[1].location = 1;
    vertex_attributes[1].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2;
    vertex_attributes[1].offset = offsetof(DebugLineVertex, prevX);
    vertex_attributes[1].buffer_slot = 0;

    // Color float4
    vertex_attributes[2].location = 2;
    vertex_attributes[2].format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT4;
    vertex_attributes[2].offset = offsetof(DebugLineVertex, r);
    vertex_attributes[2].buffer_slot = 0;

    SDL_GPUColorTargetDescription color_target = {};
    color_target.format = SDL_GetGPUSwapchainTextureFormat(m_gpu, m_window);
    color_target.blend_state.enable_blend = true;
    color_target.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
    color_target.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
    color_target.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
    color_target.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
    color_target.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
    color_target.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;

    SDL_GPUVertexBufferDescription vbo_desc = {
        .slot = 0, .pitch = sizeof(DebugLineVertex), .input_rate = SDL_GPU_VERTEXINPUTRATE_VERTEX};

    SDL_GPUGraphicsPipelineCreateInfo pipeline_info = {};
    pipeline_info.target_info.num_color_targets = 1;
    pipeline_info.target_info.color_target_descriptions = &color_target;
    pipeline_info.vertex_shader = vert_shader;
    pipeline_info.fragment_shader = frag_shader;
    pipeline_info.primitive_type = SDL_GPU_PRIMITIVETYPE_LINELIST; // every 2 vertices = one separate line
    pipeline_info.vertex_input_state.vertex_attributes = vertex_attributes;
    pipeline_info.vertex_input_state.num_vertex_attributes = 3;
    pipeline_info.vertex_input_state.vertex_buffer_descriptions = &vbo_desc;
    pipeline_info.vertex_input_state.num_vertex_buffers = 1;

    m_debug_line_pipeline = SDL_CreateGPUGraphicsPipeline(m_gpu, &pipeline_info);
    if (!m_debug_line_pipeline)
        std::cerr << "Debug line pipeline creation failed: " << SDL_GetError() << std::endl;
    SDL_ReleaseGPUShader(m_gpu, vert_shader);
    SDL_ReleaseGPUShader(m_gpu, frag_shader);
}

void RenderSystem::buildDebugLines(const GameState& game_state, UIState& ui_state)
{
    m_debug_line_vertices.clear();
    if (!ui_state.getRenderDebug())
    {
        return; // overlay is off (F3): nothing to draw
    }

    // Planets: their collision circle is simply their radius (yellow)
    const SDL_FColor PLANET_COLOR = {1.0f, 0.9f, 0.1f, 1.0f};
    for (const GravitationalBody& body : game_state.getUniverseState().macro_bodies)
    {
        addDebugCircle(body.position, body.previousPosition, body.radius, PLANET_COLOR);
    }

    // The ship: its hitbox polygon (green). One tick ago it was the same shape moved back by this tick's movement
    // (the rotation isn't interpolated, exactly like the sprite)
    const SDL_FColor SHIP_COLOR = {0.2f, 1.0f, 0.3f, 1.0f};
    const Starship& ship = game_state.getUniverseState().ship;
    std::vector<DynamoEngine::Vector2D> corners = ship.collisionPolygon();
    DynamoEngine::Vector2D movement = ship.movementThisTick();
    for (size_t i = 0; i < corners.size(); ++i)
    {
        const DynamoEngine::Vector2D& start = corners[i];
        const DynamoEngine::Vector2D& end = corners[(i + 1) % corners.size()]; // the last corner joins the first
        addDebugLine(start, start - movement, end, end - movement, SHIP_COLOR);
    }
}

void RenderSystem::addDebugLine(const DynamoEngine::Vector2D& start, const DynamoEngine::Vector2D& prev_start,
                                const DynamoEngine::Vector2D& end, const DynamoEngine::Vector2D& prev_end,
                                SDL_FColor color)
{
    m_debug_line_vertices.push_back({(float)start.x_val, (float)start.y_val, (float)prev_start.x_val,
                                     (float)prev_start.y_val, color.r, color.g, color.b, color.a});
    m_debug_line_vertices.push_back({(float)end.x_val, (float)end.y_val, (float)prev_end.x_val, (float)prev_end.y_val,
                                     color.r, color.g, color.b, color.a});
}

void RenderSystem::addDebugCircle(const DynamoEngine::Vector2D& center, const DynamoEngine::Vector2D& prev_center,
                                  double radius, SDL_FColor color)
{
    // Drawn as DEBUG_CIRCLE_SEGMENTS straight lines around the edge: plenty to look round at any zoom
    const int DEBUG_CIRCLE_SEGMENTS = 32;
    for (int i = 0; i < DEBUG_CIRCLE_SEGMENTS; ++i)
    {
        double start_angle = DynamoEngine::TWO_PI * i / DEBUG_CIRCLE_SEGMENTS;
        double end_angle = DynamoEngine::TWO_PI * (i + 1) / DEBUG_CIRCLE_SEGMENTS;
        DynamoEngine::Vector2D start_offset = {std::cos(start_angle) * radius, std::sin(start_angle) * radius};
        DynamoEngine::Vector2D end_offset = {std::cos(end_angle) * radius, std::sin(end_angle) * radius};
        addDebugLine(center + start_offset, prev_center + start_offset, center + end_offset, prev_center + end_offset,
                     color);
    }
}

void RenderSystem::uploadDebugLines(SDL_GPUCommandBuffer* cmdbuf)
{
    if (m_debug_line_vertices.empty())
        return;

    // Never copy more than the buffer holds; renderDebugLines draws m_debug_line_vertices.size(), so it stays in sync
    if (m_debug_line_vertices.size() > MAX_DEBUG_LINE_VERTICES)
    {
        printf("Too many debug line vertices %zu > %u\n", m_debug_line_vertices.size(), MAX_DEBUG_LINE_VERTICES);
        m_debug_line_vertices.resize(MAX_DEBUG_LINE_VERTICES);
    }

    // cycle = true: last frame's commands may still be reading these buffers; SDL hands us a free copy instead
    void* map = SDL_MapGPUTransferBuffer(m_gpu, m_debug_line_transfer_buffer, true);
    SDL_memcpy(map, m_debug_line_vertices.data(), m_debug_line_vertices.size() * sizeof(DebugLineVertex));
    SDL_UnmapGPUTransferBuffer(m_gpu, m_debug_line_transfer_buffer);

    SDL_GPUCopyPass* copy_pass = SDL_BeginGPUCopyPass(cmdbuf);
    SDL_GPUTransferBufferLocation src = {.transfer_buffer = m_debug_line_transfer_buffer, .offset = 0};
    SDL_GPUBufferRegion dst = {.buffer = m_debug_line_vertex_buffer,
                               .offset = 0,
                               .size = (uint32_t)(m_debug_line_vertices.size() * sizeof(DebugLineVertex))};
    SDL_UploadToGPUBuffer(copy_pass, &src, &dst, true);
    SDL_EndGPUCopyPass(copy_pass);
}

void RenderSystem::renderDebugLines(SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf,
                                    const CameraState& camera_state)
{
    if (m_debug_line_vertices.empty())
        return; // overlay off, or nothing to outline

    SDL_BindGPUGraphicsPipeline(pass, m_debug_line_pipeline);
    CameraConstants camera_constants = buildCameraConstants(camera_state, camera_state.offset);
    SDL_PushGPUVertexUniformData(cmdbuf, 0, &camera_constants, sizeof(camera_constants));

    SDL_GPUBufferBinding vbo = {.buffer = m_debug_line_vertex_buffer, .offset = 0};
    SDL_BindGPUVertexBuffers(pass, 0, &vbo, 1);
    SDL_DrawGPUPrimitives(pass, (uint32_t)m_debug_line_vertices.size(), 1, 0, 0);
}

// --------- TRAJECTORY LINES (PLANNING MODE) --------- //

// The look of the lines (tune freely). Each new piece fades in from FADE_START to LINE alpha (in the shader); while
// tracing, the newest few segments are drawn brighter on top, so the front of the line is easy to see.
static constexpr SDL_FColor TRAJECTORY_MACRO_COLOR = {0.55f, 0.75f, 1.0f, 1.0f};   // pale blue
static constexpr SDL_FColor TRAJECTORY_SHIP_COLOR = {0.4f, 1.0f, 0.5f, 1.0f};      // green, like the ship's hitbox
static constexpr SDL_FColor TRAJECTORY_STATIC_COLOR = {1.0f, 0.65f, 0.2f, 1.0f};   // amber: a static body, no forces
static constexpr SDL_FColor TRAJECTORY_PREVIEW_COLOR = {0.95f, 0.95f, 0.95f, 1.0f}; // the planet being placed
static constexpr SDL_FColor TRAJECTORY_TIP_COLOR = {1.0f, 1.0f, 1.0f, 1.0f};       // the bright front while tracing
static constexpr float TRAJECTORY_LINE_ALPHA = 0.7f;       // a settled line
static constexpr float TRAJECTORY_FADE_START_ALPHA = 0.1f; // a brand-new piece
static constexpr float TRAJECTORY_FADE_IN_SECONDS = 0.75f; // real seconds from new to settled
static constexpr int TRAJECTORY_TIP_SEGMENTS = 3;          // how much of the front is drawn bright (at most 8: see
                                                           // the MAX_TRAJECTORY_FRAME_VERTICES check)
static constexpr float TRAJECTORY_TIME_MARK_HALF_LENGTH_PX = 5.0f; // each time mark sticks out this far on screen

// The TrajectoryLine vertex shader's second uniform buffer: must match its TrajectoryConstants, field for field
struct TrajectoryShaderConstants
{
    float clock;
    float fade_in_seconds;
    float fade_start_alpha;
    float padding;
};

static TrajectoryLineVertex trajectoryVertex(const DynamoEngine::Vector2D& position,
                                             const DynamoEngine::Vector2D& previous, SDL_FColor color, float alpha,
                                             float drawn_at, bool fades_in)
{
    return {(float)position.x_val, (float)position.y_val, (float)previous.x_val, (float)previous.y_val,
            color.r,               color.g,               color.b,               alpha,
            drawn_at,              fades_in ? 1.0f : 0.0f};
}

void RenderSystem::createTrajectoryLineBuffersAndPipeline()
{
    SDL_GPUBufferCreateInfo line_vb_info = {.usage = SDL_GPU_BUFFERUSAGE_VERTEX,
                                            .size = MAX_TRAJECTORY_LINE_VERTICES * sizeof(TrajectoryLineVertex)};
    m_trajectory_line_vertex_buffer = SDL_CreateGPUBuffer(m_gpu, &line_vb_info);
    // As big as the whole line buffer, so one upload can never be too large for it (uploads are normally small)
    SDL_GPUTransferBufferCreateInfo line_tb_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD,
                                                    .size = MAX_TRAJECTORY_LINE_VERTICES * sizeof(TrajectoryLineVertex)};
    m_trajectory_line_transfer_buffer = SDL_CreateGPUTransferBuffer(m_gpu, &line_tb_info);

    SDL_GPUBufferCreateInfo frame_vb_info = {.usage = SDL_GPU_BUFFERUSAGE_VERTEX,
                                             .size = MAX_TRAJECTORY_FRAME_VERTICES * sizeof(TrajectoryLineVertex)};
    m_trajectory_frame_vertex_buffer = SDL_CreateGPUBuffer(m_gpu, &frame_vb_info);
    SDL_GPUTransferBufferCreateInfo frame_tb_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD,
                                                     .size = MAX_TRAJECTORY_FRAME_VERTICES * sizeof(TrajectoryLineVertex)};
    m_trajectory_frame_transfer_buffer = SDL_CreateGPUTransferBuffer(m_gpu, &frame_tb_info);

    // Two vertex uniform buffers: the camera (slot 0) and TrajectoryShaderConstants (slot 1)
    SDL_GPUShader* vert_shader = LoadShader(m_gpu, "Shaders/TrajectoryLine.vert", 0, 2);
    SDL_GPUShader* frag_shader = LoadShader(m_gpu, "Shaders/TrajectoryLine.frag", 0, 0);

    // The four inputs of TrajectoryLine.vert, in order
    SDL_GPUVertexAttribute vertex_attributes[4];
    vertex_attributes[0] = {.location = 0, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2,
                            .offset = offsetof(TrajectoryLineVertex, x)};
    vertex_attributes[1] = {.location = 1, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2,
                            .offset = offsetof(TrajectoryLineVertex, prevX)};
    vertex_attributes[2] = {.location = 2, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT4,
                            .offset = offsetof(TrajectoryLineVertex, r)};
    // drawn_at and fades_in sit next to each other, so the shader reads them as one float2
    vertex_attributes[3] = {.location = 3, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2,
                            .offset = offsetof(TrajectoryLineVertex, drawn_at)};

    SDL_GPUColorTargetDescription color_target = {};
    color_target.format = SDL_GetGPUSwapchainTextureFormat(m_gpu, m_window);
    color_target.blend_state.enable_blend = true;
    color_target.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
    color_target.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
    color_target.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
    color_target.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
    color_target.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
    color_target.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;

    SDL_GPUVertexBufferDescription vbo_desc = {
        .slot = 0, .pitch = sizeof(TrajectoryLineVertex), .input_rate = SDL_GPU_VERTEXINPUTRATE_VERTEX};

    SDL_GPUGraphicsPipelineCreateInfo pipeline_info = {};
    pipeline_info.target_info.num_color_targets = 1;
    pipeline_info.target_info.color_target_descriptions = &color_target;
    pipeline_info.vertex_shader = vert_shader;
    pipeline_info.fragment_shader = frag_shader;
    pipeline_info.primitive_type = SDL_GPU_PRIMITIVETYPE_LINELIST; // every 2 vertices = one separate segment
    pipeline_info.vertex_input_state.vertex_attributes = vertex_attributes;
    pipeline_info.vertex_input_state.num_vertex_attributes = 4;
    pipeline_info.vertex_input_state.vertex_buffer_descriptions = &vbo_desc;
    pipeline_info.vertex_input_state.num_vertex_buffers = 1;

    m_trajectory_line_pipeline = SDL_CreateGPUGraphicsPipeline(m_gpu, &pipeline_info);
    if (!m_trajectory_line_pipeline)
        std::cerr << "Trajectory line pipeline creation failed: " << SDL_GetError() << std::endl;
    SDL_ReleaseGPUShader(m_gpu, vert_shader);
    SDL_ReleaseGPUShader(m_gpu, frag_shader);
}

void RenderSystem::buildTrajectoryLines(const GameState& game_state, const UIState& ui_state)
{
    m_trajectory_new_vertices.clear();
    m_trajectory_frame_vertices.clear();
    const Trajectories& trajectories = game_state.getTrajectories();
    const UniverseState& universe_state = game_state.getUniverseState();
    const float zoom = game_state.getCameraState().zoom;
    const bool time_marks = ui_state.getShowTrajectoryTimeMarks();

    if (trajectories.status == TraceStatus::Idle)
    {
        m_trajectory_uploaded_vertex_count = 0; // nothing traced: draw nothing from the line buffer
        m_trajectory_points_uploaded.clear();
    }
    else
    {
        // A new trace: everything in the buffer belongs to the old one, so start it over
        if (trajectories.trace_id != m_trajectory_trace_id_uploaded || m_trajectory_points_uploaded.empty())
        {
            m_trajectory_trace_id_uploaded = trajectories.trace_id;
            m_trajectory_uploaded_vertex_count = 0;
            m_trajectory_points_uploaded.assign(trajectories.macro_body_paths.size() + 1, 0);
            m_trajectory_restart_upload = true;
        }
        const bool is_tracing = trajectories.status == TraceStatus::Tracing;

        // While there are lines, the real universe is exactly the one the trace started from (any edit wipes them),
        // so real macro body i is where path i starts
        for (size_t i = 0; i < trajectories.macro_body_paths.size() && i < universe_state.macro_bodies.size(); ++i)
        {
            const TracedPath& path = trajectories.macro_body_paths[i];
            const GravitationalBody& body = universe_state.macro_bodies[i];
            // A static body ignores every force: its line is a straight drift, so it gets its own colour
            TrajectoryStyle style = {body.isForceStatic ? TRAJECTORY_STATIC_COLOR : TRAJECTORY_MACRO_COLOR, true};
            appendTrajectorySegments(path, m_trajectory_points_uploaded[i], body.previousPosition, style,
                                     m_trajectory_new_vertices);
            m_trajectory_points_uploaded[i] = path.size();
            if (time_marks)
                appendTrajectoryTimeMarks(path, style, zoom, m_trajectory_frame_vertices);
            if (is_tracing)
                appendTrajectoryTip(path, m_trajectory_frame_vertices);
        }

        const TracedPath& ship_path = trajectories.ship_path;
        const Starship& ship = universe_state.ship;
        TrajectoryStyle ship_style = {TRAJECTORY_SHIP_COLOR, true};
        size_t& ship_points_uploaded = m_trajectory_points_uploaded.back();
        appendTrajectorySegments(ship_path, ship_points_uploaded, ship.center() - ship.movementThisTick(), ship_style,
                                 m_trajectory_new_vertices);
        ship_points_uploaded = ship_path.size();
        if (time_marks)
            appendTrajectoryTimeMarks(ship_path, ship_style, zoom, m_trajectory_frame_vertices);
        if (is_tracing)
            appendTrajectoryTip(ship_path, m_trajectory_frame_vertices);
    }

    // The planet being placed: live (no fade, no moving front), so it is simply rebuilt every frame; amber if it will
    // be static, like the traced lines
    const TracedPath& preview_path = trajectories.spawn_preview_path;
    if (!preview_path.empty())
    {
        bool will_be_static = ui_state.spawnSettings().is_force_static;
        TrajectoryStyle preview_style = {will_be_static ? TRAJECTORY_STATIC_COLOR : TRAJECTORY_PREVIEW_COLOR, false};
        appendTrajectorySegments(preview_path, 0, preview_path[0].position, preview_style,
                                 m_trajectory_frame_vertices);
        if (time_marks)
            appendTrajectoryTimeMarks(preview_path, preview_style, zoom, m_trajectory_frame_vertices);
    }
}

void RenderSystem::appendTrajectorySegments(const TracedPath& path, size_t first_point,
                                            const DynamoEngine::Vector2D& start_previous, const TrajectoryStyle& style,
                                            std::vector<TrajectoryLineVertex>& vertices_out)
{
    // The segment ending at point k runs from point k - 1 (point 0 has none). Only the very first point moves with
    // the drawn body (the shader interpolates previous -> current); every other point is fixed in the world.
    for (size_t k = std::max<size_t>(first_point, 1); k < path.size(); ++k)
    {
        const TracedPoint& from = path[k - 1];
        const TracedPoint& to = path[k];
        DynamoEngine::Vector2D from_previous = (k == 1) ? start_previous : from.position;
        vertices_out.push_back(trajectoryVertex(from.position, from_previous, style.color, TRAJECTORY_LINE_ALPHA,
                                                from.drawn_at, style.fades_in));
        vertices_out.push_back(trajectoryVertex(to.position, to.position, style.color, TRAJECTORY_LINE_ALPHA,
                                                to.drawn_at, style.fades_in));
    }
}

void RenderSystem::appendTrajectoryTimeMarks(const TracedPath& path, const TrajectoryStyle& style, float zoom,
                                             std::vector<TrajectoryLineVertex>& vertices_out)
{
    // Points are evenly spaced in time, so the marks spread out where the body is fast and bunch up where it's slow
    const size_t points_per_mark = static_cast<size_t>(
        std::lround(TRAJECTORY_TIME_MARK_SECONDS * INV_PHYSICS_TIME_STEP / TRAJECTORY_TICKS_PER_POINT));
    const double half_length = TRAJECTORY_TIME_MARK_HALF_LENGTH_PX / zoom; // world units, constant on screen
    for (size_t k = points_per_mark; k < path.size(); k += points_per_mark)
    {
        // The path's direction here (from the point before to the point after), turned a quarter turn
        const DynamoEngine::Vector2D& before = path[k - 1].position;
        const DynamoEngine::Vector2D& after = path[std::min(k + 1, path.size() - 1)].position;
        DynamoEngine::Vector2D direction = after - before;
        double length = direction.magnitude();
        if (length <= 0.0)
        {
            continue; // a body at rest has no direction to cross
        }
        DynamoEngine::Vector2D across =
            DynamoEngine::Vector2D(-direction.y_val, direction.x_val) * (half_length / length);
        DynamoEngine::Vector2D a = path[k].position + across;
        DynamoEngine::Vector2D b = path[k].position - across;
        vertices_out.push_back(trajectoryVertex(a, a, style.color, TRAJECTORY_LINE_ALPHA, path[k].drawn_at,
                                                style.fades_in));
        vertices_out.push_back(trajectoryVertex(b, b, style.color, TRAJECTORY_LINE_ALPHA, path[k].drawn_at,
                                                style.fades_in));
    }
}

void RenderSystem::appendTrajectoryTip(const TracedPath& path, std::vector<TrajectoryLineVertex>& vertices_out)
{
    // The newest segments again, on top of the line, from transparent to fully bright at the very tip
    const size_t tip_segments = static_cast<size_t>(TRAJECTORY_TIP_SEGMENTS);
    size_t first_tip_point = (path.size() > tip_segments + 1) ? path.size() - tip_segments : 1;
    for (size_t k = first_tip_point; k < path.size(); ++k)
    {
        const DynamoEngine::Vector2D& from = path[k - 1].position;
        const DynamoEngine::Vector2D& to = path[k].position;
        float from_alpha = (float)(k - first_tip_point) / TRAJECTORY_TIP_SEGMENTS;
        float to_alpha = std::min((float)(k - first_tip_point + 1) / TRAJECTORY_TIP_SEGMENTS, 1.0f);
        vertices_out.push_back(trajectoryVertex(from, from, TRAJECTORY_TIP_COLOR, from_alpha, 0.0f, false));
        vertices_out.push_back(trajectoryVertex(to, to, TRAJECTORY_TIP_COLOR, to_alpha, 0.0f, false));
    }
}

void RenderSystem::uploadTrajectoryLines(SDL_GPUCommandBuffer* cmdbuf)
{
    // Never upload more than the buffers hold (PhysicsConstants.hpp checks the worst case fits, so this is a bug guard)
    uint32_t room = MAX_TRAJECTORY_LINE_VERTICES - m_trajectory_uploaded_vertex_count;
    if (m_trajectory_new_vertices.size() > room)
    {
        printf("Too many trajectory line vertices: %zu new, room for %u\n", m_trajectory_new_vertices.size(), room);
        m_trajectory_new_vertices.resize(room);
    }
    if (m_trajectory_frame_vertices.size() > MAX_TRAJECTORY_FRAME_VERTICES)
    {
        printf("Too many per-frame trajectory vertices %zu > %u\n", m_trajectory_frame_vertices.size(),
               MAX_TRAJECTORY_FRAME_VERTICES);
        m_trajectory_frame_vertices.resize(MAX_TRAJECTORY_FRAME_VERTICES);
    }
    if (m_trajectory_new_vertices.empty() && m_trajectory_frame_vertices.empty())
    {
        return;
    }

    SDL_GPUCopyPass* copy_pass = SDL_BeginGPUCopyPass(cmdbuf);

    // Traced lines: only this frame's new segments, written right after the ones already in the buffer
    if (!m_trajectory_new_vertices.empty())
    {
        const uint32_t bytes = (uint32_t)(m_trajectory_new_vertices.size() * sizeof(TrajectoryLineVertex));
        void* map = SDL_MapGPUTransferBuffer(m_gpu, m_trajectory_line_transfer_buffer, true);
        SDL_memcpy(map, m_trajectory_new_vertices.data(), bytes);
        SDL_UnmapGPUTransferBuffer(m_gpu, m_trajectory_line_transfer_buffer);

        SDL_GPUTransferBufferLocation src = {.transfer_buffer = m_trajectory_line_transfer_buffer, .offset = 0};
        SDL_GPUBufferRegion dst = {.buffer = m_trajectory_line_vertex_buffer,
                                   .offset = (uint32_t)(m_trajectory_uploaded_vertex_count * sizeof(TrajectoryLineVertex)),
                                   .size = bytes};
        // cycle: when starting over, earlier frames may still be drawing the OLD trace from this buffer, so take a
        // fresh one (true). When appending, the old segments must stay where they are (false); that is safe because
        // we only write past the part earlier frames draw.
        SDL_UploadToGPUBuffer(copy_pass, &src, &dst, m_trajectory_restart_upload);
        m_trajectory_uploaded_vertex_count += (uint32_t)m_trajectory_new_vertices.size();
        m_trajectory_restart_upload = false;
    }

    // Per-frame lines: rebuilt every frame, so uploaded whole (cycle = true: last frame may still be reading them)
    if (!m_trajectory_frame_vertices.empty())
    {
        const uint32_t bytes = (uint32_t)(m_trajectory_frame_vertices.size() * sizeof(TrajectoryLineVertex));
        void* map = SDL_MapGPUTransferBuffer(m_gpu, m_trajectory_frame_transfer_buffer, true);
        SDL_memcpy(map, m_trajectory_frame_vertices.data(), bytes);
        SDL_UnmapGPUTransferBuffer(m_gpu, m_trajectory_frame_transfer_buffer);

        SDL_GPUTransferBufferLocation src = {.transfer_buffer = m_trajectory_frame_transfer_buffer, .offset = 0};
        SDL_GPUBufferRegion dst = {.buffer = m_trajectory_frame_vertex_buffer, .offset = 0, .size = bytes};
        SDL_UploadToGPUBuffer(copy_pass, &src, &dst, true);
    }

    SDL_EndGPUCopyPass(copy_pass);
}

void RenderSystem::renderTrajectoryLines(SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf,
                                         const CameraState& camera_state, float clock)
{
    if (m_trajectory_uploaded_vertex_count == 0 && m_trajectory_frame_vertices.empty())
        return;

    SDL_BindGPUGraphicsPipeline(pass, m_trajectory_line_pipeline);
    CameraConstants camera_constants = buildCameraConstants(camera_state, camera_state.offset);
    SDL_PushGPUVertexUniformData(cmdbuf, 0, &camera_constants, sizeof(camera_constants));
    TrajectoryShaderConstants trajectory_constants = {clock, TRAJECTORY_FADE_IN_SECONDS, TRAJECTORY_FADE_START_ALPHA,
                                                      0.0f};
    SDL_PushGPUVertexUniformData(cmdbuf, 1, &trajectory_constants, sizeof(trajectory_constants));

    if (m_trajectory_uploaded_vertex_count > 0)
    {
        SDL_GPUBufferBinding vbo = {.buffer = m_trajectory_line_vertex_buffer, .offset = 0};
        SDL_BindGPUVertexBuffers(pass, 0, &vbo, 1);
        SDL_DrawGPUPrimitives(pass, m_trajectory_uploaded_vertex_count, 1, 0, 0);
    }
    if (!m_trajectory_frame_vertices.empty())
    {
        SDL_GPUBufferBinding vbo = {.buffer = m_trajectory_frame_vertex_buffer, .offset = 0};
        SDL_BindGPUVertexBuffers(pass, 0, &vbo, 1);
        SDL_DrawGPUPrimitives(pass, (uint32_t)m_trajectory_frame_vertices.size(), 1, 0, 0);
    }
}
