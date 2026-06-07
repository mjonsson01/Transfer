// File: Transfer/src/Systems/RenderSystem.cpp

#include "Systems/RenderSystem.hpp"

// Constructor: Initializes SDL Window and GPU
RenderSystem::RenderSystem()
{
    SDL_InitSubSystem(SDL_INIT_VIDEO);
    TTF_Init();

    // 1. Window & Device Setup
    int window_flags = SDL_WINDOW_RESIZABLE | SDL_WINDOW_HIGH_PIXEL_DENSITY;
    window = SDL_CreateWindow("Transfer", SCREEN_WIDTH, SCREEN_HEIGHT, window_flags);

    SDL_GPUShaderFormat formats = SDL_GPU_SHADERFORMAT_SPIRV;
    gpu = SDL_CreateGPUDevice(formats, true, nullptr);
    if (gpu)
        SDL_ClaimWindowForGPUDevice(gpu, window);

    // 2. Resource/Font Setup
    UIFontRegular = TTF_OpenFont(Utilities::GetResourcePath("Fonts/SpaceMono-Regular.ttf").c_str(), 18);
    UIFontTitle = TTF_OpenFont(Utilities::GetResourcePath("Fonts/SpaceMono-Bold.ttf").c_str(), 32);

    createGravBodyGPUBuffer();
    createTwinklingStarGPUBuffer();
    createTwinklingStarField();
    if (gpu)
    {
        SDL_GPUCommandBuffer* initCmdBuf = SDL_AcquireGPUCommandBuffer(gpu);
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
void RenderSystem::CleanUp()
{
    if (window)
    {
        SDL_DestroyWindow(window);
        window = nullptr;
    }
    // 1. Release GPU-specific resources
    if (unifiedBodyVertexBuffer != nullptr)
        SDL_ReleaseGPUBuffer(gpu, unifiedBodyVertexBuffer);
    if (twinklingStarVertexBuffer != nullptr)
        SDL_ReleaseGPUBuffer(gpu, twinklingStarVertexBuffer);
    if (cameraUniformBuffer != nullptr)
        SDL_ReleaseGPUBuffer(gpu, cameraUniformBuffer);
    if (unifiedBodyTransferBuffer != nullptr)
        SDL_ReleaseGPUTransferBuffer(gpu, unifiedBodyTransferBuffer);
    if (twinklingStarTransferBuffer != nullptr)
        SDL_ReleaseGPUTransferBuffer(gpu, twinklingStarTransferBuffer);
    if (unifiedBodyPipeline != nullptr)
        SDL_ReleaseGPUGraphicsPipeline(gpu, unifiedBodyPipeline);
    if (twinklingStarPipeline != nullptr)
        SDL_ReleaseGPUGraphicsPipeline(gpu, twinklingStarPipeline);

    // 2. Finally, destroy the GPU device
    if (gpu != nullptr)
        SDL_DestroyGPUDevice(gpu);

    if (window)
    {
        SDL_DestroyWindow(window);
        window = nullptr;
    }
}

SDL_GPUShader* RenderSystem::LoadShader(SDL_GPUDevice* device, const char* fileName)
{
    size_t size;
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
                                    .entrypoint = "main",
                                    .format = SDL_GPU_SHADERFORMAT_SPIRV, // Windows default
                                    .stage = (std::string(fileName).find("vert") != std::string::npos)
                                                 ? SDL_GPU_SHADERSTAGE_VERTEX
                                                 : SDL_GPU_SHADERSTAGE_FRAGMENT};

    SDL_GPUShader* shader = SDL_CreateGPUShader(device, &info);
    SDL_free(code);
    return shader;
}
// --------- RENDER FULL FRAME METHOD --------- //

void RenderSystem::RenderFullFrame(GameState& gameState, UIState& UIState,
                                   const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope)
{

    SDL_GPUCommandBuffer* cmdbuf = SDL_AcquireGPUCommandBuffer(gpu);

    SceneIdentifier current_scene = UIState.getCurrentSceneID();

    if (current_scene == SceneIdentifier::GAME_SCENE)
    {
        uploadBodies(gameState, UIState, cmdbuf);
    }
    // Acquire the display target
    SDL_GPUTexture* swapchainTexture = nullptr;
    Uint32 w = 0, h = 0;
    if (!SDL_AcquireGPUSwapchainTexture(cmdbuf, window, &swapchainTexture, &w, &h))
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

        if (current_scene == SceneIdentifier::GAME_SCENE)
        {
            // Update your renderGameFrame signature to match
            renderGameFrame(gameState, UIState, allUIElementsInScope, pass, cmdbuf);
        }
        else
        {
            renderNonGameFrame(gameState, UIState, allUIElementsInScope, pass, cmdbuf);
        }

        SDL_EndGPURenderPass(pass);
    }

    SDL_SubmitGPUCommandBuffer(cmdbuf);
}

void RenderSystem::renderGameFrame(GameState& gameState, UIState& UIState,
                                   const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope,
                                   SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf)
{
    renderTwinklingStarField(pass);
    renderBodies(gameState, UIState, pass, cmdbuf);
}
void RenderSystem::renderNonGameFrame(GameState& gameState, UIState& UIState,
                                      const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope,
                                      SDL_GPURenderPass* pass, SDL_GPUCommandBuffer* cmdbuf)
{
}

void RenderSystem::uploadBodies(GameState& gameState, UIState& UIState, SDL_GPUCommandBuffer* cmdbuf)
{
    auto& particles = gameState.getParticles();
    auto& bodies = gameState.getMacroBodies();

    std::vector<UnifiedBodyVertex> vertex_data;
    vertex_data.reserve(particles.size() + bodies.size());

    for (auto& p : particles)
    {
        if (p.visible)
        {
            vertex_data.push_back(p.toUnifiedVertex());
        }
    }

    for (auto& b : bodies)
    {
        if (b.visible)
        {
            vertex_data.push_back(b.toUnifiedVertex());
        }
    }

    appendPreviewBodies(vertex_data, UIState);

    if (vertex_data.size() > MAX_UNIFIED_BODIES)
    {
        printf("Too many bodies! %zu > %d\n", vertex_data.size(), MAX_UNIFIED_BODIES);
        return;
    }
    // Copy pass
    if (!vertex_data.empty())
    {
        void* map = SDL_MapGPUTransferBuffer(gpu, unifiedBodyTransferBuffer, false);
        SDL_memcpy(map, vertex_data.data(), vertex_data.size() * sizeof(UnifiedBodyVertex));
        SDL_UnmapGPUTransferBuffer(gpu, unifiedBodyTransferBuffer);

        SDL_GPUCopyPass* copyPass = SDL_BeginGPUCopyPass(cmdbuf);

        SDL_GPUTransferBufferLocation src = {.transfer_buffer = unifiedBodyTransferBuffer, .offset = 0};
        SDL_GPUBufferRegion dst = {.buffer = unifiedBodyVertexBuffer,
                                   .offset = 0,
                                   .size = (uint32_t)(vertex_data.size() * sizeof(UnifiedBodyVertex))};

        SDL_UploadToGPUBuffer(copyPass, &src, &dst, false);
        SDL_EndGPUCopyPass(copyPass);
    }
}

void RenderSystem::renderBodies(GameState& gameState, UIState& UIState, SDL_GPURenderPass* pass,
                                SDL_GPUCommandBuffer* cmdbuf)
{
    // Quickly count how many total instances are active for drawing
    uint32_t instance_count = 0;
    for (auto& p : gameState.getParticles())
        if (p.visible)
            instance_count++;
    for (auto& b : gameState.getMacroBodies())
        if (b.visible)
            instance_count++;

    if (UIState.getMutableInputState().isPreviewingMacro)
    {
        instance_count++;
    }

    // Safety check and raw drawing
    if (instance_count > 0 && instance_count <= MAX_UNIFIED_BODIES)
    {
        SDL_BindGPUGraphicsPipeline(pass, unifiedBodyPipeline);

        SDL_GPUBufferBinding vbo = {.buffer = unifiedBodyVertexBuffer, .offset = 0};
        SDL_BindGPUVertexBuffers(pass, 0, &vbo, 1);

        SDL_DrawGPUPrimitives(pass,
                              6,              // vertices per quad
                              instance_count, // instances
                              0, 0);
    }
}

void RenderSystem::appendPreviewBodies(std::vector<UnifiedBodyVertex>& vertexData, UIState& UIState)
{
    InputState& input_state = UIState.getMutableInputState();
    if (input_state.isPreviewingMacro)
    {
        // create pseudo body and convert to unified body vertex to pass to rendering pipeline
        GravitationalBody new_preview_grav_body = {};
        new_preview_grav_body.mass = input_state.selectedMass;
        new_preview_grav_body.radius = input_state.selectedRadius;
        if (input_state.isPreviewingWithInitialVelocity)
        {
            new_preview_grav_body.position = input_state.mouseDragStartPosition;
        }
        else
        {
            new_preview_grav_body.position = input_state.mouseCurrPosition;
        }
        new_preview_grav_body.isPreview = true;

        UnifiedBodyVertex new_preview_unified_body_vertex = new_preview_grav_body.toUnifiedVertex();
        vertexData.push_back(new_preview_unified_body_vertex);
    }
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

void RenderSystem::createGravBodyGPUBuffer()
{

    SDL_GPUBufferCreateInfo vertex_buffer_info = {.usage = SDL_GPU_BUFFERUSAGE_VERTEX,
                                                  .size = MAX_UNIFIED_BODIES * sizeof(UnifiedBodyVertex)};
    unifiedBodyVertexBuffer = SDL_CreateGPUBuffer(gpu, &vertex_buffer_info);

    SDL_GPUBufferCreateInfo camera_buffer_info = {.usage = SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ,
                                                  .size = sizeof(CameraConstants)};
    cameraUniformBuffer = SDL_CreateGPUBuffer(gpu, &camera_buffer_info);

    SDL_GPUTransferBufferCreateInfo tb_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD,
                                               .size = MAX_UNIFIED_BODIES * sizeof(UnifiedBodyVertex)};
    unifiedBodyTransferBuffer = SDL_CreateGPUTransferBuffer(gpu, &tb_info);

    SDL_GPUShader* vert_shader = LoadShader(gpu, "Shaders/UnifiedGravBody.vert.spv");
    SDL_GPUShader* frag_shader = LoadShader(gpu, "Shaders/UnifiedGravBody.frag.spv");

    // 6. Define Pipeline State
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
    color_target.format = SDL_GetGPUSwapchainTextureFormat(gpu, window);
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

    // Change to TRIANGLESTRIP for your Quads!
    pipeline_info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;

    pipeline_info.vertex_input_state.vertex_attributes = vertex_attributes;
    pipeline_info.vertex_input_state.num_vertex_attributes = 8;

    SDL_GPUVertexBufferDescription vbo_desc = {
        .slot = 0, .pitch = sizeof(UnifiedBodyVertex), .input_rate = SDL_GPU_VERTEXINPUTRATE_INSTANCE};
    pipeline_info.vertex_input_state.vertex_buffer_descriptions = &vbo_desc;
    pipeline_info.vertex_input_state.num_vertex_buffers = 1;

    unifiedBodyPipeline = SDL_CreateGPUGraphicsPipeline(gpu, &pipeline_info);

    // 7. Cleanup temp shader handles
    SDL_ReleaseGPUShader(gpu, vert_shader);
    SDL_ReleaseGPUShader(gpu, frag_shader);
}

void RenderSystem::createTwinklingStarGPUBuffer()
{

    SDL_GPUBufferCreateInfo vertex_buffer_info = {.usage = SDL_GPU_BUFFERUSAGE_VERTEX,
                                                  .size = STAR_NUM * sizeof(TwinklingStarVertex)};
    twinklingStarVertexBuffer = SDL_CreateGPUBuffer(gpu, &vertex_buffer_info);

    SDL_GPUTransferBufferCreateInfo tb_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD,
                                               .size = STAR_NUM * sizeof(TwinklingStarVertex)};

    twinklingStarTransferBuffer = SDL_CreateGPUTransferBuffer(gpu, &tb_info);

    SDL_GPUShader* vert_shader = LoadShader(gpu, "Shaders/TwinklingStar.vert.spv");
    SDL_GPUShader* frag_shader = LoadShader(gpu, "Shaders/TwinklingStar.frag.spv");

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
    color_target.format = SDL_GetGPUSwapchainTextureFormat(gpu, window);
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

    // Change to TRIANGLESTRIP for your Quads!
    pipeline_info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;

    pipeline_info.vertex_input_state.vertex_attributes = vertex_attributes;
    pipeline_info.vertex_input_state.num_vertex_attributes = 5;

    SDL_GPUVertexBufferDescription vbo_desc = {
        .slot = 0, .pitch = sizeof(TwinklingStarVertex), .input_rate = SDL_GPU_VERTEXINPUTRATE_INSTANCE};
    pipeline_info.vertex_input_state.vertex_buffer_descriptions = &vbo_desc;
    pipeline_info.vertex_input_state.num_vertex_buffers = 1;

    twinklingStarPipeline = SDL_CreateGPUGraphicsPipeline(gpu, &pipeline_info);

    SDL_ReleaseGPUShader(gpu, vert_shader);
    SDL_ReleaseGPUShader(gpu, frag_shader);
}

void RenderSystem::createTwinklingStarField()
{
    twinklingStars.clear();
    twinklingStars.reserve(STAR_NUM);

    std::mt19937 rng(12345);

    std::uniform_real_distribution<float> posX(0, SCREEN_WIDTH);
    std::uniform_real_distribution<float> posY(0, SCREEN_HEIGHT);

    std::uniform_real_distribution<float> radiusDist(1.0f, 4.0f);

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

        twinklingStars.push_back(star);
    }
}

void RenderSystem::uploadTwinklingStarField(SDL_GPUCommandBuffer* cmdbuf)
{

    void* map = SDL_MapGPUTransferBuffer(gpu, twinklingStarTransferBuffer, false);

    SDL_memcpy(map, twinklingStars.data(), twinklingStars.size() * sizeof(TwinklingStarVertex));

    SDL_UnmapGPUTransferBuffer(gpu, twinklingStarTransferBuffer);

    SDL_GPUCopyPass* copyPass = SDL_BeginGPUCopyPass(cmdbuf);

    SDL_GPUTransferBufferLocation src = {.transfer_buffer = twinklingStarTransferBuffer, .offset = 0};

    SDL_GPUBufferRegion dst = {.buffer = twinklingStarVertexBuffer,
                               .offset = 0,
                               .size = (uint32_t)(twinklingStars.size() * sizeof(TwinklingStarVertex))};

    SDL_UploadToGPUBuffer(copyPass, &src, &dst, false);

    SDL_EndGPUCopyPass(copyPass);
}

void RenderSystem::renderTwinklingStarField(SDL_GPURenderPass* pass)
{
    SDL_BindGPUGraphicsPipeline(pass, twinklingStarPipeline);

    SDL_GPUBufferBinding vbo = {.buffer = twinklingStarVertexBuffer, .offset = 0};
    SDL_BindGPUVertexBuffers(pass, 0, &vbo, 1);
    SDL_DrawGPUPrimitives(pass, 6, STAR_NUM, 0, 0);
}