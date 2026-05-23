// File: Transfer/src/Systems/RenderSystem.cpp

#include "Systems/RenderSystem.hpp"

// Constructor: Initializes SDL Window and Renderer
RenderSystem::RenderSystem()
{
    SDL_InitSubSystem(SDL_INIT_VIDEO);
    TTF_Init();

    // 1. Window & Device Setup
    int window_flags = SDL_WINDOW_RESIZABLE | SDL_WINDOW_HIGH_PIXEL_DENSITY;
    window = SDL_CreateWindow("Transfer", SCREEN_WIDTH, SCREEN_HEIGHT, window_flags);

    SDL_GPUShaderFormat formats = SDL_GPU_SHADERFORMAT_SPIRV; // Windows standard
    gpu = SDL_CreateGPUDevice(formats, true, nullptr);
    if (gpu)
        SDL_ClaimWindowForGPUDevice(gpu, window);

    renderer = SDL_CreateRenderer(window, nullptr);

    // 2. Resource/Font Setup
    UIFontRegular = TTF_OpenFont(Utilities::GetResourcePath("Fonts/SpaceMono-Regular.ttf").c_str(), 18);
    UIFontTitle = TTF_OpenFont(Utilities::GetResourcePath("Fonts/SpaceMono-Bold.ttf").c_str(), 32);

    // 3. Create GPU Buffers (VRAM)typedef struct SDL_GPUBufferCreateInfo
    SDL_GPUBufferCreateInfo;

    SDL_GPUBufferCreateInfo vertex_buffer_info = {.usage = SDL_GPU_BUFFERUSAGE_VERTEX,
                                                  .size = MAX_UNIFIED_BODIES * sizeof(UnifiedBodyVertex)};
    unified_body_vertex_buffer = SDL_CreateGPUBuffer(gpu, &vertex_buffer_info);

    SDL_GPUBufferCreateInfo camera_buffer_info = {.usage = SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ,
                                                  .size = sizeof(CameraConstants)};
    camera_uniform_buffer = SDL_CreateGPUBuffer(gpu, &camera_buffer_info);

    // 4. Create Transfer Buffer (The Loading Dock)
    SDL_GPUTransferBufferCreateInfo tb_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD,
                                               .size = MAX_UNIFIED_BODIES * sizeof(UnifiedBodyVertex)};
    body_transfer_buffer = SDL_CreateGPUTransferBuffer(gpu, &tb_info);

    // 5. Shader Loading (MUST happen before Pipeline)
    SDL_GPUShader* vertShader = LoadShader(gpu, "Shaders/Body.vert.spv");
    SDL_GPUShader* fragShader = LoadShader(gpu, "Shaders/Body.frag.spv");

    // 6. Define Pipeline State
    SDL_GPUVertexAttribute vertex_attributes[8];

    // for (int i = 0; i < 8; i++)
    // {
    //     vertex_attributes[i].instance_step_rate = 1; // IMPORTANT: 1 = Once per body
    // }
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
    vertex_attributes[3].offset = offsetof(UnifiedBodyVertex, log_mass);
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
    pipeline_info.vertex_shader = vertShader;
    pipeline_info.fragment_shader = fragShader;

    // Change to TRIANGLESTRIP for your Quads!
    pipeline_info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;

    pipeline_info.vertex_input_state.vertex_attributes = vertex_attributes;
    pipeline_info.vertex_input_state.num_vertex_attributes = 8;

    SDL_GPUVertexBufferDescription vbo_desc = {
        .slot = 0, .pitch = sizeof(UnifiedBodyVertex), .input_rate = SDL_GPU_VERTEXINPUTRATE_INSTANCE};
    pipeline_info.vertex_input_state.vertex_buffer_descriptions = &vbo_desc;
    pipeline_info.vertex_input_state.num_vertex_buffers = 1;

    unified_body_pipeline = SDL_CreateGPUGraphicsPipeline(gpu, &pipeline_info);

    // 7. Cleanup temp shader handles
    SDL_ReleaseGPUShader(gpu, vertShader);
    SDL_ReleaseGPUShader(gpu, fragShader);
}
// Destructor: Cleans up SDL Window and Renderer
RenderSystem::~RenderSystem()
{
    // TTF_Quit() and SDL_QUIT() handled at the Game level
}

// --------- CLEANUP METHOD --------- //
void RenderSystem::CleanUp()
{
    clearCachedCircleTextures();
    if (renderer)
    {
        SDL_DestroyRenderer(renderer);
        renderer = nullptr;
    }
    if (window)
    {
        SDL_DestroyWindow(window);
        window = nullptr;
    }
    // 1. Release GPU-specific resources
    if (unified_body_vertex_buffer != nullptr)
        SDL_ReleaseGPUBuffer(gpu, unified_body_vertex_buffer);
    if (camera_uniform_buffer != nullptr)
        SDL_ReleaseGPUBuffer(gpu, camera_uniform_buffer);
    if (body_transfer_buffer != nullptr)
        SDL_ReleaseGPUTransferBuffer(gpu, body_transfer_buffer);
    if (unified_body_pipeline != nullptr)
        SDL_ReleaseGPUGraphicsPipeline(gpu, unified_body_pipeline);

    // 2. Finally, destroy the GPU device
    if (gpu != nullptr)
        SDL_DestroyGPUDevice(gpu);

    // 3. Clean up legacy stuff
    if (renderer)
    {
        SDL_DestroyRenderer(renderer);
        renderer = nullptr;
    }
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
    SceneIdentifier current_scene = UIState.getCurrentSceneID();
    if (current_scene == SceneIdentifier::GAME_SCENE)
    {
        renderGameFrame(gameState, UIState, allUIElementsInScope);
    }
    else
    {
        renderNonGameFrame(gameState, UIState, allUIElementsInScope);
    }

    SDL_RenderPresent(renderer);
}
void RenderSystem::renderGameFrame(GameState& gameState, UIState& UIState,
                                   const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope)
{
    // SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // black background
    // SDL_RenderClear(renderer);

    // render starry background
    // updateStars();
    // renderStars();
    // Render all preview bodies
    // renderPreviewBodies(UIState);
    // Render all the bodies (particles)
    renderBodies(gameState);

    // renderUIElements(UIState, allUIElementsInScope);

    // TEST
    // SDL_GPUCommandBuffer* cmdbuf = SDL_AcquireGPUCommandBuffer(gpu);
    // if (!cmdbuf)
    //     return;

    // SDL_GPUTexture* swapchainTexture = nullptr;
    // Uint32 w, h;

    // if (SDL_AcquireGPUSwapchainTexture(cmdbuf, window, &swapchainTexture, &w, &h))
    // {
    //     if (swapchainTexture != nullptr)
    //     {
    //         SDL_GPUColorTargetInfo color_info = {};
    //         color_info.texture = swapchainTexture;
    //         // Set this to something distinct like 1.0, 0.0, 0.0, 1.0 (Bright Red)
    //         color_info.clear_color = {1.0f, 0.0f, 0.0f, 1.0f};
    //         color_info.load_op = SDL_GPU_LOADOP_CLEAR;
    //         color_info.store_op = SDL_GPU_STOREOP_STORE;

    //         SDL_GPURenderPass* render_pass = SDL_BeginGPURenderPass(cmdbuf, &color_info, 1, nullptr);

    //         // Just clear the screen, don't bind pipelines or draw yet
    //         SDL_EndGPURenderPass(render_pass);
    //     }
    // }

    // SDL_SubmitGPUCommandBuffer(cmdbuf);
}
void RenderSystem::renderNonGameFrame(GameState& gameState, UIState& UIState,
                                      const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope)
{
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // black background
    SDL_RenderClear(renderer);
    renderUIElements(UIState, allUIElementsInScope);
}

void RenderSystem::renderPreviewBodies(UIState& UIState)
{
    InputState& input_state = UIState.getMutableInputState();
    if (input_state.isPreviewingMacro)
    {
        // create pseudo body.
        if (input_state.selectedRadius <= 1.0)
        {
            // also probably throw an error toast here somehow
            return;
        }
        GravitationalBody body = {};
        body.mass = input_state.selectedMass;
        body.radius = input_state.selectedRadius;
        // rendering preview body
        if (input_state.isPreviewingWithInitialVelocity)
        {
            body.position = input_state.mouseDragStartPosition;
        }
        else
        {
            body.position = input_state.mouseCurrPosition;
        }
        SDL_Color color = getColorForProperty(body);
        SDL_Texture* tex = circleTextureCache[static_cast<int>(body.radius)];
        SDL_SetTextureColorMod(tex, color.r, color.g, color.b);
        SDL_SetTextureAlphaMod(tex, color.a);

        float render_x = body.position.xVal;
        float render_y = body.position.yVal;
        render_x = std::round(render_x);
        render_y = std::round(render_y);
        float r = static_cast<float>(body.radius);
        SDL_FRect dstRect = {render_x - r, render_y - r, r * 2, r * 2};
        SDL_RenderTexture(renderer, tex, nullptr, &dstRect);
        if (input_state.isPreviewingWithInitialVelocity)
        {
            renderDragLine(body.position, input_state.mouseCurrPosition);
        }
    }
}
// --------- RENDER GRAVITATIONAL BODIES METHOD --------- //
void RenderSystem::renderBodies(GameState& gameState)
{
    float alpha = gameState.getAlpha();
    auto& particles = gameState.getParticles();
    auto& bodies = gameState.getMacroBodies();
    std::vector<UnifiedBodyVertex> vertex_data;
    vertex_data.reserve(gameState.getParticles().size() + gameState.getMacroBodies().size());
    for (auto& p : particles)
    {
        if (p.visible)
            vertex_data.push_back(p.toUnifiedVertex());
    }
    for (auto& b : bodies)
    {
        if (b.visible)
            vertex_data.push_back(b.toUnifiedVertex());
    }
    if (vertex_data.empty())
        return;

    bool cycle = false;
    void* map = SDL_MapGPUTransferBuffer(gpu, body_transfer_buffer, cycle);
    SDL_memcpy(map, vertex_data.data(), vertex_data.size() * sizeof(UnifiedBodyVertex));
    SDL_UnmapGPUTransferBuffer(gpu, body_transfer_buffer);

    SDL_GPUCommandBuffer* cmdbuf = SDL_AcquireGPUCommandBuffer(gpu);
    SDL_GPUCopyPass* copyPass = SDL_BeginGPUCopyPass(cmdbuf);
    SDL_GPUTransferBufferLocation src = {.transfer_buffer = body_transfer_buffer, .offset = 0};
    SDL_GPUBufferRegion dst = {.buffer = unified_body_vertex_buffer,
                               .offset = 0,
                               .size = (uint32_t)(vertex_data.size() * sizeof(UnifiedBodyVertex))};
    SDL_UploadToGPUBuffer(copyPass, &src, &dst, false);
    SDL_EndGPUCopyPass(copyPass);

    // 2. Render Pass
    SDL_GPUTexture* swapchainTexture;
    Uint32 w, h;
    if (SDL_AcquireGPUSwapchainTexture(cmdbuf, window, &swapchainTexture, &w, &h))
    {
        SDL_GPUColorTargetInfo color_info = {.texture = swapchainTexture,
                                             .clear_color = {0.0f, 0.0f, 0.0f, 1.0f},
                                             .load_op = SDL_GPU_LOADOP_CLEAR,
                                             .store_op = SDL_GPU_STOREOP_STORE};
        SDL_GPURenderPass* pass = SDL_BeginGPURenderPass(cmdbuf, &color_info, 1, nullptr);

        SDL_BindGPUGraphicsPipeline(pass, unified_body_pipeline);

        SDL_GPUBufferBinding vbo = {.buffer = unified_body_vertex_buffer, .offset = 0};
        SDL_BindGPUVertexBuffers(pass, 0, &vbo, 1);

        // printf("Total bodies: %zu\n", vertex_data.size());
        SDL_DrawGPUPrimitives(pass, 6, (uint32_t)vertex_data.size(), 0, 0);

        SDL_EndGPURenderPass(pass);
    }
    SDL_SubmitGPUCommandBuffer(cmdbuf);
}

// --------- RENDER DRAGLINES METHOD --------- //

void RenderSystem::renderDragLine(Vector2D lineStart, Vector2D lineEnd)
{
    SDL_Color color = ColorLibrary::White;
    SDL_FColor f_color = {1.0f, 1.0f, 1.0f, 1.0f};

    float dx = static_cast<float>(lineEnd.xVal - lineStart.xVal);
    float dy = static_cast<float>(lineEnd.yVal - lineStart.yVal);
    float length = static_cast<float>((lineEnd - lineStart).magnitude());
    if (length <= EPSILON)
    {
        return;
    }
    dx /= length;
    dy /= length;
    float px = -dy;
    float py = dx;
    const float thickness = 6.0f;
    const float arrow_length = 20.0f;
    const float arrow_width = 24.0f;
    float half_T = thickness / 2.0f;
    // --- IMPORTANT: shorten line so arrow has space ---
    Vector2D line_end = {lineEnd.xVal - dx * arrow_length, lineEnd.yVal - dy * arrow_length};
    SDL_Vertex line_verts[4];

    line_verts[0].position = {float(lineStart.xVal) + px * half_T, float(lineStart.yVal) + py * half_T};
    line_verts[1].position = {float(lineStart.xVal) - px * half_T, float(lineStart.yVal) - py * half_T};
    line_verts[2].position = {float(line_end.xVal) - px * half_T, float(line_end.yVal) - py * half_T};
    line_verts[3].position = {float(line_end.xVal) + px * half_T, float(line_end.yVal) + py * half_T};

    for (int i = 0; i < 4; i++)
        line_verts[i].color = f_color;
    int line_indices[6] = {0, 1, 2, 0, 2, 3};
    SDL_RenderGeometry(renderer, nullptr, line_verts, 4, line_indices, 6);

    // --- Arrowhead (triangle) ---
    float half_arrow_W = arrow_width / 2.0f;

    SDL_Vertex arrow_verts[3];

    Vector2D tip = lineEnd;
    Vector2D base = line_end;

    arrow_verts[0].position = {float(tip.xVal), float(tip.yVal)};

    arrow_verts[1].position = {float(base.xVal) + px * (arrow_width / 2), float(base.yVal) + py * (arrow_width / 2)};

    arrow_verts[2].position = {float(base.xVal) - px * (arrow_width / 2), float(base.yVal) - py * (arrow_width / 2)};

    for (int i = 0; i < 3; i++)
        arrow_verts[i].color = f_color;

    SDL_RenderGeometry(renderer, nullptr, arrow_verts, 3, nullptr, 0);
}

void RenderSystem::buildCircleTextureCache()
{
    // Safety cleanup if re-initializing
    clearCachedCircleTextures();
    int max_radius = static_cast<int>(2 * MAX_RADIUS);
    circleTextureCache.reserve(max_radius);

    for (int radius = 1; radius <= max_radius; ++radius)
    {
        int diameter = radius * 2;

        SDL_Texture* tex =
            SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_TARGET, diameter, diameter);

        if (!tex)
        {
            // SDL_Log("Failed to create circle texture: %s", SDL_GetError());
            continue;
        }

        SDL_SetTextureBlendMode(tex, SDL_BLENDMODE_BLEND);

        // --- render circle into texture once ---
        SDL_SetRenderTarget(renderer, tex);

        // transparent background
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
        SDL_RenderClear(renderer);

        // white circle (IMPORTANT: color stays neutral for modulation later)
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);

        int r2 = radius * radius;

        for (int y = 0; y < diameter; ++y)
        {
            int dy = y - radius;

            for (int x = 0; x < diameter; ++x)
            {
                int dx = x - radius;

                if (dx * dx + dy * dy <= r2)
                {
                    SDL_RenderPoint(renderer, x, y);
                }
            }
        }

        SDL_SetRenderTarget(renderer, nullptr);

        circleTextureCache[radius] = tex;
    }
}
// --------- RENDER UI ELEMENTS METHOD --------- //

void RenderSystem::renderUIElements(UIState& UIState,
                                    const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope)
{
    if (!UIState.getAllUIVisibility())
    {
        return;
    }

    for (auto& [UI_element_ID, UI_element_ptr] : allUIElementsInScope)
    {
        if (UI_element_ID == UIElementIdentifier::PLAY_GAME_BUTTON_INDEX)
        {
            UI_element_ptr->renderMe(renderer, UIState, UIFontTitle);
        }
        else
        {
            UI_element_ptr->renderMe(renderer, UIState, UIFontRegular);
        }
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

void RenderSystem::clearCachedCircleTextures()
{
    for (auto& tex : circleTextureCache)
    {
        SDL_DestroyTexture(tex);
    }
    circleTextureCache.clear();
}

// --------- TWINKLING STAR METHODS --------- //

void RenderSystem::createStarField(int numStars)
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> alphaDist(0.3f, 1.0f);
    std::uniform_real_distribution<float> speedDist(1.0f, 3.0f);
    std::uniform_int_distribution<int> xDist(0, SCREEN_WIDTH);
    std::uniform_int_distribution<int> yDist(0, SCREEN_HEIGHT);

    for (int i = 0; i < numStars; ++i)
    {
        SDL_Texture* tex = twinklingStarTextures[rand() % twinklingStarTextures.size()]; // pick texture

        float w, h;
        SDL_GetTextureSize(tex, &w, &h);

        TwinklingStar star;
        star.texture = tex;
        star.dstRect = {float(xDist(rng)), float(yDist(rng)), w, h};
        star.baseAlpha = alphaDist(rng);
        star.twinkleSpeed = speedDist(rng);
        star.currentAlpha = star.baseAlpha;

        twinklingStars.push_back(star);
    }
}

void RenderSystem::updateStars()
{
    float time = SDL_GetTicks() / 1000.0f; // seconds
    for (auto& star : twinklingStars)
    {
        star.currentAlpha = star.baseAlpha + 0.3f * sinf(time * star.twinkleSpeed);
        star.currentAlpha = std::clamp(star.currentAlpha, 0.0f, 1.0f);
    }
}

void RenderSystem::renderStars()
{
    for (auto& star : twinklingStars)
    {
        Uint8 alpha = static_cast<Uint8>(star.currentAlpha * 255);
        SDL_SetTextureAlphaMod(star.texture, alpha);
        SDL_RenderTexture(renderer, star.texture, nullptr, &star.dstRect);
    }
}

void RenderSystem::createStarTextures()
{
    const int maxStarRadius = 4;
    const int minStarRadius = 1;

    for (int r = minStarRadius; r <= maxStarRadius; ++r)
    {
        SDL_Texture* sourceTex = circleTextureCache[r];
        float w, h;
        SDL_GetTextureSize(sourceTex, &w, &h);

        // Create a NEW unique texture for this star size
        SDL_Texture* starTex =
            SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_TARGET, (int)w, (int)h);

        // Ensure transparency is enabled for the new texture
        SDL_SetTextureBlendMode(starTex, SDL_BLENDMODE_BLEND);

        // Copy the circle to the new star texture
        SDL_SetRenderTarget(renderer, starTex);
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0); // Clear with transparent
        SDL_RenderClear(renderer);
        SDL_RenderTexture(renderer, sourceTex, nullptr, nullptr);

        // Reset render target to the screen
        SDL_SetRenderTarget(renderer, nullptr);

        // Now this texture belongs ONLY to the star system
        twinklingStarTextures.push_back(starTex);
    }
}