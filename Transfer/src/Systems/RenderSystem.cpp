// File: Transfer/src/Systems/RenderSystem.cpp

#include "Systems/RenderSystem.h"

// Constructor: Initializes SDL Window and Renderer
RenderSystem::RenderSystem()
{
    // Add audio later
    SDL_InitSubSystem(SDL_INIT_VIDEO);
    TTF_Init();
    SDL_SetHint(SDL_HINT_RENDER_DRIVER, "opengl");
    int desired_x_resolution = SCREEN_WIDTH;
    int desired_y_resolution = SCREEN_HEIGHT;
    int no_flags = 0;

    // Assign window pointer to instance
    window = SDL_CreateWindow("Transfer", desired_x_resolution, desired_y_resolution, no_flags);

    // Assign renderer pointer to instance. No specific driver, so nullptr.
    renderer = SDL_CreateRenderer(window, nullptr);
    // SDL_SetRenderVSync(renderer, 1);
    if (TTF_Init() == false)
    {
        SDL_Log("Failed to initialize SDL_ttf: %s", SDL_GetError());
        return;
    }
    std::string fontPathRegular = Utilities::GetResourcePath("Fonts/SpaceMono-Regular.ttf");
    std::string fontPathTitle = Utilities::GetResourcePath("Fonts/SpaceMono-Bold.ttf");

    UIFontRegular = TTF_OpenFont(fontPathRegular.c_str(), 18);
    UIFontTitle = TTF_OpenFont(fontPathTitle.c_str(), 24);
    buildCircleTextureCache();
    createStarTextures();
    createStarField(STAR_NUM);
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
}

// --------- RENDER FULL FRAME METHOD --------- //

void RenderSystem::RenderFullFrame(GameState& gameState, UIState& UIState,
                                   const std::unordered_map<UIElementIdentifier, UIElement*>& allUIElementsInScope)
{
    Uint64 start = SDL_GetPerformanceCounter();
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // black background
    SDL_RenderClear(renderer);

    // render starry background
    updateStars();
    renderStars();

    // Render all preview bodies
    renderPreviewBodies(UIState);
    // Render all the bodies (particles)
    renderBodies(gameState);

    renderUIElements(UIState, allUIElementsInScope);

    // Display the frame
    Uint64 beforePresent = SDL_GetPerformanceCounter();
    SDL_RenderPresent(renderer);
    Uint64 end = SDL_GetPerformanceCounter();
    double render_ms = (beforePresent - start) * 1000.0 / SDL_GetPerformanceFrequency();
    double present_ms = (end - beforePresent) * 1000.0 / SDL_GetPerformanceFrequency();

    // SDL_Log("PRESENT STALL: %f ms", present_ms);
    // SDL_Log("Render STALL: %f ms", render_ms);
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

    // 1. Group points by color to minimize state changes
    // Using a map of SDL_Color to a vector of points
    // Note: If you have a fixed set of colors, a fixed-size array is even faster
    std::map<Uint64, std::vector<SDL_FPoint>> coloredBatches;

    for (const auto& particle : particles)
    {
        if (!particle.visible)
            continue;

        // Linear interpolation for smooth movement
        float render_x = particle.previousPosition.xVal * (1.0f - alpha) + particle.position.xVal * alpha;
        float render_y = particle.previousPosition.yVal * (1.0f - alpha) + particle.position.yVal * alpha;

        SDL_Color c = getColorForProperty(particle);

        // Pack color into a Uint32 for easy mapping
        Uint32 colorKey = (c.r << 24) | (c.g << 16) | (c.b << 8) | c.a;
        coloredBatches[colorKey].push_back({render_x, render_y});
    }

    // 2. Render each batch in a single call
    for (auto& [colorKey, points] : coloredBatches)
    {
        Uint8 r = (colorKey >> 24) & 0xFF;
        Uint8 g = (colorKey >> 16) & 0xFF;
        Uint8 b = (colorKey >> 8) & 0xFF;
        Uint8 a = colorKey & 0xFF;

        SDL_SetRenderDrawColor(renderer, r, g, b, a);
        SDL_RenderPoints(renderer, points.data(), (int)points.size());
    }

    // 3. Render Macro Bodies (Keep these as textures since they are few and large)
    for (auto& body : gameState.getMacroBodies())
    {
        if (!body.visible)
            continue;

        SDL_Color color = getColorForProperty(body);
        SDL_Texture* tex = circleTextureCache[static_cast<int>(body.radius)];
        SDL_SetTextureColorMod(tex, color.r, color.g, color.b);
        SDL_SetTextureAlphaMod(tex, color.a);

        float render_x = body.previousPosition.xVal * (1.0f - alpha) + body.position.xVal * alpha;
        float render_y = body.previousPosition.yVal * (1.0f - alpha) + body.position.yVal * alpha;
        float r = static_cast<float>(body.radius);
        SDL_FRect dstRect = {render_x - r, render_y - r, r * 2, r * 2};

        SDL_RenderTexture(renderer, tex, nullptr, &dstRect);
    }
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
    int max_radius = static_cast<int>(MAX_RADIUS);
    circleTextureCache.reserve(max_radius);

    for (int radius = 1; radius <= max_radius; ++radius)
    {
        int diameter = radius * 2;

        SDL_Texture* tex =
            SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_TARGET, diameter, diameter);

        if (!tex)
        {
            SDL_Log("Failed to create circle texture: %s", SDL_GetError());
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
        return;

    for (auto& pair : allUIElementsInScope)
    {
        if (pair.first == PLAY_GAME_BUTTON_INDEX)
        {
            pair.second->renderMe(renderer, UIState, UIFontTitle);
        }
        else
            pair.second->renderMe(renderer, UIState, UIFontRegular);
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
    // if (absMass > MAX_MASS)
    //     return SDL_Color{0, 0, 0, 255};

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

// Error in here somewhere. funky
// {
//     double mass = body.mass;
//     double absMass = std::abs(mass);

//     // 1. Event Horizon check
//     if (absMass > MAX_MASS)
//         return SDL_Color{0, 0, 0, 255};

//     // 2. High-Compression Scaling
//     // We use a very low exponent (0.07) so that small masses
//     // stay white/pale longer before turning deep blue/red.
//     static const double exponent = 0.07;
//     static const double maxScaled = std::pow(MAX_MASS, exponent);

//     // t goes from 0.0 (massless) to 1.0 (at MAX_MASS)
//     double t = std::clamp(std::pow(absMass, exponent) / maxScaled, 0.0, 1.0);

//     // 3. Start at White (255, 255, 255)
//     // As 't' increases, we subtract color from the channels we DON'T want.
//     Uint8 r, g, b;
//     Uint8 dropOff = static_cast<Uint8>(255 * t);

//     if (mass < 0)
//     {
//         // NEGATIVE: Skew toward RED
//         // We keep Red high, and pull Green/Blue down toward 0
//         r = 255;
//         g = 255 - dropOff;
//         b = 255 - dropOff;
//     }
//     else
//     {
//         // POSITIVE: Skew toward BLUE
//         // We keep Blue high, and pull Red/Green down toward 0
//         r = 255 - dropOff;
//         g = 255 - dropOff;
//         b = 255;
//     }

//     Uint8 opacity = (body.isMacroGhost || !body.isCollidable) ? 175 : 255;
//     return SDL_Color{r, g, b, opacity};
// }
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

// void RenderSystem::createStarTextures()
// {
//     const int maxStarRadius = 3;
//     const int minStarRadius = 1;

//     for (int r = minStarRadius; r <= maxStarRadius; ++r)
//     {
//         SDL_Color color = {255, 255, 255, 255};
//         SDL_Texture* tex = circleTextureCache[r];
//         SDL_SetTextureColorMod(tex, color.r, color.g, color.b);
//         SDL_SetTextureAlphaMod(tex, color.a);
//         twinklingStarTextures.push_back(tex);
//     }
// }

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