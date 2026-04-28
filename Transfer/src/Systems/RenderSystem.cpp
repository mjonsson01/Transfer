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

void RenderSystem::RenderFullFrame(GameState& gameState, UIState& UIState, const std::vector<UIElement*>& allUIElementsInScope)
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
    // Render Particles
    for (auto& particle : gameState.getParticles())
    {
        if (!particle.visible)
        {
            continue; // don't render invisible particles
        }
        SDL_Color color = getColorForProperty(particle);
        SDL_Texture* tex = circleTextureCache[static_cast<int>(particle.radius)];
        SDL_SetTextureColorMod(tex, color.r, color.g, color.b);
        SDL_SetTextureAlphaMod(tex, color.a);

        Vector2D current_position = particle.position;
        Vector2D previous_position = particle.previousPosition;
        // Alpha Interpolation causes particle flickers for small particles.
        // Remove for now. Figure out dynamical fix later. Leave interpolation
        // on for slow mo.

        float render_x, render_y;
        if (particle.radius <= 1.0 && gameState.getToggleSlow() == false)
        {
            render_x = particle.position.xVal;
            render_y = particle.position.yVal;
            render_x = std::round(render_x);
            render_y = std::round(render_y);
        }
        else
        {
            render_x = previous_position.xVal * (1.0f - alpha) + current_position.xVal * alpha;
            render_y = previous_position.yVal * (1.0f - alpha) + current_position.yVal * alpha;
        }
        float r = static_cast<float>(particle.radius);
        SDL_FRect dstRect = {render_x - r, render_y - r, r * 2, r * 2};

        SDL_RenderTexture(renderer, tex, nullptr, &dstRect);
    }

    // Render Macro Bodies
    for (auto& body : gameState.getMacroBodies())
    {
        if (!body.visible)
        {
            continue;
        }
        SDL_Color color = getColorForProperty(body);
        SDL_Texture* tex = circleTextureCache[static_cast<int>(body.radius)];
        SDL_SetTextureColorMod(tex, color.r, color.g, color.b);
        SDL_SetTextureAlphaMod(tex, color.a);

        Vector2D current_position = body.position;
        Vector2D previous_position = body.previousPosition;
        float render_x = previous_position.xVal * (1.0f - alpha) + current_position.xVal * alpha;
        float render_y = previous_position.yVal * (1.0f - alpha) + current_position.yVal * alpha;
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

    // SDL_Log("Circle texture cache built: %d radii", max_radius);
}
// --------- RENDER UI ELEMENTS METHOD --------- //

void RenderSystem::renderUIElements(UIState& UIState, std::vector<UIElement*> allUIElementsInScope)
{
    if (!UIState.getAllUIVisibility())
        return;

    for (auto& element : allUIElementsInScope)
    {
        if (element->getUIElementType() == PLAY_GAME_BUTTON_INDEX)
        {
            element->renderMe(renderer, UIState, UIFontTitle);
        }
        else
            element->renderMe(renderer, UIState, UIFontRegular);
    }
}

// Smooth interpolation color lookup function

// --------- RENDER UTILITY HELPERS --------- //
SDL_Color RenderSystem::getColorForProperty(const GravitationalBody& body)
// NEED TO OPTIMIZE OUT CONSTANT ACCESS
{
    // 1. Use Logarithmic scaling so small changes in light bodies are visible
    // mass + 1.0 avoids log(0)
    double logMass = std::log10(body.mass + 1.0);
    double logMax = std::log10(MAX_MASS + 1.0);
    double t = std::clamp(logMass / logMax, 0.0, 1.0);
    Uint8 opacity = 255;

    // Ghost mode render
    if (body.isMacroGhost || !body.isCollidable)
    {
        opacity = 175;
    }
    // 2. Multi-stop gradient (Blue -> Cyan -> Green -> Yellow -> Red)
    Uint8 r, g, b;

    if (t < 0.25)
    { // Blue to Cyan
        r = 0;
        g = static_cast<Uint8>(t * 4.0 * 255);
        b = 255;
    }
    else if (t < 0.5)
    { // Cyan to Green
        r = 0;
        g = 255;
        b = static_cast<Uint8>((0.5 - t) * 4.0 * 255);
    }
    else if (t < 0.75)
    { // Green to Yellow
        r = static_cast<Uint8>((t - 0.5) * 4.0 * 255);
        g = 255;
        b = 0;
    }
    else
    { // Yellow to Red
        r = 255;
        g = static_cast<Uint8>((1.0 - t) * 4.0 * 255);
        b = 0;
    }

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
    const int maxStarRadius = 3;
    const int minStarRadius = 1;

    for (int r = minStarRadius; r <= maxStarRadius; ++r)
    {
        SDL_Color color = {255, 255, 255, 255};
        SDL_Texture* tex = circleTextureCache[r];
        SDL_SetTextureColorMod(tex, color.r, color.g, color.b);
        SDL_SetTextureAlphaMod(tex, color.a);
        twinklingStarTextures.push_back(tex);
    }
}
