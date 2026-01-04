// File: Transfer/src/Systems/RenderSystem.cpp

#include "Systems/RenderSystem.h"

// Constructor: Initializes SDL Window and Renderer
RenderSystem::RenderSystem()
{
    // Add audio later
    SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS);
    TTF_Init();

    int desired_x_resolution = SCREEN_WIDTH;
    int desired_y_resolution = SCREEN_HEIGHT;
    int no_flags = 0;

    // Assign window pointer to instance
    window = SDL_CreateWindow("Transfer", desired_x_resolution, desired_y_resolution, no_flags);

    // Assign renderer pointer to instance. No specific driver, so nullptr.
    renderer = SDL_CreateRenderer(window, nullptr);
    if (TTF_Init() == false) {
        SDL_Log("Failed to initialize SDL_ttf: %s", SDL_GetError());
        return;
    }
    const char* basePath = SDL_GetBasePath(); // returns the folder where the executable lives
    std::string fontPath = std::string(basePath) + "Assets/Fonts/SpaceMono-Bold.ttf";

    UIFont = TTF_OpenFont(fontPath.c_str(), 18);
    if (!UIFont) {
        SDL_Log("Failed to load font: %s", SDL_GetError());
        return;
    }
    createStarTextures();
    createStarField(STAR_NUM);   
}

// Destructor: Cleans up SDL Window and Renderer
RenderSystem::~RenderSystem()
{

    if (renderer) {
        SDL_DestroyRenderer(renderer);
        renderer = nullptr;
    }
    if (window) {
        SDL_DestroyWindow(window);
        window = nullptr;
    }
    // TTF_Quit() and SDL_QUIT() handled at the Game level
}

// --------- CLEANUP METHOD --------- //
void RenderSystem::CleanUp()
{
    clearCachedCircleTextures();
}


// --------- TOTAL ENERGY CALCULATION METHOD --------- //

void RenderSystem::RenderFullFrame(GameState& state, UIState& UIState)
{
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // black background
    SDL_RenderClear(renderer);

    // render starry background
    updateStars();
    renderStars();

    // Render all the bodies (particles)
    RenderSystem::renderBodies(state);

    // Render UI Elements
    uiSystem.RenderUIElements(renderer, UIState, UIFont);

    // Display the frame
    SDL_RenderPresent(renderer);
}

// --------- RENDER GRAVITATIONAL BODIES METHOD --------- //
void RenderSystem::renderBodies(GameState& state)
{
    float alpha = state.getAlpha();
    // Render Particles
    for (auto& particle : state.getParticles()) {
        SDL_Color color = getColorForMass(particle.mass);
        SDL_Texture* tex = getCircleTexture(static_cast<int>(particle.radius), color);
        Vector2D current_position = particle.position;
        Vector2D previous_position = particle.previousPosition;
        // Alpha Interpolation causes particle flickers for small particles. Remove for now. Figure out dynamical fix later. 
        // Leave interpolation on for slow mo.
        float render_x, render_y;
        if (particle.radius <= 1.0 && state.getToggleSlow() == false) {
            render_x = particle.position.xVal;
            render_y = particle.position.yVal;
            render_x = std::round(render_x);
            render_y = std::round(render_y);
        }
        else {
            render_x = previous_position.xVal * (1.0f - alpha) + current_position.xVal * alpha;
            render_y = previous_position.yVal * (1.0f - alpha) + current_position.yVal * alpha;
        }
        float r = static_cast<float>(particle.radius);
        SDL_FRect dstRect = {
            render_x - r,
            render_y - r,
            r * 2,
            r * 2
        };

        SDL_RenderTexture(renderer, tex, nullptr, &dstRect);
    }

    // Render Macro Bodies
    for (auto& body : state.getMacroBodies()) {
        SDL_Color color = getColorForMass(body.mass);
        SDL_Texture* tex = getCircleTexture(static_cast<int>(body.radius), color);
        Vector2D current_position = body.position;
        Vector2D previous_position = body.previousPosition;
        float render_x = previous_position.xVal * (1.0f - alpha) + current_position.xVal * alpha;
        float render_y = previous_position.yVal * (1.0f - alpha) + current_position.yVal * alpha;
        float r = static_cast<float>(body.radius);
        SDL_FRect dstRect = {
            render_x - r,
            render_y - r,
            r * 2,
            r * 2
        };

        SDL_RenderTexture(renderer, tex, nullptr, &dstRect);
    }
    
}

// --------- RENDER INPUT ARTIFACTS METHOD --------- //
// Renders user input artifacts like drag lines. TBI
void RenderSystem::renderInputArtifacts(GameState& state)
{

}

// Smooth interpolation color lookup function

// --------- RENDER UTILITY HELPERS --------- //
SDL_Color RenderSystem::getColorForMass(double mass)
{
    // Simple mapping: lighter masses are blue, heavier masses are red
    Uint8 r = static_cast<Uint8>(std::min((mass / MAX_MASS) * 255, 255.0)); // assuming max mass of 1e25 for scaling
    Uint8 g = 0;
    Uint8 b = static_cast<Uint8>(255 - r);
    Uint8 a = 255; // fully opaque
    // if (mass >= MAX_MASS-1000.0) return ColorLibrary::Black;
    // if (mass <= MAX_MASS/10000.0) return ColorLibrary::White;
    
    return SDL_Color{ r, g, b, a };
}

// Paired below

// static SDL_Color HSVtoRGB(double h, double s, double v, Uint8 a = 255)
// {
//     h = std::fmod(h, 360.0);
//     if (h < 0) h += 360.0;

//     double c = v * s;
//     double x = c * (1.0 - std::fabs(std::fmod(h / 60.0, 2.0) - 1.0));
//     double m = v - c;

//     double r1=0, g1=0, b1=0;
//     if      (h < 60)  { r1 = c; g1 = x; b1 = 0; }
//     else if (h < 120) { r1 = x; g1 = c; b1 = 0; }
//     else if (h < 180) { r1 = 0; g1 = c; b1 = x; }
//     else if (h < 240) { r1 = 0; g1 = x; b1 = c; }
//     else if (h < 300) { r1 = x; g1 = 0; b1 = c; }
//     else              { r1 = c; g1 = 0; b1 = x; }

//     Uint8 r = (Uint8)std::clamp((r1 + m) * 255.0, 0.0, 255.0);
//     Uint8 g = (Uint8)std::clamp((g1 + m) * 255.0, 0.0, 255.0);
//     Uint8 b = (Uint8)std::clamp((b1 + m) * 255.0, 0.0, 255.0);
//     return SDL_Color{ r, g, b, a };
// }

// SDL_Color RenderSystem::getColorForMass(double mass)
// {
//     // Normalize mass to [0,1]
//     double t = std::clamp(mass / MAX_MASS, 0.0, 1.0);

//     // Optional: log-scale for better contrast if masses span many orders of magnitude
//     // double t = std::clamp(std::log10(mass / MIN_MASS) / std::log10(MAX_MASS / MIN_MASS), 0.0, 1.0);

//     // Map t to hue. 0..360 sweeps full rainbow.
//     // If you want low mass blue and high mass red, use 240 -> 0:
//     double hue = 240.0 * (1.0 - t); // blue (240) -> red (0)

//     double sat = 1.0;
//     double val = 1.0;
//     return HSVtoRGB(hue, sat, val, 255);
// }


// --------- CIRCLE TEXTURE CACHE METHODS --------- //
// Helper to correctly destroy the circle texture cache.
void RenderSystem::clearCachedCircleTextures()
{
    for (auto& pair : circleTextureCache) {
        SDL_DestroyTexture(pair.second);
    }
    circleTextureCache.clear();
}

// Helper to get the texture for a given radius/color pair
SDL_Texture* RenderSystem::getCircleTexture(int radius, SDL_Color color)
{
    CircleKey key{ radius, color };

    auto it = circleTextureCache.find(key);
    if (it != circleTextureCache.end())
        return it->second;

    SDL_Texture* tex = createCircleTexture(radius, color);
    circleTextureCache[key] = tex;
    return tex;
}

// Helper to create a circle (grav body) texture and places in the cache.
SDL_Texture* RenderSystem::createCircleTexture(int radius, SDL_Color color)
{
    int diameter = radius * 2;
    SDL_Texture* tex = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_TARGET, diameter, diameter);
    SDL_SetTextureBlendMode(tex, SDL_BLENDMODE_BLEND);

    SDL_SetRenderTarget(renderer, tex);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0); // transparent background
    SDL_RenderClear(renderer);

    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);

    // Simple naive n^2 fill circle -- fast enough since only executed once. May further optimize to
    // Bresenham's circle algorithm later
    for (int w = 0; w < diameter; w++)
    {
        for (int h = 0; h < diameter; h++)
        {
            int dx = radius - w;
            int dy = radius - h;
            if ((dx*dx + dy*dy) <= (radius*radius))
                SDL_RenderPoint(renderer, w, h);
        }
    }
    SDL_SetRenderTarget(renderer, nullptr);
    return tex;
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
        star.dstRect = { float(xDist(rng)), float(yDist(rng)), w, h };
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
        SDL_Texture* tex = createCircleTexture(r, color);
        twinklingStarTextures.push_back(tex);
    }
}
