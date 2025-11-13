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
    // SDL_Log("Loaded font from path: %s", fontPath.c_str());
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
    // TTF_Quit(); // Handled at the Game level.
}

void RenderSystem::CleanUp()
{
    clearCachedCircleTextures();
}

// Top level rendering method, calls all subordinates in order
void RenderSystem::RenderFullFrame(GameState& state, UIState& UIState)
{
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // black background
    SDL_RenderClear(renderer);

    // render starry background
    updateStars();
    renderStars();
    // Render all the bodies
    RenderSystem::renderBodies(state);

    // Render any input artifacts (like drag lines) TBI
    // RenderSystem::renderInputArtifacts(state);
    // Add new rendering system functions later

    // Render UI Elements
    uiSystem.RenderUIElements(renderer, UIState, UIFont);

    // Display the frame
    SDL_RenderPresent(renderer);
}

// Renders all gravitational bodies in the game state 
void RenderSystem::renderBodies(GameState& state)
{
    float alpha = state.getAlpha();
    for (auto& body : state.getBodies()) {
        SDL_Color color = getColorForMass(body.getMass());
        SDL_Texture* tex = getCircleTexture(body.getRadius(), color);
        Vector2D currPosition = body.getPosition();
        Vector2D prevPosition = body.getPrevPosition();
        float renderX = prevPosition.x_val * (1.0f - alpha) + currPosition.x_val * alpha;
        float renderY = prevPosition.y_val * (1.0f - alpha) + currPosition.y_val * alpha;

        SDL_FRect dstRect = {
            renderX - body.getRadius(),
            renderY - body.getRadius(),
            body.getRadius() * 2,
            body.getRadius() * 2
        };

        SDL_RenderTexture(renderer, tex, nullptr, &dstRect);
    }
}

// Renders user input artifacts like drag lines. TBI
void RenderSystem::renderInputArtifacts(GameState& state)
{

}

// Smooth interpolation color lookup function
SDL_Color RenderSystem::getColorForMass(double mass)
{
    // Simple mapping: lighter masses are blue, heavier masses are red
    Uint8 r = static_cast<Uint8>(std::min((mass / MAX_MASS) * 255, 255.0)); // assuming max mass of 1e25 for scaling
    Uint8 g = 0;
    Uint8 b = static_cast<Uint8>(255 - r);
    Uint8 a = 255; // fully opaque

    return SDL_Color{ r, g, b, a };
}

// Helper to render the Frame Rate Counter if enabled.
// void RenderSystem::renderFrameRateCounter(float fps)
// {
//     // Convert FPS to string
//     std::string fpsText = "FPS: " + std::to_string(static_cast<int>(fps));
//     SDL_Surface* textSurface = TTF_RenderText_Blended(FPSFont, fpsText.c_str(), fpsText.length(), ColorLibrary::White);
//     if (!textSurface) {
//         SDL_Log("Text surface creation failed: %s", SDL_GetError());
//         return;
//     }

//     SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
//     if (!textTexture) {
//         SDL_Log("Text texture creation failed: %s", SDL_GetError());
//         return;
//     }
    
//     SDL_FRect dstRect = {10.0f, 10.0f, static_cast<float>(textSurface->w), static_cast<float>(textSurface->h)};
//     SDL_RenderTexture(renderer, textTexture, nullptr, &dstRect);
//     SDL_DestroySurface(textSurface);
//     SDL_DestroyTexture(textTexture);
// }

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
    // Bresenham's circle algorithm just for the memes.
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

void RenderSystem::createStarField(int numStars)
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> alphaDist(0.3f, 1.0f);
    std::uniform_real_distribution<float> speedDist(1.0f, 3.0f);
    std::uniform_int_distribution<int> xDist(0, SCREEN_WIDTH);
    std::uniform_int_distribution<int> yDist(0, SCREEN_HEIGHT);

    for (int i = 0; i < numStars; ++i)
    {
        SDL_Texture* tex = starTextures[rand() % starTextures.size()]; // pick texture

        float w, h;
        SDL_GetTextureSize(tex, &w, &h);

        Star star;
        star.texture = tex;
        star.dstRect = { float(xDist(rng)), float(yDist(rng)), w, h };
        star.baseAlpha = alphaDist(rng);
        star.twinkleSpeed = speedDist(rng);
        star.currentAlpha = star.baseAlpha;

        stars.push_back(star);
    }
}

void RenderSystem::updateStars()
{
    float time = SDL_GetTicks() / 1000.0f; // seconds
    for (auto& star : stars)
    {
        star.currentAlpha = star.baseAlpha + 0.3f * sinf(time * star.twinkleSpeed);
        star.currentAlpha = std::clamp(star.currentAlpha, 0.0f, 1.0f);
    }
}

void RenderSystem::renderStars()
{
    for (auto& star : stars)
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
        starTextures.push_back(tex);
    }
}
