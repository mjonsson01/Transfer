// Build (from the repo root, macOS; any .cpp under Transfer/src except main.cpp):
//   clang++ -std=c++20 -O1 -w -ITransfer/src -IThirdParty/SDL3/include -IThirdParty/SDL3_ttf/include Notes/drafts/Harnesses/ghost_harness.cpp $(find Transfer/src -name '*.cpp' ! -name main.cpp) -LThirdParty/SDL3/lib -LThirdParty/SDL3_ttf/lib -lSDL3 -lSDL3_ttf -Wl,-rpath,$PWD/ThirdParty/SDL3/lib -Wl,-rpath,$PWD/ThirdParty/SDL3_ttf/lib -o /tmp/ghost && /tmp/ghost
// Two overlapping particles closing head-on, mass ratio 100 (>= 8). Ghost = they pass through; fixed = they bounce.
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "Systems/PhysicsSystem.hpp"
#include <cstdio>

static GravitationalBody makeParticle(double x, double vx, double mass)
{
    GravitationalBody p;
    p.mass = mass; p.invMass = 1.0 / mass; p.radius = 2.0;
    p.position = {x, 0.0}; p.previousPosition = p.position; p.velocity = {vx, 0.0};
    p.isParticle = true; p.isFragment = true; p.isAccretable = true; p.isCollidable = true;
    return p;
}

int main()
{
    GameState game_state; UIState ui_state; PhysicsSystem physics;
    game_state.getParticlesMutable().push_back(makeParticle(0.0, +100.0, 1.0));   // light, moving right
    game_state.getParticlesMutable().push_back(makeParticle(3.5, -100.0, 100.0)); // heavy, moving left
    for (int tick = 0; tick < 120; ++tick)
        physics.UpdateSystemFrame(game_state, ui_state);
    auto& ps = game_state.getParticles();
    printf("count=%zu\n", ps.size());
    for (auto& p : ps) printf("mass=%6.1f x=%8.2f vx=%8.2f\n", p.mass, p.position.x_val, p.velocity.x_val);
    bool ghosted = ps.size() == 2 && ps[0].position.x_val > ps[1].position.x_val;
    printf(ghosted ? "RESULT: GHOSTED THROUGH\n" : "RESULT: BOUNCED\n");
}
