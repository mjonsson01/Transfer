// Build (from the repo root, macOS):
//   clang++ -std=c++20 -O1 -w -ITransfer/src -IThirdParty/SDL3/include -IThirdParty/SDL3_ttf/include Notes/drafts/Harnesses/ship_physics_harness.cpp $(find Transfer/src -name '*.cpp' ! -name main.cpp) -LThirdParty/SDL3/lib -LThirdParty/SDL3_ttf/lib -lSDL3 -lSDL3_ttf -Wl,-rpath,$PWD/ThirdParty/SDL3/lib -Wl,-rpath,$PWD/ThirdParty/SDL3_ttf/lib -o /tmp/ship_physics_harness && /tmp/ship_physics_harness
// Ship physics: thrust feel, and a circular orbit around one planet (one-way gravity + Verlet).
#define private public // harness only: lets the test place the ship directly
#include "Player/Starship.hpp"
#undef private
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "Systems/PhysicsSystem.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
int main()
{
    { // thrust: hold W for 120 ticks (1 s) from rest, no planets
        GameState gs; UIState ui; PhysicsSystem p;
        DEPRECATED_InputState& in = ui.getMutableDEPRECATED_InputState();
        in.isRequestingThrust = true; in.positiveThrust = true;
        for (int t = 0; t < 120; ++t) p.UpdateSystemFrame(gs, ui);
        Starship& s = gs.getPlayerMutable().starship;
        printf("thrust 1 s from rest: speed %.1f px/s (expect ~1200), direction (%.3f, %.3f) (expect (0, -1): nose up)\n",
               s.velocity().magnitude(), s.velocity().normalize().x_val, s.velocity().normalize().y_val);
    }
    { // orbit: planet M at the origin, ship centre at r = 1000 to the right, circular speed straight down
        GameState gs; UIState ui; PhysicsSystem p;
        GravitationalBody planet; planet.mass = 1e9; planet.invMass = 1e-9; planet.radius = 100; planet.isMacro = true;
        planet.position = {0, 0}; planet.previousPosition = planet.position;
        gs.getMacroBodiesMutable().push_back(planet);
        Starship& s = gs.getPlayerMutable().starship;
        double r = 1000.0, eps = (s.halfSize() + planet.radius) / 2.0;
        double v = std::sqrt(planet.mass * r * r / std::pow(r * r + eps * eps, 1.5)); // circular speed with softening
        s.m_position = DynamoEngine::Vector2D(r, 0) - DynamoEngine::Vector2D(s.halfSize(), s.halfSize());
        s.m_prev_position = s.m_position; s.m_velocity = {0, v};
        double period = 2 * M_PI * r / v; int ticks = int(2 * period * 120);
        double rmin = 1e18, rmax = 0;
        for (int t = 0; t < ticks; ++t)
        {
            p.UpdateSystemFrame(gs, ui);
            double d = s.center().magnitude(); rmin = std::min(rmin, d); rmax = std::max(rmax, d);
        }
        printf("orbit: speed %.1f px/s, 2 orbits = %d ticks, distance stayed in [%.2f, %.2f] (start 1000) -> %.3f%% wobble\n",
               v, ticks, rmin, rmax, 100 * (rmax - rmin) / r);
        printf("planet moved to (%.4f, %.4f) (two-way: a tiny wobble, ship/planet mass = 1e-6)\n",
               gs.getMacroBodies()[0].position.x_val, gs.getMacroBodies()[0].position.y_val);
    }
    { // momentum: a LIGHT planet (mass 100000, only 100x the ship) and the ship pull on each other from rest
        GameState gs; UIState ui; PhysicsSystem p;
        GravitationalBody planet; planet.mass = 1e5; planet.invMass = 1e-5; planet.radius = 20; planet.isMacro = true;
        planet.position = {0, 0}; planet.previousPosition = planet.position;
        gs.getMacroBodiesMutable().push_back(planet);
        Starship& s = gs.getPlayerMutable().starship;
        s.m_position = DynamoEngine::Vector2D(400, 0) - DynamoEngine::Vector2D(s.halfSize(), s.halfSize());
        s.m_prev_position = s.m_position;
        for (int t = 0; t < 120; ++t) p.UpdateSystemFrame(gs, ui);
        const GravitationalBody& pl = gs.getMacroBodies()[0];
        DynamoEngine::Vector2D momentum = s.velocity() * s.mass() + pl.velocity * pl.mass;
        printf("after 1 s: ship v=(%.3f, %.3f) planet v=(%.5f, %.5f) total momentum=(%.2e, %.2e) (expect ~0: started at rest)\n",
               s.velocity().x_val, s.velocity().y_val, pl.velocity.x_val, pl.velocity.y_val, momentum.x_val, momentum.y_val);
    }
}
