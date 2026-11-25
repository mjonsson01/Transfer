// Build + run: see README.md in this folder.
// Ship physics: thrust feel, and a circular orbit around one planet (one-way gravity + Verlet).
#define private public // harness only: lets the test place the ship directly
#include "Player/Starship.hpp"
#undef private
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "Systems/PhysicsSystem.hpp"
#include "HarnessInput.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
int main()
{
    { // thrust: hold W for 120 ticks (1 s) from rest, no planets
        GameState gs; UIState ui; PhysicsSystem p;
        for (int t = 0; t < 120; ++t) p.advanceUniverseState(gs.getUniverseStateMutable(), shipControls(1, 0)); // W held
        Starship& s = gs.getUniverseStateMutable().ship;
        printf("W held 1 s from rest: throttle %.3f (expect 1.000), speed %.1f px/s (expect ~600: the throttle ramps "
               "0 -> 1 over that second, so the average is half of 1200), direction (%.3f, %.3f) (expect (0, -1))\n",
               s.throttle(), s.velocity().magnitude(), s.velocity().normalize().x_val, s.velocity().normalize().y_val);
        double speed_after_ramp = s.velocity().magnitude();
        for (int t = 0; t < 120; ++t) p.advanceUniverseState(gs.getUniverseStateMutable(), ShipControls()); // released
        printf("released for 1 s: throttle %.3f (expect 1.000: it holds), speed gained %.1f px/s (expect ~1200)\n",
               s.throttle(), s.velocity().magnitude() - speed_after_ramp);
        double speed_before_s = s.velocity().magnitude();
        for (int t = 0; t < 240; ++t) p.advanceUniverseState(gs.getUniverseStateMutable(), shipControls(-1, 0)); // S
        printf("S held 2 s: throttle %.3f (expect 0.000, never below), speed change %.1f px/s (expect ~+600 while it "
               "ramped down, then 0: no reverse thrust)\n",
               s.throttle(), s.velocity().magnitude() - speed_before_s);
    }
    { // orbit: planet M at the origin, ship centre at r = 1000 to the right, circular speed straight down
        GameState gs; UIState ui; PhysicsSystem p;
        GravitationalBody planet; planet.mass = 1e9; planet.invMass = 1e-9; planet.radius = 100; planet.isMacro = true;
        planet.position = {0, 0}; planet.previousPosition = planet.position;
        gs.getUniverseStateMutable().macro_bodies.push_back(planet);
        Starship& s = gs.getUniverseStateMutable().ship;
        double r = 1000.0, eps = (s.halfSize() + planet.radius) / 2.0;
        double v = std::sqrt(planet.mass * r * r / std::pow(r * r + eps * eps, 1.5)); // circular speed with softening
        s.m_position = DynamoEngine::Vector2D(r, 0) - DynamoEngine::Vector2D(s.halfSize(), s.halfSize());
        s.m_prev_position = s.m_position; s.m_velocity = {0, v};
        double period = 2 * M_PI * r / v; int ticks = int(2 * period * 120);
        double rmin = 1e18, rmax = 0;
        for (int t = 0; t < ticks; ++t)
        {
            p.advanceUniverseState(gs.getUniverseStateMutable(), ShipControls());
            double d = s.center().magnitude(); rmin = std::min(rmin, d); rmax = std::max(rmax, d);
        }
        printf("orbit: speed %.1f px/s, 2 orbits = %d ticks, distance stayed in [%.2f, %.2f] (start 1000) -> %.3f%% wobble\n",
               v, ticks, rmin, rmax, 100 * (rmax - rmin) / r);
        printf("planet moved to (%.4f, %.4f) (two-way: a tiny wobble, ship/planet mass = 1e-6)\n",
               gs.getUniverseState().macro_bodies[0].position.x_val, gs.getUniverseState().macro_bodies[0].position.y_val);
    }
    { // momentum: a LIGHT planet (mass 100000, only 100x the ship) and the ship pull on each other from rest
        GameState gs; UIState ui; PhysicsSystem p;
        GravitationalBody planet; planet.mass = 1e5; planet.invMass = 1e-5; planet.radius = 20; planet.isMacro = true;
        planet.position = {0, 0}; planet.previousPosition = planet.position;
        gs.getUniverseStateMutable().macro_bodies.push_back(planet);
        Starship& s = gs.getUniverseStateMutable().ship;
        s.m_position = DynamoEngine::Vector2D(400, 0) - DynamoEngine::Vector2D(s.halfSize(), s.halfSize());
        s.m_prev_position = s.m_position;
        for (int t = 0; t < 120; ++t) p.advanceUniverseState(gs.getUniverseStateMutable(), ShipControls());
        const GravitationalBody& pl = gs.getUniverseState().macro_bodies[0];
        DynamoEngine::Vector2D momentum = s.velocity() * s.mass() + pl.velocity * pl.mass;
        printf("after 1 s: ship v=(%.3f, %.3f) planet v=(%.5f, %.5f) total momentum=(%.2e, %.2e) (expect ~0: started at rest)\n",
               s.velocity().x_val, s.velocity().y_val, pl.velocity.x_val, pl.velocity.y_val, momentum.x_val, momentum.y_val);
    }
}
