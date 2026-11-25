// Build + run: see README.md in this folder.
// Ship collisions: bounce off a planet, shatter on a hard hit, plow through debris, rest on a planet, momentum.
#define private public // harness only: place the ship directly
#include "Player/Starship.hpp"
#undef private
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "DynamoEngine/Physics/Collision2D.hpp"
#include "Systems/PhysicsSystem.hpp"
#include "HarnessInput.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
using V = DynamoEngine::Vector2D;

static GravitationalBody planet(V pos, double mass, double radius)
{
    GravitationalBody b; b.mass = mass; b.invMass = 1.0 / mass; b.radius = radius; b.isMacro = true;
    b.isCollidable = true; b.isShatterable = true; b.isAccretable = true; b.position = pos; b.previousPosition = pos;
    return b;
}
static void placeShip(Starship& s, V center, V velocity) { s.m_position = center - V(s.halfSize(), s.halfSize()); s.m_prev_position = s.m_position; s.m_velocity = velocity; }
static double worstOverlap(const GameState& gs)
{
    double worst = 0; auto poly = gs.getUniverseState().ship.collisionPolygon();
    for (auto& b : gs.getUniverseState().macro_bodies) { auto c = DynamoEngine::circleVsConvexPolygon(b.position, b.radius, poly); if (c.touching) worst = std::max(worst, c.depth); }
    for (auto& b : gs.getUniverseState().particles) { auto c = DynamoEngine::circleVsConvexPolygon(b.position, b.radius, poly); if (c.touching) worst = std::max(worst, c.depth); }
    return worst;
}

int main()
{
    { // 0. energy audit: ship (1000) head-on into one particle (mass 3000) at rest, no planets = no gravity
        GameState gs; UIState ui; PhysicsSystem p;
        GravitationalBody q; q.mass = 3000; q.invMass = 1.0 / 3000; q.radius = 5; q.isParticle = true; q.isCollidable = true;
        q.position = {0, -120}; q.previousPosition = q.position; gs.getUniverseStateMutable().particles.push_back(q);
        Starship& s = gs.getUniverseStateMutable().ship; placeShip(s, {0, 0}, {0, -400});
        auto kinetic = [&]() { double k = 0.5 * s.mass() * s.m_velocity.squareMagnitude(); for (auto& b : gs.getUniverseState().particles) k += 0.5 * b.mass * b.velocity.squareMagnitude(); return k; };
        double before = kinetic(), reported = 0;
        for (int t = 0; t < 120; ++t) { p.advanceUniverseState(gs.getUniverseStateMutable(), ShipControls()); reported += s.lastImpactEnergy(); }
        double mu = 1.0 / (1.0 / 1000 + 1.0 / 3000), expect = 0.5 * mu * 400 * 400 * (1 - 0.2 * 0.2);
        printf("0 energy: kinetic before %.0f after %.0f -> lost %.0f; reported %.0f; formula 1/2 mu v^2 (1-e^2) = %.0f\n",
               before, kinetic(), before - kinetic(), reported, expect);
    }
    { // 1. slow head-on hit on a heavy planet (G=1, mass 1e6: weak pull at this range), below shatter speed
        GameState gs; UIState ui; PhysicsSystem p;
        gs.getUniverseStateMutable().macro_bodies.push_back(planet({0, -300}, 1e6, 100)); // planet above the ship's nose
        Starship& s = gs.getUniverseStateMutable().ship; placeShip(s, {0, 0}, {0, -200});
        double vin = s.m_velocity.y_val, before = 0; int hit_tick = -1;
        for (int t = 0; t < 240; ++t) { before = s.m_velocity.y_val; p.advanceUniverseState(gs.getUniverseStateMutable(), ShipControls()); if (hit_tick < 0 && s.m_velocity.y_val > 0) hit_tick = t; }
        printf("1 bounce: in %.0f px/s, out %.1f px/s (expect about +40 or a bit less: 0.2 x 200, reversed, minus the planet pulling it back), planets left %zu, overlap now %.3f\n",
               vin, s.m_velocity.y_val, gs.getUniverseState().macro_bodies.size(), worstOverlap(gs));
    }
    { // 2. fast hit: no ship shattering any more -> the planet survives, the ship bounces
        GameState gs; UIState ui; PhysicsSystem p;
        gs.getUniverseStateMutable().macro_bodies.push_back(planet({0, -300}, 1e6, 100));
        Starship& s = gs.getUniverseStateMutable().ship; placeShip(s, {0, 0}, {0, -600});
        for (int t = 0; t < 60; ++t) p.advanceUniverseState(gs.getUniverseStateMutable(), ShipControls());
        printf("2 fast hit: planets left %zu (expect 1), particles %zu (expect 0), ship vy %.1f (expect > 0: bounced), overlap %.3f\n",
               gs.getUniverseState().macro_bodies.size(), gs.getUniverseState().particles.size(), s.m_velocity.y_val, worstOverlap(gs));
    }
    { // 3. plow through a debris cloud: momentum of ship + debris is conserved (no planets = no gravity)
        GameState gs; UIState ui; PhysicsSystem p;
        releaseSpawn(p, gs, ui, true, 60, 8e5, {0, -400}); // 800 particles of mass 1000 each, somewhere in the world...
        V mean(0, 0); for (auto& b : gs.getUniverseState().particles) mean += b.position; mean /= double(gs.getUniverseState().particles.size());
        for (auto& b : gs.getUniverseStateMutable().particles) { b.position += V(0, -400) - mean; b.previousPosition = b.position; } // ...moved to (0,-400)
        Starship& s = gs.getUniverseStateMutable().ship; placeShip(s, {0, 0}, {0, -300});
        auto momentum = [&]() { V m = s.m_velocity * s.mass(); for (auto& b : gs.getUniverseState().particles) m += b.velocity * b.mass; return m; };
        V before = momentum();
        double worst_during = 0;
        for (int t = 0; t < 360; ++t) { p.advanceUniverseState(gs.getUniverseStateMutable(), ShipControls()); worst_during = std::max(worst_during, worstOverlap(gs)); }
        V after = momentum();
        printf("3 debris: momentum before (%.1f, %.1f) after (%.1f, %.1f), ship v (%.1f, %.1f), worst ship overlap seen %.3f px\n",
               before.x_val, before.y_val, after.x_val, after.y_val, s.m_velocity.x_val, s.m_velocity.y_val, worst_during);
    }
    { // 5. why scenario 4 shattered before: impact speed of a ship dropped 7 px above that planet
        double g = 1e9 / (260.0 * 260.0);
        printf("5 drop: surface gravity of a 1e9 planet at the ship's centre ~%.0f px/s^2; falling 7 px gives ~%.0f px/s (shatter at %.0f)\n",
               g, std::sqrt(2 * g * 7), MIN_SHATTER_SPEED);
    }
    { // 4. resting: ship dropped from just above a heavy planet, gravity keeps pulling it in for 4 s
        GameState gs; UIState ui; PhysicsSystem p;
        GravitationalBody pl = planet({0, 0}, 1e9, 200); pl.isForceStatic = true; // shatterable, but the ship can't shatter anything now
        gs.getUniverseStateMutable().macro_bodies.push_back(pl);
        Starship& s = gs.getUniverseStateMutable().ship; placeShip(s, {0, -260}, {0, 0}); // tail fin ~7 px above the surface
        double worst = 0, last_speed = 0;
        for (int t = 0; t < 480; ++t) { p.advanceUniverseState(gs.getUniverseStateMutable(), ShipControls()); worst = std::max(worst, worstOverlap(gs)); last_speed = s.m_velocity.magnitude(); }
        printf("4 resting: planets left %zu, worst overlap %.3f px, ship speed after 4 s %.2f px/s, ship centre y %.2f\n",
               gs.getUniverseState().macro_bodies.size(), worst, last_speed, s.center().y_val);
    }
}
