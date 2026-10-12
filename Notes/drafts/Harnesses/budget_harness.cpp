// Build + run: see README.md in this folder.
// Particle budget: spawns are refused using live + 800 per macro, so physics can never pass MAX_LIVE_PARTICLES.
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "Systems/PhysicsSystem.hpp"
#include <algorithm>
#include <cstdio>

static void spawn(PhysicsSystem& physics, GameState& gs, UIState& ui, bool cluster, double x, double radius)
{
    DEPRECATED_InputState& in = ui.getMutableDEPRECATED_InputState();
    in.selectedRadius = radius;
    in.selectedMass = 1000.0;
    in.mouseCurrPosition = {x, 0.0};
    in.isCreatingParticleCluster = cluster;
    in.isCreatingMacro = !cluster;
    physics.updateGravBodyInstantiations(gs, ui);
}

int main()
{
    {
        GameState gs;
        UIState ui;
        PhysicsSystem p;
        for (int i = 0; i < 30; ++i)
            spawn(p, gs, ui, true, i * 500.0, 50.0);
        printf("30 clusters               -> particles=%zu (expect 20000)\n", gs.getUniverseState().particles.size());
    }
    {
        GameState gs;
        UIState ui;
        PhysicsSystem p;
        for (int i = 0; i < 100; ++i)
            spawn(p, gs, ui, false, i * 500.0, 50.0);
        printf("100 macro spawns, empty   -> macros=%zu (expect 25)\n", gs.getUniverseState().macro_bodies.size());
    }
    {
        GameState gs;
        UIState ui;
        PhysicsSystem p;
        for (int i = 0; i < 12; ++i)
            spawn(p, gs, ui, true, i * 500.0, 50.0); // 9600 particles
        for (int i = 0; i < 100; ++i)
            spawn(p, gs, ui, false, 1e5 + i * 500.0, 50.0);
        printf("9600 particles + macros   -> macros=%zu (expect 13: 9600 + 13*800 = 20000)\n",
               gs.getUniverseState().macro_bodies.size());
    }
    {
        GameState gs;
        UIState ui;
        PhysicsSystem p;
        for (int i = 0; i < 100; ++i)
            spawn(p, gs, ui, false, i * 90.0, 50.0); // neighbours overlap
        auto& macros = gs.getUniverseStateMutable().macro_bodies;
        for (size_t i = 0; i < macros.size(); ++i)
            macros[i].velocity = {(i % 2 == 0) ? 3000.0 : -3000.0, 0.0};
        size_t max_particles = 0, max_bodies = 0;
        for (int t = 0; t < 600; ++t)
        {
            p.advanceUniverseState(gs.getUniverseStateMutable(), ui);
            max_particles = std::max(max_particles, gs.getUniverseState().particles.size());
            max_bodies = std::max(max_bodies, gs.getUniverseState().particles.size() + gs.getUniverseState().macro_bodies.size());
        }
        printf("smash all 25 macros       -> max particles=%zu, max bodies=%zu, macros left=%zu (expect <= 20000)\n",
               max_particles, max_bodies, gs.getUniverseState().macro_bodies.size());
    }
}
