// Build (from the repo root, macOS; any .cpp under Transfer/src except main.cpp):
//   clang++ -std=c++20 -O1 -w -ITransfer/src -IThirdParty/SDL3/include -IThirdParty/SDL3_ttf/include Notes/drafts/Harnesses/budget_harness.cpp $(find Transfer/src -name '*.cpp' ! -name main.cpp) -LThirdParty/SDL3/lib -LThirdParty/SDL3_ttf/lib -lSDL3 -lSDL3_ttf -Wl,-rpath,$PWD/ThirdParty/SDL3/lib -Wl,-rpath,$PWD/ThirdParty/SDL3_ttf/lib -o /tmp/budget && /tmp/budget
// Particle budget: spawns are refused using live + 800 per macro, so physics can never pass MAX_LIVE_PARTICLES.
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "Systems/PhysicsSystem.hpp"
#include <algorithm>
#include <cstdio>

static void spawn(PhysicsSystem& physics, GameState& gs, UIState& ui, bool cluster, double x, double radius)
{
    DEPRECATED_InputState& in = ui.getMutableDEPRECATED_InputState();
    in.selectedRadius = radius; in.selectedMass = 1000.0; in.mouseCurrPosition = {x, 0.0};
    in.isCreatingParticleCluster = cluster; in.isCreatingMacro = !cluster;
    physics.UpdateGravBodyInstantiations(gs, ui);
}

int main()
{
    { GameState gs; UIState ui; PhysicsSystem p;
      for (int i = 0; i < 30; ++i) spawn(p, gs, ui, true, i * 500.0, 50.0);
      printf("30 clusters               -> particles=%zu (expect 20000)\n", gs.getParticles().size()); }
    { GameState gs; UIState ui; PhysicsSystem p;
      for (int i = 0; i < 100; ++i) spawn(p, gs, ui, false, i * 500.0, 50.0);
      printf("100 macro spawns, empty   -> macros=%zu (expect 25)\n", gs.getMacroBodies().size()); }
    { GameState gs; UIState ui; PhysicsSystem p;
      for (int i = 0; i < 12; ++i) spawn(p, gs, ui, true, i * 500.0, 50.0);  // 9600 particles
      for (int i = 0; i < 100; ++i) spawn(p, gs, ui, false, 1e5 + i * 500.0, 50.0);
      printf("9600 particles + macros   -> macros=%zu (expect 13: 9600 + 13*800 = 20000)\n", gs.getMacroBodies().size()); }
    { GameState gs; UIState ui; PhysicsSystem p;
      for (int i = 0; i < 100; ++i) spawn(p, gs, ui, false, i * 90.0, 50.0); // neighbours overlap
      auto& macros = gs.getMacroBodiesMutable();
      for (size_t i = 0; i < macros.size(); ++i) macros[i].velocity = {(i % 2 == 0) ? 3000.0 : -3000.0, 0.0};
      size_t max_particles = 0, max_bodies = 0;
      for (int t = 0; t < 600; ++t)
      {
          p.UpdateSystemFrame(gs, ui);
          max_particles = std::max(max_particles, gs.getParticles().size());
          max_bodies = std::max(max_bodies, gs.getParticles().size() + gs.getMacroBodies().size());
      }
      printf("smash all 25 macros       -> max particles=%zu, max bodies=%zu, macros left=%zu (expect <= 20000)\n",
             max_particles, max_bodies, gs.getMacroBodies().size()); }
}
