// Build (from the repo root, macOS):
//   clang++ -std=c++20 -O2 -w -INotes/drafts/Harnesses -ITransfer/src -IThirdParty/SDL3/include -IThirdParty/SDL3_ttf/include Notes/drafts/Harnesses/grid_equivalence_harness.cpp Notes/drafts/Harnesses/OldParticleGrid.cpp $(find Transfer/src -name '*.cpp' ! -name main.cpp) -LThirdParty/SDL3/lib -LThirdParty/SDL3_ttf/lib -lSDL3 -lSDL3_ttf -Wl,-rpath,$PWD/ThirdParty/SDL3/lib -Wl,-rpath,$PWD/ThirdParty/SDL3_ttf/lib -o /tmp/grid_eq && /tmp/grid_eq
// The new column-stencil grid must report exactly the same unordered pairs as the old 5-cell grid.
#include "OldParticleGrid.hpp"
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "Systems/PhysicsSystem.hpp"
#include "Utilities/Physics/UniformParticleGrid.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <random>
#include <utility>

template <typename Grid> static std::vector<std::pair<size_t, size_t>> pairs(const std::vector<GravitationalBody>& ps, double& ms)
{
    Grid g; std::vector<size_t> c; std::vector<std::pair<size_t, size_t>> out;
    auto t0 = std::chrono::steady_clock::now();
    g.build(ps);
    for (size_t i = 0; i < ps.size(); ++i) { g.queryCandidates(i, ps, c); for (size_t j : c) out.push_back({std::min(i, j), std::max(i, j)}); }
    ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    std::sort(out.begin(), out.end());
    return out;
}

static bool check(const char* name, const std::vector<GravitationalBody>& ps)
{
    double old_ms, new_ms;
    auto a = pairs<OldParticleGrid>(ps, old_ms), b = pairs<UniformParticleGrid>(ps, new_ms);
    bool dup = std::adjacent_find(b.begin(), b.end()) != b.end();
    bool same = (a == b);
    printf("%-30s n=%6zu pairs old=%7zu new=%7zu  %s%s  old=%.2f ms new=%.2f ms\n", name, ps.size(), a.size(), b.size(),
           same ? "IDENTICAL" : "DIFFERENT", dup ? " (DUPLICATES!)" : "", old_ms, new_ms);
    return same && !dup;
}

int main()
{
    bool ok = true;
    std::mt19937 rng(7);
    for (int scene = 0; scene < 4; ++scene) // random scatter, mixed sizes, around the origin (negative coordinates too)
    {
        std::uniform_real_distribution<double> pos(-300.0 * (scene + 1), 300.0 * (scene + 1));
        std::uniform_real_distribution<double> rad(0.5, scene == 3 ? 12.0 : 3.0);
        std::vector<GravitationalBody> ps(20000);
        for (auto& p : ps) { p.position = {pos(rng), pos(rng)}; p.radius = rad(rng); }
        char name[64]; snprintf(name, sizeof(name), "random scatter %d", scene);
        ok &= check(name, ps);
    }
    { GameState gs; UIState ui; PhysicsSystem p; DEPRECATED_InputState& in = ui.getMutableDEPRECATED_InputState();
      for (int i = 0; i < 25; ++i) { in.selectedRadius = 50; in.selectedMass = 1000; in.mouseCurrPosition = {(i % 5) * 110.0 - 200, (i / 5) * 110.0 - 200}; in.isCreatingParticleCluster = true; p.UpdateGravBodyInstantiations(gs, ui); }
      ok &= check("packed clusters", gs.getParticles());
      for (int t = 0; t < 240; ++t) p.UpdateSystemFrame(gs, ui);
      ok &= check("clusters after 2 s", gs.getParticles()); }
    printf(ok ? "ALL IDENTICAL\n" : "MISMATCH\n");
}
