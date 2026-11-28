// File: Transfer/src/Core/DeterminismCheck.cpp
#include "Core/DeterminismCheck.hpp"

// Custom Imports
#include "Core/GameState.hpp"
#include "Core/UIState.hpp"
#include "DynamoEngine/Math/DeterministicRandom.hpp"
#include "Systems/Physics/PhysicsSystem.hpp"

// Standard Library Imports
#include <cinttypes>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
constexpr int TICKS_PER_SCENE = 1200;   // 10 game seconds
constexpr int TICKS_PER_CHECKSUM = 300; // a checksum line every 2.5 game seconds

// Mixes one double's exact bits into a checksum
void mixDouble(uint64_t& checksum, double value)
{
    checksum = DynamoEngine::combineHash(checksum, DynamoEngine::doubleBits(value));
}

// Every body's position, velocity, mass and radius, and the ship's state, in order
uint64_t universeChecksum(const UniverseState& universe)
{
    uint64_t checksum = 0;
    for (const std::vector<GravitationalBody>* list : {&universe.macro_bodies, &universe.particles})
    {
        checksum = DynamoEngine::combineHash(checksum, list->size());
        for (const GravitationalBody& body : *list)
        {
            mixDouble(checksum, body.position.x_val);
            mixDouble(checksum, body.position.y_val);
            mixDouble(checksum, body.velocity.x_val);
            mixDouble(checksum, body.velocity.y_val);
            mixDouble(checksum, body.mass);
            mixDouble(checksum, body.radius);
        }
    }
    if (universe.has_ship)
    {
        mixDouble(checksum, universe.ship.center().x_val);
        mixDouble(checksum, universe.ship.center().y_val);
        mixDouble(checksum, universe.ship.velocity().x_val);
        mixDouble(checksum, universe.ship.velocity().y_val);
        mixDouble(checksum, universe.ship.throttle());
    }
    return checksum;
}

// Spawns one planet (is_cluster false) or one debris cluster through the real spawn path
void spawn(PhysicsSystem& physics, UniverseState& universe, const SpawnSettings& settings, bool is_cluster,
           DynamoEngine::Vector2D position, DynamoEngine::Vector2D velocity)
{
    PlayerInput input;
    input.commands.spawn = true;
    input.commands.spawn_aim.is_cluster = is_cluster;
    input.commands.spawn_aim.world_position = position;
    input.commands.spawn_aim.with_velocity = (velocity.x_val != 0.0 || velocity.y_val != 0.0);
    input.commands.spawn_aim.world_velocity = velocity;
    physics.updateGravBodyInstantiations(universe, input, settings);
}

// Spawn spots for the clusters: hashed from the scene number, so they're the same everywhere (no std:: random)
DynamoEngine::Vector2D spawnSpot(int scene, int cluster)
{
    uint64_t key = DynamoEngine::combineHash(DynamoEngine::combineHash(2026, static_cast<uint64_t>(scene)),
                                             static_cast<uint64_t>(cluster));
    return {DynamoEngine::hashToRange(DynamoEngine::combineHash(key, 1), -800.0, 800.0),
            DynamoEngine::hashToRange(DynamoEngine::combineHash(key, 2), -800.0, 800.0)};
}

// Sets up scene `scene` (planets now; clusters arrive during the run, see clustersFor)
void setUpScene(int scene, PhysicsSystem& physics, UniverseState& universe, SpawnSettings& settings)
{
    universe.has_ship = (scene == 5);
    universe.has_collision_box = (scene != 4 && scene != 5);
    universe.collision_box_min_corner = {-1000.0, -1000.0};
    universe.collision_box_max_corner = {1000.0, 1000.0};
    settings.is_shatterable = false;
    settings.is_force_static = false;

    switch (scene)
    {
    case 1: // bounce debris piling onto a static planet (dense contacts)
        settings.radius = 150.0;
        settings.mass = 2e8;
        settings.is_force_static = true;
        spawn(physics, universe, settings, false, {0.0, 0.0}, {0.0, 0.0});
        break;
    case 2: // accretable debris and three moving planets (merging, planet-debris contacts)
    case 3: // the same with one static planet
        settings.radius = 60.0;
        settings.mass = 5e7;
        settings.is_force_static = (scene == 3);
        spawn(physics, universe, settings, false, {-300.0, 200.0}, {0.0, 0.0});
        settings.is_force_static = false;
        spawn(physics, universe, settings, false, {400.0, -300.0}, {-120.0, 60.0});
        spawn(physics, universe, settings, false, {100.0, 600.0}, {80.0, -150.0});
        break;
    case 4: // shatterable planets smashing into each other (random fragment speeds and sizes)
        settings.radius = 80.0;
        settings.mass = 5e7;
        settings.is_shatterable = true;
        spawn(physics, universe, settings, false, {-400.0, 0.0}, {900.0, 30.0});
        spawn(physics, universe, settings, false, {400.0, 20.0}, {-900.0, 0.0});
        settings.radius = 50.0;
        settings.mass = 2e7;
        spawn(physics, universe, settings, false, {0.0, 500.0}, {0.0, -1200.0});
        break;
    case 5: // the ship turning and thrusting through debris near a planet (ship contacts, sin/cos of its rotation)
        settings.radius = 120.0;
        settings.mass = 1e8;
        spawn(physics, universe, settings, false, {0.0, -600.0}, {0.0, 0.0});
        universe.ship.moveBy({-500.0, 300.0});
        break;
    default: // 6: negative mass: a negative planet and negative debris among positive ones
        settings.radius = 120.0;
        settings.mass = 1e8;
        spawn(physics, universe, settings, false, {0.0, 0.0}, {0.0, 0.0});
        settings.radius = 40.0;
        settings.mass = -2e6;
        spawn(physics, universe, settings, false, {300.0, 0.0}, {0.0, 150.0});
        break;
    }
    settings.is_force_static = false;
    settings.is_shatterable = false;
}

// The debris clusters scene `scene` spawns: one every 10 ticks, `count` of them
void clustersFor(int scene, SpawnSettings& settings, int& count_out)
{
    settings.radius = 96.0;
    settings.mass = (scene == 6) ? -1e6 : 1e6;
    settings.setFlag(SpawnFlag::Bounce, scene == 1 || scene == 5);
    count_out = (scene == 4) ? 1 : 5;
}

void printLine(FILE* file, const std::string& line)
{
    std::printf("%s\n", line.c_str());
    if (file != nullptr)
    {
        std::fprintf(file, "%s\n", line.c_str());
    }
}
} // namespace

int runDeterminismCheck()
{
    FILE* file = std::fopen("determinism_check.txt", "w");
    printLine(file, "Transfer determinism check: these lines must be identical on every machine and build");
    uint64_t overall = 0;
    for (int scene = 1; scene <= 6; ++scene)
    {
        GameState game_state;
        UIState ui_state;
        PhysicsSystem physics;
        UniverseState& universe = game_state.getUniverseStateMutable();
        SpawnSettings& settings = ui_state.spawnSettings();
        setUpScene(scene, physics, universe, settings);
        int cluster_count = 0;
        clustersFor(scene, settings, cluster_count);

        for (int tick = 1; tick <= TICKS_PER_SCENE; ++tick)
        {
            if (tick % 10 == 0 && tick / 10 <= cluster_count)
            {
                spawn(physics, universe, settings, true, spawnSpot(scene, tick / 10), {0.0, 0.0});
            }
            // The ship (scene 5): full throttle for the first second, then turning left and right in 2 s turns
            ShipControls controls;
            if (scene == 5)
            {
                controls.thrust = (tick <= 120) ? 1 : 0;
                controls.rotation = ((tick / 240) % 2 == 0) ? 1 : -1;
            }
            physics.advanceUniverseState(universe, controls);

            if (tick % TICKS_PER_CHECKSUM == 0)
            {
                uint64_t checksum = universeChecksum(universe);
                overall = DynamoEngine::combineHash(overall, checksum);
                char line[160];
                std::snprintf(line, sizeof(line), "scene %d tick %4d: %3zu planets %5zu particles  %016" PRIX64, scene,
                              tick, universe.macro_bodies.size(), universe.particles.size(), checksum);
                printLine(file, line);
            }
        }
    }
    char line[64];
    std::snprintf(line, sizeof(line), "OVERALL %016" PRIX64, overall);
    printLine(file, line);
    if (file != nullptr)
    {
        std::fclose(file);
    }
    return 0;
}
