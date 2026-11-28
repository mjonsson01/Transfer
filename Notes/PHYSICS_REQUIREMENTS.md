# Physics requirements (Transfer)

Started 2026-11-27 from Marco's answers. This is what every physics change gets graded against. The numbers are
starting guesses, to be tuned by looking at the game; what matters is that every solver is measured the same way.

## The game, in Marco's words

A physics puzzle game with a "highly overcompetent" physics engine: collisions should be extremely pretty and dense.
If the particle count gets high enough, it should also be a passable fluid approximation. The main planets dictate
most of the gameplay; debris is an obstacle and a level object to design around.

## Decisions (2026-11-27)

1. **Debris is gameplay.** It damages the ship, feels gravity fields, and is part of puzzle solutions (negative mass
   and the merge "rail gun" boost). So debris must be exactly as repeatable as planets.
2. **Physics stays on the CPU, identical on Mac and Windows.** GPU math generally differs between vendors and drivers,
   so the GPU is for rendering, not for gameplay physics.
3. **Levels have a fixed amount of rubble.** The designer chooses it; the sandbox gets a setting outside the levels.
   Particle counts never adapt to the player's machine, because rubble volume changes the outcome.
4. **A slow machine plays in slow motion, never differently.** Same ticks, same results, just slower. Thread count
   and render settings may follow the machine; physics content may not.
5. **Target: a 4-core laptop handles at least 50k particles.**
6. **Materials must be extensible:** start with one or two, and keep adding over time (e.g. hot magma from collisions
   that cools and gets more viscous; gas clouds and gas planets later).

## Requirements

### A. Repeatable (because it's a puzzle game)
- A1. Same setup gives a bit-identical outcome: every run, any thread count (1, 2, 4, 8), Mac and Windows.
- A2. A settled scene stays settled: no creeping, no spontaneous spin while the player plans.
- A3. Never explodes: no NaN, no energy from nowhere, no overlap correction flinging particles.
- A4. (Later) Planning-mode prediction matches what happens. Today it predicts planets and the ship only, with no
  collisions; with debris as gameplay this needs a design decision.

### B. Overcompetent, pretty, dense
- B1. Dense clumps hold their shape: overlap too small to see.
- B2. Piles settle within a few seconds instead of boiling forever.
- B3. Impacts look right for the material (rubble bounces, sand doesn't).
- B4. Fast debris doesn't pass through things at normal game speeds.

### C. Performance
- C1. 50k particles at 60 FPS on a 4-core laptop at 1x speed. At 1x a frame runs 2 ticks, so the physics budget is
  about 2 ms per tick, worst case (a deep pile), with 4 threads.
- C2. Thread count is one number; more cores only make it faster, never different.
- C3. The screen stays smooth when physics can't keep up, and the sim speed slider keeps working. DONE 2026-11-27:
  physics runs on its own thread (Core/PhysicsThread.hpp), the main thread renders the newest published copy and
  never waits for ticks. When ticks can't keep up, owed game time beyond MAX_PHYSICS_BACKLOG_TICKS is dropped: the
  simulation runs slower than the slider ("Current Sim Speed" next to the FPS) while the frame rate holds. Results don't
  change (same ticks; checked bit-identical). Player changes are applied between ticks as commands.

### D. Materials (design for growth)
- D1. A material table: adding a material means adding an entry (and only sometimes new behavior).
- D2. Per-particle state that can change: temperature now, phase later (molten / solid).
- D3. Contact behavior comes from the two particles' materials: restitution, friction, cohesion (stickiness), and
  their temperature-dependent versions (magma: runny when hot, stiff when cool).
- D4. Gas is a different kind of interaction (pressure, no packing) and comes later, on the same grid.

## Starting parameters (v0, 2026-11-27): build around these, then tune by feel

Marco can't pick numbers before feeling the gameplay, so these are deliberate guesses to build and measure against.
Change them freely after playing; the scorecard thresholds follow.

**Scale**
| | v0 | today |
|---|---|---|
| Debris in a typical level | 5k - 15k | |
| Debris in a big set piece | 30k | |
| Guaranteed on the reference 4-core laptop, at 1x | 50k | |
| Sandbox maximum | 100k | MAX_LIVE_PARTICLES 20k |
| Debris radius | 2 - 8 px | a radius-96 cluster gives ~2.5 px fragments |
| Sim speed slider | 0.1x - 3x (requested; see C3) | 0.1x - 3x |

**Feel**
| | v0 |
|---|---|
| Settle time after the last impact | <= 3 s |
| Overlap in a settled pile | mean <= 2% of touching distance, worst <= 10% (at 2.5 px radius, 10% is 0.5 px: invisible) |
| Touching neighbors per particle in a pile | <= 6 |
| Debris passing through the ship or a planet | never, below 10,000 px/s |
| Debris passing through debris | allowed above ~1,000 px/s relative speed (invisible inside a cloud) |
| Strongest pull a pile must hold against | 5,000 px/s^2 (today's deep wall scene) |
| Spin of a settled pile | < 1e-4 deg/s |

**Materials** (two to start; the rest show the table can grow)
| Material | Restitution | Friction | Cohesion | Notes |
|---|---|---|---|---|
| Rubble (default debris) | 0.3 | 0.6 | 0 | today's debris is effectively restitution 1.0 (ELASTIC_LOSS_FACTOR 0.994): a big feel change, tune first |
| Bouncy (the Bounce flag) | 0.9 | 0.2 | 0 | |
| Sand (later) | 0.05 | 0.8 | 0 | |
| Ice (later) | 0.5 | 0.05 | 0 | |
| Magma (later) | 0 | hot 0.05 -> cool 0.9 | hot low -> cool high (fuses) | made by violent impacts, cools over time; parameters follow temperature |

**Performance**
| | v0 |
|---|---|
| Debris collisions, 50k deep pile, 4 threads | <= 2 ms/tick |
| Whole tick, same scene | <= 3 ms |
| Reference machine | 4 threads on Marco's Mac until a real 4-core laptop is chosen |
| Per particle, per tick | 4 cores x 2 ms / 50k = 160 ns of CPU time. Today single-threaded: 282 ns in a crushed pile (14.4k: 4.06 ms), 88 ns spread out (80k: 7 ms) |

## Determinism rules (engineering side of A1)

What "identical on Mac and Windows" requires of the code. All four were put in place on 2026-11-27.
- **Fixed tick, fixed order.** Every loop over bodies runs in an order that doesn't depend on threads, timing or
  memory addresses.
- **No unseeded randomness in physics.** Physics randomness is HASHED from `UniverseState::random_seed`, the tick
  (`tick_count`) and the body (DynamoEngine/Math/DeterministicRandom.hpp, `bodyRandomKey` in PhysicsSystem.cpp), never
  drawn from a generator. `randomDouble` (time-seeded, used by shatters and cluster spawns) is gone.
- **No standard-library transcendental functions in physics.** `sqrt`, `floor`, `frexp`, `ldexp`, `remainder` are exact
  everywhere; `sin`, `cos`, `pow`, `cbrt` are not, so physics uses `deterministicSin/Cos/Cbrt`
  (DynamoEngine/Math/DeterministicMath.hpp): fragment layout, ship rotation, stacked-pair directions, accretion radius.
- **Same floating-point code generation.** CMakeLists.txt sets `-ffp-contract=off` (clang/gcc) and `/fp:precise`
  (MSVC) for every target. Never `-ffast-math`. Harness builds use the same flag. Measured: without it every checksum
  of the check changes on the Mac; with it a -O2 build and the Release build match exactly.
- **The check:** `TransferGame --determinism-check` (TransferGame.app/Contents/MacOS/TransferGame on the Mac) plays
  six fixed scenes in ~1 s and writes determinism_check.txt. The Mac's output is saved as
  Notes/determinism_check_mac.txt: the Windows build must print exactly the same lines. (Not yet run on
  Windows as of 2026-11-27.)
- Still to watch: thread count (add 1/2/4/8 to the check when threads arrive) and any new physics code using std::
  math or randomness.

## Scorecard scenes (to build as one harness)

| Scene | Checks | Starting threshold |
|---|---|---|
| Deep wall pile, 14.4k and 50k (1e10 planet outside the box) | B1, B2, C1 | settles within 5 s (mean movement < 0.01 diameters/tick); mean overlap < 2%, worst < 10%; touching per particle <= 6 (the most equal discs can touch without overlapping) |
| The "590" case (4800 debris, planet 590 from the wall) | B2 | settles |
| Two clusters fired together, empty space | A3 | momentum constant (relative 1e-9); kinetic energy never rises |
| Pile on a static planet (spin_harness) | A2 | spin below 1e-4 deg/s |
| Any scene, threads 1/2/4/8, Mac vs Windows | A1 | identical checksums |
| Timing | C1 | <= 2 ms/tick worst case at 50k on 4 threads (reference laptop to be chosen) |

## Chosen direction

**Soft Step** (Erin Catto, Box2D v3; see his 2024 "Solver2D" post) on our sequential solver, then cell coloring for
threads. It keeps the one-pair-at-a-time quality, is deterministic regardless of thread count, runs on the CPU, and
its contacts take material parameters (restitution, friction, cohesion, stiffness) naturally. Stages, each measured
with the scorecard: accumulated clamped impulses -> warm starting (needs particle IDs that survive the sort by cell)
-> soft contacts + relax pass -> restitution step -> friction -> substeps if still needed -> cell coloring + threads.
The determinism fixes above are cheap and can go first.

Rejected (2026-11-27): gather/Jacobi as the main solver (collapses deep piles, see CLAUDE_NOTES); GPU physics
(breaks cross-platform determinism); adapting particle counts to the machine (changes outcomes).

## Open questions
- Rotating, shaped bodies (chunky asteroids, rubble pieces) as a separate body type later. Soft Step was built for
  rigid bodies, so its contacts extend to rotation; disc debris stays as it is.
- Planning-mode prediction with debris (A4).
- Ship damage model from debris impacts.
- Reference "4-core laptop" for C1.
- First materials to ship (rock/rubble only? sand?).
