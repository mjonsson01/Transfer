# Transfer: Claude handoff notes

**Fresh Claude session, on any device: read this whole file, then `Notes/REWORK.txt`, before doing anything.**
This folder is the only memory that travels between machines (Claude's own memory is per-device).
Last full rewrite: 2026-09-27, on macOS, branch `EngineSep` at `959b486` ("Finished initial full implementation of visor_view").
Updated 2026-09-29 on Windows, at `7f852af` (everything committed and pushed; build/scripts rework done on both platforms).
Updated 2026-10-09 on macOS: EngineSep merged to master (PR 31), Mac packaging merged (PR 32/33, master at `a1cb3a7`).
The folder was renamed `CLAUDE/` -> `Notes/`. Current working branch: `Bugfixes` (REWORK items first, then new features).

---

## 1. Working with Marco (READ FIRST)

**Rules**
- The codebase is **READ-ONLY**. Write only inside `Notes/` unless Marco *explicitly* says to change repo files ("apply these changes", "finish applying").
  Permission is per request, not standing.
- Never hand over a patch to apply blindly. Marco said: "I won't apply a patch blindly like that. Walk me through it."
  Explain every change: what, why, and how it connects to what he already knows.
- Keep `Notes/REWORK.txt` current. Format: `[BUG]` / `[REFACTOR]` / `[PERFORMANCE UPDATE]`, then function, then file, then about one sentence.
  When Marco reports an item fixed, DELETE it (reword if only partly fixed). Also prune items whose code no longer exists.
- Before presenting code, verify it in a scratch copy of the repo: rsync it (excluding build dirs, Notes and .git), configure with
  `-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=<repo>/build-tests/_deps/googletest-src`, build the game and tests, and run TidyEngine.
  Offer drafts; only write to the repo when told.
  - **Re-copy `Transfer/Assets/Shaders/*.spv` right before any visual test.** They're gitignored build output that only MakeTransfer
    rebuilds; a stale copy drew every text quad as a solid white box (2026-09-28, wrongly blamed on an SDL change at first).
  - Windows: the scratchpad path is too long for MSBuild (MAX_PATH 260) -> `subst Q: <scratchpad>` and build in `Q:\...`, then `subst Q: /D`.
    A scratch copy also needs its own `ThirdParty\` (copy the real one; don't rebuild shadercross for a scratch test).
  - Claude can check visuals on Windows: launch the exe, PowerShell `CopyFromScreen` on its window rect, Read the PNG
    (script pattern: `$p.MainWindowHandle` + `GetWindowRect`; it sometimes grabs the wrong window, so check the title bar).
  - The FontAtlas tests check metrics only, not pixels: passing tests do NOT prove text renders.
- Ask when a design choice is genuinely his (widget behavior, look, placement). Don't guess at it.

**Who he is and how he learns**
- Marco is a self-described novice C++ developer. He built Transfer himself; it's a passion project he intends to **ship on Steam**,
  and he won't settle for incomplete. Judge advice by what gets it shipped.
- It's an educational project: quiz him on his own code, explain the *why*, and grade answers honestly
  (e.g. "mostly right, but ..."). He likes being quizzed and asks for brutal grading.
- **Walkthrough mode** (used for the visor dropdown, and he loved it): design the header together, then HE types the code
  (typing, not pasting, is how he learns it). Give code in parts with "things to notice", review what he typed,
  then give him tests to write from a table of name / steps / expectations. Ask one quiz per step.
- He pushes back, and he's often right. Examples: "why nointerpolation on a uint?", "why is the bool unused?", "closing on
  mouse-exit is a hidden assumption". Concede clearly when he is right.

**Code style he wants**
- Readability first: plain, explicit code over clever code. No variadic templates or perfect forwarding. Explain any new
  syntax piece by piece.
- KISS: avoid `std::optional` where a named sentinel works (`static constexpr int NO_OPTION = -1;`).
- Named aliases for callbacks (`using UIAction = std::function<void()>;`), named setters (`setOnClick`), and
  **create → configure → hand over** (`make_unique`, setters, then `addChild(std::move(ptr))`).
- Comments explain the non-obvious *why*, and a rule shared between two places is written on both sides (e.g. the VisorView ⇄ shader numbers).
- One source of truth: derive values instead of storing copies (option rows are computed from `m_rect`; the dropdown pulls
  its selection from `CameraState::visor_view` every frame).
- His files are NOT clang-formatted wholesale. Don't reformat lines you didn't change; keep diffs minimal.

**Naming (clang-tidy `readability-identifier-naming`, config in `Transfer/src/.clang-tidy`)**
- Types, namespaces and enum values: PascalCase. Methods and functions: camelCase. Getters have no `get` prefix (`rect()`, `value()`).
- Locals, params and public struct fields: snake_case. Private/protected members: `m_snake_case`. Constants: UPPER_SNAKE.
- Old game code still uses older styles (`gameState`, `getCameraStateMutable`). Only new/touched code follows the rules.
- Files: `// File: <path>` header, `#pragma once`, then import groups `// Custom Imports`, `// SDL Imports`, `// Standard Library Imports`.
- `misc-include-cleaner` is on with `MissingIncludes: false` (unused includes flagged; missing ones are not; no IWYU pragmas).

---

## 2. Repo, build, test, tools

- Repo: `github.com/mjonsson01/Transfer`. master holds the finished rework (EngineSep merged 2026-09-29); work happens on topic branches merged by PR (now `Bugfixes`).
  C++20, CMake with GLOB_RECURSE (new files are picked up automatically), SDL3 + SDL_GPU + SDL3_ttf.
- `Transfer/src/` = game + `DynamoEngine/` (STATIC lib target, linked PUBLIC to SDL). **The engine must never include game code**;
  the test target builds the engine alone, which enforces this.
- Engine warnings: `-Wall -Wextra -Wimplicit-fallthrough -Wno-unused-parameter` (MSVC `/W4 /wd4100`), PRIVATE to the engine.
  Why unused-parameter is off: overrides must match base signatures (e.g. `onMouseReleased(pos, released_inside)`).
- **Scripts live in `Scripts/Windows/*.bat` and `Scripts/Apple/*.sh`** (moved 2026-09-29; `Scripts/README.md` is the user guide).
  Same four on both platforms, same arguments; each `cd`s to the repo root itself (`%~dp0..\..` / `$(dirname "$0")/../..`).
  `.gitattributes` forces `*.sh` LF and `*.bat` CRLF. Mac scripts keep their executable bit only because they were `git mv`'d.
  - `SetupDependencies` (once after cloning; safe to re-run, skips what `VERSION.txt` says is installed) fills gitignored `ThirdParty/`:
    - Windows: SDL3 3.4.16 + SDL3_ttf 3.2.2 prebuilt VC zips, SHA-256 pinned (curl, tar, PowerShell Get-FileHash).
    - Mac: SDL3 + SDL3_ttf built from source at pinned git COMMITS (tags release-3.4.16 = fa2c02bb..., release-3.2.2 = a1ce3670...),
      SDL_ttf with `SDLTTF_VENDORED=ON` (FreeType/HarfBuzz are git submodules: the 1.5 MB source tarball does NOT contain them).
    - Both: `ThirdParty/ShaderCross` = SDL_shadercross built from source at commit 1ff05bec... (no releases/tags exist), vendored
      (builds DXC + SPIRV-Cross: slow, git + Python needed), `cmake --install` -> `bin/`, plus SDL3.dll / libSDL3*.dylib copied in.
      Mac passes `CMAKE_INSTALL_RPATH=@executable_path/../lib` (upstream sets none). Windows builds in `ThirdParty\_src\sc`
      (short on purpose: MAX_PATH) with `git -c core.longpaths=true`.
    - Stable SDL only: an ODD SDL3 minor (3.5.x) is a prerelease. Marco's old `C:\SDL3` was a 3.5.0 dev snapshot.
  - `MakeTransfer [debug|release|clean]` lints (SKIP_TIDY=1 skips), **compiles HLSL** with `ThirdParty/ShaderCross/bin/shadercross`
    (mac → .msl, Windows → .spv, into Transfer/Assets/Shaders), then builds `build/`. **Only this recompiles shaders.**
    The .bat compiles in a for loop and `goto :shaderFailed` on the first error: `exit /b 1` INSIDE a for loop lost its exit code
    (returned 0 under `cmd /c`; verified 2026-09-29). Plain `if` blocks are fine. Same pattern for SetupDependencies' tool check.
  - `RunTests [ctest args]`: googletest v1.17.0 via FetchContent, option `TRANSFER_BUILD_TESTS`, own dir `build-tests/`.
    The .bat passes `--config`/`-C Debug` for multi-config (Visual Studio) generators.
  - `TidyEngine`: clang-tidy over every engine .cpp/.hpp (34 files, all clean). SDL headers from `ThirdParty/`. The mac version needs
    `-isysroot $(xcrun --show-sdk-path)` and `-resource-dir $(clang -print-resource-dir)`. The .bat uses VS Code's bundled clang-tidy.
    You can pass a specific binary: `CLANG_TIDY=...`.
  - Verified 2026-09-29 on Windows, run from `C:\`: Tidy clean, 120/120, MakeTransfer compiles all 10 shaders + builds.
    **Windows shadercross source build VERIFIED 2026-09-29** (Marco ran SetupDependencies.bat; installs ThirdParty\ShaderCross at
    1ff05bec...). It needed two lessons: (1) MAX_PATH: DXC's MSBuild tlog paths passed 260 chars under the repo path even from
    `_src\sc`, so the .bat now `subst`s a free letter T:-Z: onto `ThirdParty\_src` for the build (unmapped on every exit path, sources
    deleted through the short drive). (2) A leftover tree from a failed run broke the next one: `rmdir /s /q` silently left files
    (~27 idle MSBuild node-reuse processes held them), then `git submodule update` refused non-empty folders and the script still
    reached configure (errors about CheckAtomic/GetHostTriple). Fix for the user: `taskkill /im MSBuild.exe /f`, delete `_src`.
    Script hardening NOT yet done: fail if rmdir leaves anything; check external\DirectXShaderCompiler\CMakeLists.txt after submodules.
    Mac: first run (2026-09-29) froze the Mac: `cmake --build --parallel` with no number = UNLIMITED `make -j` with the Makefiles
    generator -> ~50 clang processes, memory 100% (MSBuild's /m is CPU-bounded, so Windows never showed it). FIXED 2026-09-29:
    `--parallel "$BUILD_JOBS"` = min(CPU cores, RAM GB / 2), at least 1 (sysctl hw.ncpu / hw.memsize).
  - **Mac VALIDATED 2026-09-29 (reported by Marco; Claude saw no Mac output):** SetupDependencies.sh (with the job limit; it printed
    12 jobs), RunTests.sh and MakeTransfer.sh all work. Activity Monitor showed MORE than 12 `clang` processes during the build: likely
    driver + `clang -cc1` pairs or link steps (make's -j12 limits jobs, not processes); never confirmed with `ps`. The number that
    matters is Memory Pressure staying out of red.
- Tests: top-level `Tests/` mirrors `src/` (Tests/DynamoEngine/{Input,Rendering,UI,UI/Widgets}); include roots are `Transfer/src` and `Tests`.
  `Tests/TestingUtilities/VectorsNear.hpp` is an AssertionResult helper. `EXPECT_DEBUG_DEATH` is used for assert paths. The CMake define
  `TRANSFER_TEST_FONT_PATH` points FontAtlas tests at `Transfer/Assets/Fonts/SpaceMono-Regular.ttf`. **120/120 pass at handoff.**
  Known nit: Test_SceneManager.cpp lives in Tests/DynamoEngine/UI/ but its header comment says .../Scenes/.
- `.vscode/c_cpp_properties.json` is GENERATED by CMake on configure (option TRANSFER_GENERATE_VSCODE_CONFIG); hand edits are lost.
- Docs: `PythonScripts/create_docs.py` regenerates `Documentation/class_map.md` + `structure.txt` (stale since the UI rework; rerun it).
- `.gitignore` no longer ignores `*.md` (checked 2026-09-29), so Notes/ and Scripts/README.md are tracked normally.
- SDL in CMake (both platforms since 2026-09-29): `find_package(SDL3/SDL3_ttf CONFIG PATHS ${TRANSFER_THIRD_PARTY_DIR}/... NO_DEFAULT_PATH)`
  (cache var, default `<repo>/ThirdParty`) -> imported targets `SDL3::SDL3` / `SDL3_ttf::SDL3_ttf`. No Homebrew anywhere any more.
  Old dependency locations are gone on Windows (checked 2026-09-29): `C:\SDL3`, `..\..\SDL3_TTF`, `LocalShaderCross/`, `LocalSpirvCross/`
  (their .gitignore lines too). On the Mac, Homebrew's sdl3/sdl3_ttf are no longer used (uninstalling is optional).
  - Windows details:
    DLLs are copied next to the game AND TransferTests with `$<TARGET_RUNTIME_DLLS:...>` + `COMMAND_EXPAND_LISTS` (needs CMake 3.21).
    The test copy must be added BEFORE `gtest_discover_tests`: discovery runs the exe post-build, and a missing DLL = exit 0xc0000135.
    `SDL3_INCLUDE_DIRS` / `SDL3_DLLS` no longer exist.
  - The old `CLAUDE/drafts/` folder (stress harness, sub-pixel fade diff) was never committed and is GONE (not found on the Mac
    2026-10-09). The descriptions in section 5 are all that's left of them.

---

## 3. Architecture (current, after the UI/scene rework)

**Game loop** (`Core/Game.cpp`, `Game::Run`): frame_seconds → updateFPS → `ProcessInput(frame_seconds)` → UpdateInstantiations
(spawns from input) → PlayAudio → fixed-step physics (120 Hz accumulator × time scale, only while `scenes.currentScene().settings().runs_simulation`)
→ alpha → RenderFrame → **`scenes.applyPendingSwitch()`** (deferred scene switches happen here) → limitFrameRate (60 FPS).
`Game` owns GameState, UIState, InputSystem, PhysicsSystem, RenderSystem, AudioSystem, `DynamoEngine::SceneManager scenes`.

**Input** (`Systems/InputSystem`): SDL events → `DynamoEngine::SDLInputIntake` (free fn `translateSDLEvent`) → `DynamoEngine::InputState`
(edges vs levels, `beginInputFrame`, focus-lost clears held keys) → **UI first**: `updateSceneUI` does setWindowSize, `updateElements(dt)` and
`processInput` (the order matters), giving a `UIInputResult`. Then the game side: `UIInputConsumed = pointer_captured`,
`translateGameInputs(game_state, ui_state, scenes)` / `translateMenuInputs`, writing into **`Core/DEPRECATED_InputState`**. That struct is still
the bridge to physics and rendering (spawn flags, slider values, preview ghost) and goes away with the action map (E).
Key bindings: Esc = pause/resume (requestSwitch), Backspace/Delete = clear, WASD = starship, left/right release = spawn macro/cluster,
Shift = initial velocity, middle mouse = pan, wheel = zoom, **Tab = cycle visor view** (`updateVisor` in the sim branch).

**Scenes**: engine `SceneId = uint32_t`, `Scene` (SceneSettings{runs_simulation, draws_world} with presets menu()/simulation(), owns a
`UIRoot`, virtual onEnter/onExit/update), `SceneManager` (addScene with duplicate check, requestSwitch with unknown-id check: log + assert,
ignored in Release; deferred applyPendingSwitch, which also calls `ui().cancelPointerInput()` on the scene being left).
Game: `Scenes/TransferScenes.hpp/.cpp`: `namespace TransferScene { enum : DynamoEngine::SceneId { StartMenu, Game, Pause, TestVisual }; }`
(a plain enum on purpose, so it needs no casts), `addTransferScenes(scenes, ui_state, game_state)`, and one `buildXScene` per scene.
Lambdas capture Game members by reference, which is safe because they outlive the scenes.

**Engine UI API** (`DynamoEngine/UI/`): this is what the game builds with.
- `UIElement`: owns children via unique_ptr (`addChild` returns `UIElement&`); `setPlacement(UIPlacement{align (9 spots + Fill), size, margin})`;
  `updateLayout`, `rect()`, virtual `containsPoint` (left/top edges inside); `setVisible`; `setLayer(HUD < Menu < Overlay)` (inherited);
  `setZIndex` (siblings only); `requestSound(UISound)` bubbles up to the root. Overridables: `update(dt)`, `draw(builder) const`
  (self only), `onMousePressed` (return true = "this press is mine"), `onMouseDragged`, `onMouseReleased(pos, released_inside)`,
  `onMouseEntered`, **`onMouseHover(pos)`** (every frame while hovered), `onMouseExited`.
- `UIRoot`: UI scale = (window height / `UI_REFERENCE_HEIGHT` 720) × player scale. It works in UI space: `uiSpaceSize()` and `screenToUISpace()`.
  It rebuilds layout and draw order every frame. Draw order = tree order (parent first, siblings by z), then stable-sorted by layer. Hit test =
  reverse draw order. `processInput` handles hover, offers a press top-down until one returns true, then capture/drag/release (left button only).
  `cancelPointerInput`, `setSoundHandler`, `topmostElementAt`, `elementsInDrawOrder`.
- Widgets (`UI/Widgets/`):
  - `UILabel` (setText / setTextSource pull).
  - `UIButton` (setOnClick; clicks on release inside; gray 128 / hovered 108 / held 88).
  - `UISlider` (label, `SliderMapping{to_value, to_position}` with `linear()`; ticks by knob position, 30 per track; `setValue` doesn't call back;
    the knob highlights only while hovered or dragged).
  - `UIRow` / `UIColumn` (size to their children).
  - **`UIDropdown`** (Marco typed it). Constructed as `(title, option_labels)`.
    - The button reads "Title: Choice"; it **selects on press**.
    - `setOnOptionChosen(int)` is called only when the choice CHANGES. Re-choosing plays a click but reports nothing.
    - `setSelectionSource(int())` pulls the selection every frame; `selectOption` doesn't call back.
    - While open, `containsPoint` is true everywhere: a press outside closes the list and is still captured, so no spawn happens behind it.
    - It sets the Overlay layer in its constructor. Option rows are computed from `m_rect`, never stored.
    - It does NOT close on mouse-exit (Marco rejected that as a hidden assumption).
- Rendering support: `Rendering/FontAtlas` (`buildAtlas(font, font_size, pixel_scale)` bakes at the real pixel size; metrics are in UI points;
  1024² atlas), `UI/UIGeometryBuilder` (addRect, addText with pixel-snapped start, addTextCentered), `Rendering/UIVertex` (modes None/Solid/Textured).
  RenderSystem re-bakes the atlas when `uiScale × SDL_GetWindowPixelDensity` changes and pushes `uiSpaceSize` to the UI shader.

**Game scenes' UI**
- StartMenu: "Play Game" button.
- Game:
  - FPS label (top-left).
  - Slider row (bottom-center): speed (linear), radius (linear), and mass (game-side `massMapping()`, a signed log curve).
    The sliders write into DEPRECATED_InputState.
  - **Visor dropdown, top-right**, wired to `camera_state.visor_view`.
- Pause: "Resume" button.
- TestVisual: runs the simulation and draws nothing.

**Visor / view modes**
- `Core/CameraState.hpp`: `enum class VisorView : uint32_t { Realistic=0, Mass=1, Charge=2, Temperature=3 }`. The numbers must match the
  shader's `VIEW_*` constants. There is also `inline VisorView nextVisorView(VisorView)` (a switch) and the field `CameraState::visor_view`.
- `buildCameraConstants` sends `static_cast<uint32_t>(visor_view)` as `viewMode`.
- In `UnifiedGravBody.vert.hlsl`, the value passes through `VertexOutput` (the last field, TEXCOORD14) to `.frag.hlsl`.
- The frag branches: Mass keeps the real coloring. Realistic, Charge and Temperature are PLACEHOLDERS (green / red / blue for now), and
  anything else is debug magenta.
- The start view is currently `Realistic` (green placeholder); Marco may want `Mass` at startup.

**Physics** (unchanged by the rework): `GameState` holds `macroBodies` and `particles` (`std::vector<GravitationalBody>`).
- Integration: velocity Verlet KDK (half kick + drift → forces → half kick), with collisions handled before integration.
- Gravity: macro–macro and particle–macro only, with Plummer softening.
- Particle broadphase: `UniformParticleGrid`.
- Shattering: Vogel-spiral fragments, 800 fragments.

---

## 4. Design intent confirmed by Marco (don't argue these away)
- Accretion merges keep the heavier body's position; only velocity is updated (momentum-conserving).
- `DEFAULT_FRAGMENT_COUNT = 800` is fixed on purpose. Culling particles with radius < 1 is intentional.
- Gentle accretion should crumble macro bodies into fragments rather than pop.
- Pause is a full separate scene (no scene stack for now).
- Mouse is the primary input; keyboard focus (G) is planned. NO controller support planned. Text input will be needed eventually.
- World-space UI wanted later (H): speed readouts on velocity vectors (constant on-screen size), hover readouts for macro bodies,
  predictive dotted trajectories.
- Visor styling: "match the buttons now, style later".

---

## 5. Status and roadmap

**Done, all merged to master:**
- A: input.
- Engine split: `DynamoEngine` static lib, googletest, TidyEngine.
- B1–B5: engine UI (geometry, UIElement, UIRoot, input dispatch, widgets).
- S: SceneManager.
- B5b/B5c: game switched to scenes + widgets; old UISystem, old scenes and Entities/UIElements deleted.
- Crisp scaled text.
- Knob-only slider hover.
- **D: visor dropdown + Tab cycling + shader view modes.**
- **Build rework (2026-09-28/29, verified on Windows AND Mac):** all scripts in `Scripts/{Windows,Apple}` (+ `Scripts/README.md`);
  `SetupDependencies` pins and installs SDL3 3.4.16, SDL3_ttf 3.2.2 and shadercross into gitignored `ThirdParty/`; CMake finds SDL
  only there on both platforms; DLLs copied next to the game and tests; `.gitattributes` for script line endings.
- **Mac packaging (2026-09-29 -> 2026-10-09, PR 32/33; confirmed working by Marco):** Release `TransferGame.app` gets
  `INSTALL_RPATH @executable_path/../Frameworks` + `BUILD_WITH_INSTALL_RPATH`, SDL dylibs copied into Contents/Frameworks, then
  `codesign --force --deep --sign -` as the LAST post-build step. Lesson: copy each dylib TO its install name
  (`IMPORTED_SONAME_RELEASE`, e.g. libSDL3_ttf.0.dylib), not the real file's name (`IMPORTED_LOCATION_RELEASE`): the first
  version kept the real name and dyld failed "Library not loaded: @rpath/...SDL3_ttf...".
  Checked by Claude 2026-10-09 with otool/codesign: the exe's ONLY LC_RPATH is @executable_path/../Frameworks (no ThirdParty path),
  both dylibs load via @rpath, `codesign -v --deep` passes. Still arm64-only with `minos 26.0` (see REWORK BUILD / TOOLING).
  Sharing: Mac zip (`ditto -c -k --keepParent`) runs after the user allows it once (ad-hoc signed = Gatekeeper warning);
  Windows zip = build\Release minus DynamoEngine.lib (exe, SDL3.dll, SDL3_ttf.dll, Resources\), SmartScreen warns once.
- **Bug batch 1 (2026-10-09, branch `Bugfixes`, typed by Marco in walkthrough mode, tested in-game):**
  (1) `Game::Run` clamps physics_time_accumulator to `MAX_PHYSICS_STEPS_PER_FRAME` (8, PhysicsConstants.hpp) x PHYSICS_TIME_STEP
  right after the time-scaled `+=` (stops the spiral of death; excess time is dropped). (2) `getCollisionInfo`: should_blow_up uses
  approaching_speed = -normal_speed (normal points a -> b), so separating overlaps no longer shatter; abs_normal_speed is still used
  by the MAX_ACCRETION_COLLISION_SPEED check. (3) All three collision passes `break` out of the inner loop once the OUTER body is
  isMarkedForDeletion (no double absorb / double shatter). Marco also renamed gameState -> game_state, originalBody -> original_body
  etc. in PhysicsSystem.hpp/.cpp in the same change.
- **Particle ghosting fix (2026-10-09, `Bugfixes`, typed by Marco, tested in-game):** particle pairs with mass ratio >= 8 phased
  through each other (gentle-contact branch -> handleAccretion -> `if (heavier.isParticle) return;` = no response at all). Fix: in
  handleDynamicCollision, `else if (heavier.isParticle) handleElasticCollisions(lighter, heavier);` BETWEEN can_crumble and
  handleAccretion (so a shatterable macro lighter than a particle still crumbles). handleAccretion's early return stays as the guard.
  Verified with a harness (mass 1 vs 100, head-on): before = swapped sides untouched; after = bounce matching textbook elastic x 0.994.
  Marco's wording: particles can't ACCRETE (absorb); they can BE accreted. He was right when Claude claimed otherwise.
- **Particle budget + render overflow safety (2026-10-09, `Bugfixes`, typed by Marco, tested in-game: hit the max, spawns refused,
  no warnings).** Design is Marco's (he rejected 'bounce at the limit' and 'cap macros' as wrong):
  - `potentialParticleCount` = live particles (incl. pending) + macros x DEFAULT_FRAGMENT_COUNT (FLAT 800, not survivableFragmentCount:
    accretion grows a macro's exact count after budgeting; it can't grow the flat one). `hasRoomForParticles(gs, n)` = potential + n
    <= MAX_LIVE_PARTICLES (20000). ONLY player spawns check it: createMacroBody asks for 800, createParticleCluster for its exact
    survivableFragmentCount (all or nothing). Shatter/crumble have NO checks (old elastic 'Hard safety net' deleted) = physics runs in
    full and still can't pass the limit. Empty scene -> 25 planets max (Marco OK'd the overhead; raise MAX_LIVE_PARTICLES if tight).
  - The guarantee RELIES ON: particles never shatter (isShatterable false), particle-particle accretion off, promotion off. Changing
    any of those breaks the bound -> RenderSystem prints `WARNING: ... body buffers ... Is the particle budget still being enforced?`
  - MAX_MACRO_BODIES removed (implied by the budget). `INITIAL_UNIFIED_BODY_CAPACITY = MAX_LIVE_PARTICLES + 1` (+1 = spawn preview;
    budget caps ALL bodies at MAX because each macro costs 800 but is 1 body). Buffers reserved at full size at startup via
    `createUnifiedBodyBuffers(capacity)`; uploadUnifiedBodies grows them (doubling from max(capacity, INITIAL)) ONLY in an emergency,
    with the warning. renderBodies draws unifiedBodyVertices.size() (no recount).
  - UI / starship / velocity-vector uploads: if over their MAX, print + resize() the vector (their draws use .size()). They don't grow:
    their sizes only change with code, so overflow = a bug.
  - Harnesses (real PhysicsSystem, no GPU) in `Notes/drafts/Harnesses/` (ghost_harness, budget_harness; build line in each header).
    GPU growth was checked with a throwaway hacked copy (tiny initial capacity + 30000 injected particles).
  - Open UX question (Marco's call): a refused spawn is silent; feedback (sound / flash / 'limit reached' text) later. Level design
    will also set per-level body limits (Marco, 2026-10-09).
- **GPU cycling (2026-10-09, `Bugfixes`, typed by Marco, checked in-game):** every PER-FRAME upload (bodies, starship, UI, stars,
  velocity vectors) now passes `cycle = true` to both SDL_MapGPUTransferBuffer and SDL_UploadToGPUBuffer, so the CPU never rewrites
  memory a still-running earlier frame references. rebuildFontAtlas stays `false` (brand-new transfer buffer + texture each time =
  never bound). Rule: cycle when rewriting a resource an earlier, possibly unfinished command referenced.
- **First profile (2026-10-09, macOS `sample` on the Release .app, heavy scene):** ~75% of the main thread idle in limitFrameRate
  (SDL_Delay + busy-wait), real work ~2.6 ms/frame of 16.7. Physics 18% of samples, of which 82% UniformParticleGrid::queryCandidates;
  render ~5% (122 samples of memmove inlined in RenderFullFrame = the 2.4 MB star-field upload); gravity only 39 samples.
  `sample <process name> 10 -file profile.txt` needs the game RUNNING (a path to the .app fails with 'unable to read input graph').
  Mac has only Command Line Tools (no Xcode/Instruments/xctrace). Sampling hides single-frame spikes -> Tracy or worst-frame timers.
- **Star field + particle grid speedups (2026-10-09, `Bugfixes`, typed by Marco):**
  - Star field uploaded ONCE in the RenderSystem constructor (the existing init command buffer), no longer per frame (2.4 MB/frame).
    It's static: twinkle = time uniform in the frag shader, parallax = camera constants.
  - UniformParticleGrid::queryCandidates: column stencil. Own cell + (0,1) [next row, cy+1 = BELOW on screen: y grows downward,
    Marco corrected Claude's 'above'] are walked with one firstEntryAtOrAfter + forward scan; next column (1,-1),(1,0),(1,1) =
    3 consecutive keys = one search + scan. 2 binary searches per particle instead of 10. Relies on packCell's layout (comment on
    both sides). Verified IDENTICAL pair sets vs the old grid on 6 scenes (`Notes/drafts/Harnesses/grid_equivalence_harness.cpp` +
    OldParticleGrid.*), 2.2-3.2x faster (e.g. 5.5 -> 2.2 ms/tick at 20000 particles). Benchmark showed queries were ~90% of grid time,
    candidates only 0-3/particle (so the giant-cell suspect was NOT the cause in tested scenes).
  - Key aliasing beyond +-1e6 cells (>= 2M px from origin, ~10-15 fully zoomed-out screens) only adds candidates the distance check
    rejects; never misses a real neighbour. Not a bug. Idea for Marco: despawn far-away particles to return budget.
- **Ship sprite (2026-10-09, `Bugfixes` after commit 58742d5, typed by Marco in walkthrough mode, seen working in-game):**
  Marco's Procreate drawing exported as `Transfer/Assets/Visual/Ships/FutureShip.png` (2048x2048 RGBA, transparent, nose up,
  centred; the .procreate source stays OUT of Assets/ because everything there is bundled). RenderSystem::loadSpriteTexture(path
  under Assets/): SDL_LoadPNG (built into SDL 3.4, no SDL_image) -> SDL_ConvertSurface RGBA32 -> SDL_PremultiplySurfaceAlpha ->
  texture R8G8B8A8, usage SAMPLER|COLOR_TARGET (COLOR_TARGET is REQUIRED by SDL_GenerateMipmapsForGPUTexture: verified by SDL's
  debug assert), num_levels = halvings + 1 (2048 -> 12) -> upload level 0 -> generate mipmaps outside any pass. Returns nullptr +
  logs if missing; renderStarship then skips (Marco's choice). m_sprite_sampler: linear + mipmap LINEAR + clamp, max_lod = 1000
  (default 0 = never leaves level 0). Starship::buildGeometry = one textured square (6 verts, u right / v down, v=0 = nose) of side
  shipSize (Marco set 150; the 50 placeholder looked tiny), centred where the old box was. Shaders: vert passes uv; frag samples
  t0/s0 space2. Pipeline: frag 1 sampler, blend src = ONE (premultiplied). Draw order now stars -> bodies -> ship -> arrow -> UI
  (Marco: ship above bodies; he plans physics so they never overlap). Old buildGeometry kept commented out by Marco (to delete).
  Open quiz: what the ship would look like zoomed in vs out without the mipmap generation call.
- **Ship physics, DONE (2026-10-09, `Bugfixes`, typed by Marco step by step, 'works great' in-game).** Marco's decisions: convex POLYGON hitbox, BOUNCE on impact (pushed
  out = never overlaps, loses energy), ONE-WAY gravity (planets pull the ship, not back), rotation player-only (no spin physics).
  Plan: (1) DynamoEngine/Physics/Collision2D: `CircleContact circleVsConvexPolygon(center, radius, polygon)` -> {touching,
  normal (polygon -> circle, unit), depth}; inside-centre case = nearest edge's outward normal, depth r + d; outside r - d;
  exactly touching = not touching; outward normals flipped using the corner average (any winding). DONE: typed by Marco,
  Tests/DynamoEngine/Physics/Test_Collision2D.cpp 9 tests (Claude wrote the last 3 at Marco's request and FIXED his
  squareReversed(), which wasn't reversed: the winding mutation then went uncaught). Mutation-checked, TidyEngine 36 files clean.
  DECISION CHANGE (same day): gravity is TWO-WAY after all (Marco: 'I do want actual mass'); ship DRY_MASS = 1000 (m_mass,
  fuel goes on top later). Marco first thought the ship needed mass to orbit: it doesn't (a = GM/r^2, its mass cancels) --
  mass matters for pulling planets back, collisions and thrust. LONG-TERM (Marco): Tsiolkovsky rocket equation with fuel mass
  for main-engine and RCS burns: simulate it (fuel -= mass_flow*dt; accel = thrust_force / current mass), the equation itself
  only for UI (delta-v remaining). Seam = Starship::applyThrust. RCS ROTATION would reverse 'player-only rotation' (needs
  angular velocity/torque/inertia): a conscious later decision.
  FINAL DESIGN (supersedes parts of the plan below): (3) debug overlay done (F3, green ship polygon + yellow planet circles,
  DebugLine shaders listed in BOTH MakeTransfer scripts). (4) Starship: m_mass = DRY_MASS 1000; applyThrust sets
  m_thrust_acceleration = nose x THRUST_ACCELERATION (1200 px/s^2 = old feel); KDK Verlet applyVelocityVerletPhase1/2 inside the
  bodies' phases; PhysicsSystem::updateShipGravity (in updateAllForces) = TWO-WAY ship<->planet forces (same softening rule as
  calculateGravity), setGravityForce -> a = F/m + thrust. Debris and ship don't attract. (5) handleShipCollisions (in
  handleCollisions, before the pending-fragment flush): planets + particles (bounding-radius pre-check), resolveShipContact =
  push apart split by inverse mass (force-static = 0), impulse with SHIP_RESTITUTION = 0.2 only when closing, returns absorbed
  energy 1/2 mu v^2 (1-e^2) -> ship.setLastImpactEnergy (hook for future DAMAGE; Marco's plan). NO ship-caused shattering
  (Marco: planets will be scaled way up; strong gravity made any fall exceed MIN_SHATTER_SPEED). Marco renamed all Starship
  members to m_ style (m_position, m_prev_position, m_rotation, m_ship_size, m_velocity + velocity() getter).
  Harnesses in Notes/drafts/Harnesses: hitbox_, ship_physics_ (thrust/orbit/momentum), ship_collision_ (energy audit exact,
  bounce, no shatter, debris momentum, resting); they use `#define private public` to place the ship.
  (2) ship shape: 6 points traced from the PNG's convex hull, as fractions of the image x shipSize: (32,1464) (959,366)
  (1094,366) (2021,1464) (1242,1744) (811,1744) of 2048 px (96.4% of hull area). (3) debug lines + F3 toggle wiring the unused
  UIState::renderDebug (default VIEW_DEBUG): ship polygon green, planet circles yellow, no particles; own LINELIST pipeline with
  prev+current positions (interpolated like the sprite). (4) ship mass, one-way gravity, Verlet, thrust = accel x dt.
  (5) collisions vs macros + particles: push out + bounce (ELASTIC_LOSS_FACTOR).
- **Delete keys (2026-10-09, `Bugfixes`, typed by Marco):** Delete/Backspace (a Mac's 'delete' key IS SDL_SCANCODE_BACKSPACE)
  = undo the newest spawn, repeatable; Shift + Delete = clear everything (the old behaviour). PhysicsSystem::removeNewestSpawn:
  highest macroIdentifier still carried by any body, std::erase_if on macros + particles (immediately, so it works with time
  stopped); no history stored. Clusters now get a spawn ID too (fixed the shared -1 shader seed). Harness:
  Notes/drafts/Harnesses/undo_spawn_harness.cpp. Max planet radius was raised by Marco ('bigger planets'), nothing else needed.
- **UICheckbox (2026-10-09, `Bugfixes`, header + .cpp typed by Marco; Test_UICheckbox.cpp (9 tests) and the TestVisual demo
  checkbox written by Claude at Marco's request):** engine widget = UIButton look/feel (128/108/88 gray, click on release
  inside) + m_is_checked; ToggledAction(bool new_state), CheckedSource pulled in update() (dropdown pattern), setChecked doesn't
  report. Box = 60% of row height, white fill inset 20%, label after one gap. New UISound::Checkbox (Marco's name) mapped in
  TransferScenes queueUISound. AudioSystem auto-loads every .wav in Assets/SoundEffects by file STEM. Mutation-checked, tidy
  38 files clean. NEXT (phase 2): persistent spawn-settings struct (populateGravBodyProperties... hard-codes flags; the
  isCreating* flags are reset by resetTransientFlags after every spawn), a panel of checkboxes bound via source + action.
- Physics fix (2026-09-27, applied by Claude at Marco's request; committed in `26a094a`): shatter fragments go into
  `PhysicsSystem::m_pending_fragments` (via a `fragments_out` parameter on substituteWithParticles / ...FromImpact) and are appended
  to `particles` after all three collision passes (`handleCollisions`); `createParticleCluster` passes `particles` directly (no loop
  running, must appear with time stopped). `liveParticleCount()` = particles + pending for both MAX_LIVE_PARTICLES checks.
  The `reserve(+32000)` hack and MAX_SIMULTANEOUS_SHATTERS_PER_TICK are gone. Verified with an ASan crash test
  (`CLAUDE/drafts/PendingFragments/stress_harness.cpp`, since lost). Marco first tried a member-only variant and
  flushed inside substituteWithParticles, which re-created the bug. Lesson: only empty the waiting list where no loop over particles is running.
- Sub-pixel fade (2026-09-27) REJECTED by Marco for now ("not sure it's necessary, needlessly complex"); he's restoring the original. The draft (`CLAUDE/drafts/SubPixelFade/subpixel_fade.diff`) is lost; this description is enough to redo it.
  CameraConstants `_padding1` -> `pixelDensity` (set from SDL_GetWindowPixelDensity; same name in all 5 cbuffers). Vertex shader: d = 2*r*zoom*density;
  d < 1 -> draw a 1-px quad, `nointerpolation float coverage = d*d` (TEXCOORD15, last in both structs; MSL locn8 flat); frag skips the circle
  discard when coverage < 1 (square = no flicker) and alpha *= coverage. CPU: uploadUnifiedBodies skips d < 0.0627 (= sqrt(1/255), fade is
  invisible); renderBodies draws unifiedBodyVertices.size() (fixes the recount REWORK item). Physics stays tied to MIN_PARTICLE_RADIUS:
  canBreakIntoSurvivableFragments (0.75R >= MIN) guards shatter / mutual shatter / crumble / createParticleCluster (replaces the magic 1.0);
  cleanupParticles radius check kept as a commented safety net. Measured: impact shatters of tiny bodies did NOT lose mass today
  (FromImpact's radius growth rescues the single fragment, 0/500 trials) -- the guard matters for raising MIN (load balancing) and for the
  cluster path (radius 1.2 click spawned 1 particle that died next tick). Behavior change: bodies with 0.75R < MIN now bounce instead of
  becoming one particle. Visual result NOT verified by Claude (can't see screen): Marco must zoom out on a shatter (Retina + Windows 1x).

**Next candidates** (let Marco pick):
- Build/scripts follow-ups (the rework itself is DONE, see the done list): the [BUG]/[REFACTOR] items in REWORK's BUILD / TOOLING
  section (Windows setup: silent rmdir failure, unchecked submodules, implicit job count; Mac arch / minimum macOS version;
  Windows VC++ runtime).
  Git gotcha: restaging on Windows drops the Mac scripts' executable bit (100755 -> 100644); fix with
  `git update-index --chmod=+x Scripts/Apple/*.sh` before committing.
- Visor styling: row gaps, a selected-row look, a real "Realistic" coloring, then the actual Charge/Temperature views.
- E: action map, i.e. rebindable keys that replace hard-coded scancodes. It will also retire DEPRECATED_InputState. All key bindings
  live in `translateGameInputs` on purpose, to make E easy.
- G: keyboard focus (Tab will then conflict in menus; the visor hotkey only runs in simulation scenes, on purpose).
- H: world-anchored UI.
- R1–R4: rendering into the engine.
  - R1 GPU upload helper (fixes cycle=false).
  - R2 pipeline helper (about -300 lines).
  - R3 engine UIRenderer.
  - R4 RenderContext.
  - SDLContext lifetime fix (SDL_Quit runs before members die; see REWORK).
- REWORK.txt: ~33 items (physics bugs, render risks, audio, build/tooling).

---

## 6. Learning / quiz log (what to revisit)
- Round 1 (2026-09-24, about 46/100). Strong: loop order, softening, momentum conservation. Weak:
  - velocity Verlet KDK order and why symplectic matters
  - where alpha comes from
  - `std::vector` reallocation invalidating references (important)
  - UniformParticleGrid half-stencil (another session wrote it)
  - Tends to answer "what it should do" rather than "what the code does".
- Since then he has understood well:
  - edges vs levels
  - deferred scene switching
  - why `[&x]` vs `[x]` matters (by-value capture freezes the value)
  - `explicit` on one-argument constructors
  - `static_cast` direction (after one slip)
  - `enum class` vs a plain enum in a namespace (type safety vs generic SceneId)
  - override signatures must match the base (unused params)
  - `containsPoint` decides whether a press reaches an element; `onMousePressed`'s return decides whether it stops there
  - ordering bugs (`optionAt` before `m_is_open = false`)
- 2026-09-28/29 build work: understood why the test exe died with 0xc0000135 (DLL not found) and why Windows vs Mac `--parallel`
  differed (MSBuild caps at CPU count, bare `make -j` is unlimited). Quiz "why cap jobs at RAM/2 too": mostly right (10 jobs x ~1 GB
  > 8 GB); missed that the OS needs RAM too, and that CPU count bounds speed while RAM bounds safety (take the min). Still open, never
  answered: why `NO_DEFAULT_PATH` (a leftover SDL elsewhere, e.g. found via PATH-derived prefixes, could win over ThirdParty); what
  happens if the DLL copy goes AFTER gtest_discover_tests; why `if errorlevel 1` rather than `if %ERRORLEVEL% NEQ 0` inside a ()
  block (%...% is expanded once, when the whole block is parsed).
- Test-writing habits taught: use fixture helpers instead of hand-computed coordinates; `EXPECT_EQ(a, b)` rather than `EXPECT_TRUE(a == b)`
  (the failure prints both values); one assertion per promise in the test name; break the code on purpose and watch the test go red.
- Open physics question he raised: runaway spinning of particle clumps. Hypotheses: sequential impulse order bias, elastic re-bounce
  without an approaching check, positional correction adding energy.
- 2026-09-27: With MIN_PARTICLE_RADIUS temporarily at 40, fragments looked like self-attracting "grav bodies". Particle-particle gravity is
  commented out, and particles CANNOT merge (handleAccretion returns early when heavier.isParticle; REWORK [BUG]; Marco caught Claude's wrong
  "they merge" claim). Real cause: handleElasticCollisions' dynamic branch has NO approaching check (the static branch has `if (v_n < 0.0)`),
  so born-overlapping fragments get their separating normal velocity swapped back toward each other every tick -> they stick/jiggle = looks
  like attraction; ratio >= 8 pairs ghost through. Fix next session = one approaching check (REWORK handleElasticCollisions [BUG]); also the
  prime suspect for the open "runaway spinning clumps" question.
- 2026-10-09 quizzes (Bugfixes): spiral of death ~5/10 (said the backlog is 'spread over' 8 steps -- it's DROPPED; 0.5s = 60 steps
  not 50; missed the snowball). Approaching speed ~8/10 (called a dot product 'in the x direction' -- it's a scalar; caught Claude's
  quiz missing radii). isMarkedForDeletion break 10/10. Dispatcher vs handleAccretion 8/10. +1 in buffer size 5/10 (said 0-based
  indexing; it's the preview body -- counts never need +1). Budget arithmetic 9/10. Doubling loop 6/10 (counted 3 runs, it's 2; the
  start value isn't a run). Marco corrected Claude twice: 'particles can't accrete' (absorb) is right; bouncing at the limit is wrong.
- 2026-09-28 FIXED (applied by Claude at Marco's request; committed in `4737715`): handleElasticCollisions dynamic branch now only exchanges normal
  velocities when closing_speed = v_smaller_n - v_larger_n > 0 (normal points smaller -> larger); positional correction still always runs.
  Verified: overlapping separating pair keeps -50/+50 (before: flipped every tick); head-on pair bounces once; 120/120; ASan stress clean
  with MIN=1. (MIN_PARTICLE_RADIUS is back to 1.0 as of `4737715`; at 40, cleanupParticles also deleted any particle under radius 40.) New REWORK [BUG]: getCollisionInfo's should_blow_up uses |normal speed|.

