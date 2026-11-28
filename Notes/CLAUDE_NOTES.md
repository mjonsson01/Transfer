# Transfer: Claude handoff notes

**Fresh Claude session, on any device: read this whole file, then `Notes/REWORK.txt`, before doing anything.**
This folder is the only memory that travels between machines (Claude's own memory is per-device).
Last full rewrite: 2026-09-27, on macOS, branch `EngineSep` at `959b486` ("Finished initial full implementation of visor_view").
Updated 2026-09-29 on Windows, at `7f852af` (everything committed and pushed; build/scripts rework done on both platforms).
Updated 2026-10-09 on macOS: EngineSep merged to master (PR 31), Mac packaging merged (PR 32/33, master at `a1cb3a7`).
The folder was renamed `CLAUDE/` -> `Notes/`. Current working branch: `Bugfixes` (REWORK items first, then new features).
Updated 2026-10-10 on macOS, end of a long session: Bugfixes merged to master up to PR 35 (`ba10f43`). UNCOMMITTED at handoff:
the cluster spawn preview (FragmentLayout.hpp/.cpp are new, untracked -- `git add` them) + REWORK edits. 139/139 tests, tidy 38
clean. Section 5's done list has one entry per feature from 2026-10-09/10 (newest physics ones near the top of 'Done').
Since 2026-10-09 Marco sometimes asks Claude to APPLY changes directly ("apply it"); the default is still walkthrough
(he types). Verify in a scratch copy first either way. Tip: compile Transfer/src to .o once and link harnesses against it.
Updated 2026-10-11: the cluster preview is committed (`d41ec3d`). New branch `IntegrationLines`: planning mode + trajectory
tracer, DESIGN AGREED, nothing typed yet (see 'IN PROGRESS' at the top of section 5's 'Next candidates').

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

**SOURCE LAYOUT (reorganized 2026-11-25 by Claude, Marco's choice of layout; older notes use the OLD paths):**
```
Transfer/src/
  main.cpp
  Core/              Game, GameState, UIState                       (only the top level of the game)
  Systems/
    Physics/         PhysicsSystem, TrajectoryPrediction.cpp, VelocityVerletIntegration, UniverseState,
                     Trajectories, SpawnSettings, GravitationalBody(+Pair), FragmentLayout, UniformParticleGrid,
                     PhysicsConstants
    Rendering/       RenderSystem, GPUTypes, CameraData (the shaders' CameraConstants cbuffer), TwinklingStars
    Input/           InputSystem, PlayerInput
    Audio/           AudioSystem, MusicModeEnum, SoundEffect
  Camera/            CameraState, CameraTransform (screen <-> world)
  Player/            Player, Starship
  Scenes/            TransferScenes
  Constants/         EngineConstants, GameSystemConstants
  Utilities/         CustomMathUtilities, SystemPathUtility
  HLSL/, DynamoEngine/   (unchanged)
```
Old -> new: Core/{CameraState} -> Camera/; Core/{PlayerInput} -> Systems/Input/; Core/{SpawnSettings,Trajectories,
UniverseState} -> Systems/Physics/; Entities/Physics/* -> Systems/Physics/; Entities/Sound/* -> Systems/Audio/;
Entities/VisualElements/TwinklingStars -> Systems/Rendering/; Systems/XSystem -> Systems/X/XSystem (Render ->
Rendering/); Utilities/Constants/{Engine,GameSystem}Constants -> Constants/, PhysicsConstants -> Systems/Physics/;
Utilities/Physics/* -> Systems/Physics/; Utilities/Rendering/{CameraData,GPUTypes} -> Systems/Rendering/,
CameraTransform -> Camera/; Utilities/Math/ and Utilities/System/ flattened into Utilities/. Done with git mv + a
script rewriting every include (src, Tests, harnesses); verified: deterministic harnesses identical, Release build
0 warnings, 148/148 tests, TidyEngine clean. Documentation/class_map.md is stale (rerun PythonScripts/create_docs.py).

**Game loop** (`Core/Game.cpp`, `Game::Run`): frame_seconds → updateFPS → `ProcessInput(frame_seconds)` → UpdateInstantiations
(spawns from input) → PlayAudio → fixed-step physics (120 Hz accumulator × time scale, only while `scenes.currentScene().settings().runs_simulation`)
→ alpha → RenderFrame → **`scenes.applyPendingSwitch()`** (deferred scene switches happen here) → limitFrameRate (60 FPS).
`Game` owns GameState, UIState, InputSystem, PhysicsSystem, RenderSystem, AudioSystem, `DynamoEngine::SceneManager scenes`.

**Input** (`Systems/InputSystem`): SDL events → `DynamoEngine::SDLInputIntake` (free fn `translateSDLEvent`) → `DynamoEngine::InputState`
(edges vs levels, `beginInputFrame`, focus-lost clears held keys) → **UI first**: `updateSceneUI` does setWindowSize, `updateElements(dt)` and
`processInput` (the order matters), giving a `UIInputResult`. Then the game side (a release on the UI, pointer_captured, spawns nothing):
`translateGameInputs(game_state, ui_state, scenes, player_input, pointer_captured)` / `translateMenuInputs(scenes)`, writing into
**`Core/PlayerInput`** (owned by Game, passed explicitly; since 2026-11-25, replacing DEPRECATED_InputState): ShipControls,
is_aiming_spawn + SpawnProperties (world space, computed once in InputSystem::aimSpawn), one-frame PlayerCommands cleared every frame.
Key bindings: Esc = pause/resume (requestSwitch), Delete = undo newest spawn, Shift+Delete = clear, W/S = throttle up/down
(holds), A/D = turn, P = planning, I = trace (planning), F3 = hitboxes, F4 = profiler, left/right release = spawn macro/cluster,
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
- `UIRoot`: UI scale = `UIRoot::fitScale(window)` × player scale, fitScale = min(width / `UI_REFERENCE_WIDTH` 1280, height / `UI_REFERENCE_HEIGHT` 720)
  (since 2026-11-26; was height only), so the 1280 x 720 design space ALWAYS fits: extra space appears at the sides (ultrawide) or top/bottom
  (16:10, Mac fullscreen). fitScale is static so world overlays can use it too. It works in UI space: `uiSpaceSize()` and `screenToUISpace()`.
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
    The sliders write into SpawnSettings (radius, mass) and UIState::setSimulationSpeed (since 2026-11-25).
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
- Particle broadphase: `UniformParticleGrid`: a HASHED grid since 2026-11-26 (cells folded into 2 slots per particle,
  built by counting sort, cell lookup = one read of m_slot_start; entries carry their cell so slot-sharing cells are
  skipped). queryCandidates(index, out): the particle list parameter was dropped (cells are stored at build).
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
  centered; the .procreate source stays OUT of Assets/ because everything there is bundled). RenderSystem::loadSpriteTexture(path
  under Assets/): SDL_LoadPNG (built into SDL 3.4, no SDL_image) -> SDL_ConvertSurface RGBA32 -> SDL_PremultiplySurfaceAlpha ->
  texture R8G8B8A8, usage SAMPLER|COLOR_TARGET (COLOR_TARGET is REQUIRED by SDL_GenerateMipmapsForGPUTexture: verified by SDL's
  debug assert), num_levels = halvings + 1 (2048 -> 12) -> upload level 0 -> generate mipmaps outside any pass. Returns nullptr +
  logs if missing; renderStarship then skips (Marco's choice). m_sprite_sampler: linear + mipmap LINEAR + clamp, max_lod = 1000
  (default 0 = never leaves level 0). Starship::buildGeometry = one textured square (6 verts, u right / v down, v=0 = nose) of side
  shipSize (Marco set 150; the 50 placeholder looked tiny), centered where the old box was. Shaders: vert passes uv; frag samples
  t0/s0 space2. Pipeline: frag 1 sampler, blend src = ONE (premultiplied). Draw order now stars -> bodies -> ship -> arrow -> UI
  (Marco: ship above bodies; he plans physics so they never overlap). Old buildGeometry kept commented out by Marco (to delete).
  Open quiz: what the ship would look like zoomed in vs out without the mipmap generation call.
- **Ship physics, DONE (2026-10-09, `Bugfixes`, typed by Marco step by step, 'works great' in-game).** Marco's decisions: convex POLYGON hitbox, BOUNCE on impact (pushed
  out = never overlaps, loses energy), ONE-WAY gravity (planets pull the ship, not back), rotation player-only (no spin physics).
  Plan: (1) DynamoEngine/Physics/Collision2D: `CircleContact circleVsConvexPolygon(center, radius, polygon)` -> {touching,
  normal (polygon -> circle, unit), depth}; inside-center case = nearest edge's outward normal, depth r + d; outside r - d;
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
  m_thrust_acceleration = nose x THRUST_ACCELERATION (1200 px/s^2 = old feel); KDK Verlet applyVelocityVerletFirstHalf/2 inside the
  bodies' phases; PhysicsSystem::addGravityBetweenShipAndBodies (in updateAllForces) = TWO-WAY ship<->planet forces (same softening rule as
  addGravityBetweenGravBodyPair), setGravityForce -> a = F/m + thrust. Debris and ship don't attract. (5) handleShipCollisions (in
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
  report. Box = 60% of row height, white fill inset 20%, label after one gap. New UISound::StandardClick (Marco's name) mapped in
  TransferScenes queueUISound. AudioSystem auto-loads every .wav in Assets/SoundEffects by file STEM. Mutation-checked, tidy
  38 files clean. NEXT (phase 2): persistent spawn-settings struct (populateGravBodyProperties... hard-codes flags; the
  isCreating* flags are reset by resetTransientFlags after every spawn), a panel of checkboxes bound via source + action.
- **Spawn panel + Static (2026-10-09, `Bugfixes`, typed by Marco):** Core/SpawnSettings.hpp (persistent, NOT reset per spawn;
  owned by UIState::spawnSettings()), passed into createMacroBody/createParticleCluster/populate... Panel = UIRow of
  UICheckbox at BottomCenter (margin 20) under the sliders; SLIDER_ROW_BOTTOM_MARGIN is DERIVED (panel margin + checkbox
  height + 50 gap) so they can't overlap; each checkbox = setCheckedSource + setOnToggled on a SpawnSettings field
  (~5-6 fit across 16:9 at 200 wide). Marco's MEANING of isForceStatic: no forces act on it; it may shatter, absorb and be
  absorbed (other flags + normal mass ratios decide); in momentum exchanges it is INFINITELY heavy: bounces/ship = wall;
  absorbing keeps its velocity; a body that absorbs a static one takes its velocity AND becomes static (option A). Static
  drifts if spawned with velocity. Spawned clusters: every particle static (loop in createParticleCluster). Shatter/crumble
  debris of a static body is NOT static (Marco tried it, disliked the behaviour, reverted). Flags stay orthogonal: one
  checkbox = one flag (isBounce '&&' bug logged in REWORK). Harness: Notes/drafts/Harnesses/static_spawn_harness.cpp
  (its case 3 expects static debris -> now expect 0 static). CheckboxToggle.wav raised +12 dB (peak -1.5 dBFS).
- **Spawn flags round 2 (2026-10-09/10, `Bugfixes`):** ONLY PLANETS can be static (Marco: static particle clouds behaved badly);
  SpawnSettings field renamed is_force_static. handleElasticCollisions static branch uses RELATIVE velocity
  (dyn.v - stat.v) and pushes the other body fully out (drifting statics passed through resting bodies). isBounce check is
  now `||` (one bounce body is enough; checked before shatter, so bounce bodies never shatter). New SpawnSettings
  is_shatterable / is_accretable (default true = old hard-coded) / is_bounce (false); clusters copy accretable + bounce to
  every particle (not static, not shatterable). Static vs static (and not merging): bounce AS EQUAL MASSES (half the push,
  swap normal velocities x ELASTIC_LOSS_FACTOR), both stay static (Marco's choice; they used to pass through).
  Spawn panel = UIRow of two titled sections (UIColumn: UILabel title 27 tall + UIRow of 200x36 checkboxes):
  'Global Instantiation Properties' (Accretable, Bounce) and 'Macro Body Instantiation Properties' (Static, Shatterable);
  helpers addSpawnCheckbox(row, label, bool&) and makeSpawnSection(title, row). Checkbox labels get 164 of 200 points
  (18 pt SpaceMono ~11 pt/char, 'Shatterable Macro' = 187 did not fit). UIColumn lives in UIRow.hpp.
  Harnesses: spawn_flags_, static_bounce_, static_spawn_, layout_ (prints frame 1 and frame 3).
- **Nested layout fix (2026-10-10, APPLIED BY CLAUDE at Marco's explicit one-off request, with a test):** new
  `virtual void UIElement::updateSize() {}`; UIRow/UIColumn override it: children's updateSize first, then fit around them;
  their updateLayout now starts with updateSize(). Before, a nested row/column reported LAST frame's size (0 on frame 1),
  one frame of lag per nesting level (the spawn panel sat below the screen for 2 frames). Test
  UIRow.NestedContainersAreRightAfterOneUpdate (Test_UIWidgets.cpp): row > column(title + row) > boxes, checked after ONE
  update; failed before the fix (outer row 90 wide instead of 300), passes after, fails again if the children-first loop is
  removed. 139/139, tidy 38 clean. Small known cost: a nested container's size is recomputed once per ancestor level
  (tiny UI trees, irrelevant).
- **Shatter order + stacked bodies (2026-10-10, `Bugfixes`):** (a) Shatter decided PER BODY (typed by Marco):
  lighter_shatters = lighter.isShatterable; heavier_shatters = heavier.isShatterable && ratio <= MUTUAL_SHATTER_...; shatter if
  shouldBlowUp && (either). Before, everything was gated on the lighter one, and for EQUAL masses 'lighter' = spawned second
  (pickMassPair returns {a, b} on a tie), so a shatterable set-up body hit by a non-shatterable one just bounced. A surviving
  non-shatterable planet then absorbs the debris it flies through (accretion, ratio ~800). (b) Exactly stacked bodies (same
  click spot, or two identical clusters = same Vogel positions) stuck forever: every handleElasticCollisions branch returned
  on distance 0. Helper `directionFromTo(from, to)` = unit vector, or a RANDOM unit vector when the points coincide; used in
  all three branches (APPLIED BY CLAUDE at Marco's request; Marco had started it -> Claude's first script inserted a duplicate
  helper, then cleaned it up and finished his half-typed both-static branch). Mostly-overlapped bodies were always fine
  (planets separate in ~10 ticks; 1-px-offset identical clusters take ~2 s: crowding). Harnesses: stacked_spawn_,
  overlap_spawn_, shatter_order_. Tip: compile Transfer/src once to .o files and link harnesses (per-harness full
  recompiles timed out at 600 s).
- **Heavy bodies never separated (2026-10-10):** the dynamic push-apart in handleElasticCollisions was guarded by
  `totalInvMass > EPSILON` (1e-8): every pair heavier than ~2e8 each was SKIPPED (stacked 3e8 planets stuck at full overlap).
  Marco changed it to `> 0`; Claude proposed the same-sign version (shares by inverse mass only when invMass_a * invMass_b > 0,
  else 0.5/0.5) so negative-mass pairs also separate without the shares exploding when inverse masses nearly cancel.
  Lesson: a fixed epsilon only makes sense at a known scale (masses span 1..1e10). Marco's decision: very heavy planets that
  touch may shatter right after separating (two 1e9 r=60 planets gain ~414 px/s closing speed in one tick from gravity, >
  MIN_SHATTER_SPEED) -- fine (option B). Ideas kept for later: shatter threshold relative to mutual escape speed; resting
  contact bounces forever (ELASTIC_LOSS_FACTOR 0.994). Harnesses: placement_separation_, heavy_overlap_.
- **Cluster spawn preview (2026-10-10, applied by Marco, works in-game):** new free functions in
  Utilities/Physics/FragmentLayout.hpp/.cpp: survivableFragmentCount (was a PhysicsSystem member) + layOutAndPopulateFragments (the Vogel
  spiral, no velocities, body untouched). substituteWithParticles = layOutAndPopulateFragments + random velocities + mark deleted;
  RenderSystem::appendPreviewBodies draws layOutAndPopulateFragments(preview_body) for a cluster, so the preview is EXACTLY the real
  cluster (harness: 0 difference). DEPRECATED_InputState::isPreviewingCluster (right button, on world; left wins if both);
  shift adds the velocity arrow for either. Right-click had NEVER previewed (checked history back to April) -- Marco thought
  it regressed. EngineConstants MAX_PREVIEW_BODIES = 800, INITIAL_UNIFIED_BODY_CAPACITY = MAX_LIVE_PARTICLES +
  MAX_PREVIEW_BODIES; PhysicsConstants static_asserts DEFAULT_FRAGMENT_COUNT <= MAX_PREVIEW_BODIES. randomDouble is seeded
  from std::random_device: harness output involving fragments/stacking varies run to run (compare ranges, not exact text).
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

**IN PROGRESS (2026-10-11, branch `IntegrationLines`): planning mode + trajectory tracer. Design agreed in discussion; next
step = walkthrough, starting with the tracer class header.**
- Keys (no checkboxes, Marco dropped them): **P** toggles PLANNING mode (game scene only). **I**, only while planning, toggles
  tracing: start, or clear the lines. I outside planning = ignored.
- Planning = forced pause. Effective time scale in Game::Run = planning ? 0 : slider value. The slider is never overwritten, so
  leaving planning resumes at whatever the slider shows (dragging it while planning changes the speed you return to; Marco OK'd).
  No save/restore, no UISlider value source needed. Ship controls run inside the physics tick, so the ship is frozen too.
- Tracer = a separate copy of the state (macros + ship + debris) stepped with the REAL functions in the same order
  (Verlet phase 1, updateAllForces, phase 2): NO collisions, no cleanup, no player input. Contacts phase through (Marco: contact
  visualisation later). Same code + same order = the line is exactly what the game does while nothing un-modelled happens
  (incl. the skipped first half kick of a body spawned while paused). Idea: the copy is a GameState with the same vectors; a new
  public PhysicsSystem step does only the three gravity/Verlet calls.
- Debris is IN the copy (it pulls planets: addGravityBetweenGravBodyPair(particle, macro) is two-way) but NEVER gets lines (Marco
  2026-10-11: particle paths dropped, 20000 x points is too much data). Lines = macro bodies + ship only.
  Particles don't collide in the prediction -> least accurate next to a fresh shatter.
  Shatter randomness (randomDouble, shared generator) is why the prediction must NOT run collisions.
- Drifting static bodies get a line (free: Verlet drifts them). Rewind = separate later feature (needs collision history).
- Growth: the line is TRACED, not instant: accumulator like Game::Run, TRAJECTORY_TRACE_SPEED (game s per real s, tunable
  constant, Marco wants 'as fast as looks smooth') + a max steps per frame cap (spiral-of-death lesson). Record a point every
  few ticks. Horizon: long, still a constant to tune (memory needs SOME max; stops when reached).
- Wipe: any edit of the real state while tracing/done wipes the lines and WAITS for I again. Plan: detect by comparing the real
  macros + ship (exact ==, safe because nothing changes while paused) and particle count to the snapshot the trace started
  from, instead of calling invalidate() from every edit site (covers future velocity dragging automatically).
- Look: transparent front, solid back (alpha per vertex). SDL_GPU lines are 1 px; thicker = quads later.
  Status readout label (e.g. 'Integrating...').
- Future (Marco): drag initial velocity vectors while planning, guided by the lines.
- Decided 2026-10-11: each new piece of line FADES IN by age (finished lines end fully solid), and while tracing the TIP
  stays brighter than the trail so the front is obvious. Leaving planning (P) makes the lines vanish. Readout = UILabel
  top-left under FPS: 'PLANNING: press I to integrate' (no lines) / 'PLANNING: Integrating...' (tracing) / 'PLANNING'
  (finished); empty outside planning.
- ARCHITECTURE (2026-10-11, after 3 rejected variants: Game-owned tracer / tracer in PhysicsSystem holding a GameState /
  tracer in GameState = recursion). Root cause: GameState is a bag of everything and physics takes GameState&.
  PHASE 0 = refactor first (Marco: 'as early as we can'): a plain copyable struct **UniverseState** (Marco's name; not
  'World', it's a space sim) = what physics simulates: macro bodies, particles, the ship, the spawn-ID counter (proposed).
  GameState owns a UniverseState AND a Trajectories (separately). Physics functions take UniverseState&, not GameState&.
  NO separate PredictionSystem (Marco: prediction = advancing a subset of physics -> physics owns it; Claude conceded).
  NAMES (Marco, explicit on purpose; he finds the codebase's names blur together):
    PhysicsSystem::advanceUniverseState(universe, controls)   = one real tick (replaces UpdateSystemFrame), calls:
      handleCollisionEvents -> applyShipInputPhysics -> integrateUniversePhysics -> cleanupBodies
    integrateUniversePhysics = Verlet first half (half kick + drift), all forces, Verlet second half. Shared by the real
      tick and the prediction (sub-function names proposed by Claude, pending: integrateVelocityVerletFirstHalf,
      calculateNetForces, integrateVelocityVerletSecondHalf).
    advanceIntegratedTrajectories = the prediction's wrapper (Marco): integrator only, makes 'just the integrator' vs
      'the full system' obvious at the call site.
  FILES: Physics owns it all, but the new code goes OUT of PhysicsSystem.cpp (Marco). Claude's proposal (pending): the
    integration rules use no PhysicsSystem member data -> free functions in their own file; a TrajectoryPredictor class in
    its own file (prediction universe, start snapshot, accumulator, recording), owned by PhysicsSystem, calls them directly.
  Player stays (Marco: progression, stats, later) but NO LONGER holds the Starship: the ship moves into UniverseState.
  UniverseState = { macro_bodies, particles, ship, newest_spawn_id }. Phase 0 split into 0a (UniverseState switch) and
  0b (extract + rename integration), each verified BIT-IDENTICAL (fixed no-shatter scene, N ticks) before/after.
  Logged in REWORK 2026-10-11: PhysicsSystem file split, GameState grouping, naming/hierarchy pass, setUniverseState
  (Marco: starting from ANY universe state is baseline physics, do later).
- PHASE 0 PROGRESS: harness Notes/drafts/Harnesses/universe_refactor_harness.cpp (fixed deterministic scene, 1200 ticks,
  %a hexfloats) + universe_refactor_baseline.txt (output of the PRE-refactor code, d41ec3d). After every step: build the
  tree, run with -DAFTER_REFACTOR, diff against the baseline = must be empty. Only the adapter block may change.
  0a-1 DONE: Core/UniverseState.hpp (typed by Marco). 0a-2 DONE (typed by Marco, harness IDENTICAL): GameState owns
  `UniverseState m_universe_state` with getUniverseState()/getUniverseStateMutable() (Marco kept GameState's get-style,
  consistent within the class), old body/ID getters forward to it (temporary), Player no longer holds the ship
  (8 call sites -> getUniverseStateMutable().ship). 0a-3 DONE (typed by Marco, ALL groups at once, harness IDENTICAL):
  every private physics function takes UniverseState& universe_state; Marco also renamed UpdateSystemFrame ->
  updateSystemFrame and CleanUp -> cleanUp (Audio/Physics/Render). The spawn path (createMacroBody,
  createParticleCluster, populateGravBodyPropertiesFromDEPRECATED_InputState) still takes GameState, only for the camera.
  0a-4 DONE (typed by Marco, harness IDENTICAL, real CMake Debug build 0 warnings): forwarding getters gone (except the
  unused incrementMaxIDInstantiated/getMaxIDInstantiated: delete), spawn path takes (UniverseState&, const CameraState&),
  advanceUniverseState(UniverseState&, UIState&), createParticle pruned. Marco is ALSO renaming game code to the tidy rules
  as he goes (RenderSystem: m_window, m_gpu, m_unified_body_vertices, renderFullFrame...). His rename was textual: it
  also hit comments ('release the m_window from m_gpu') and a LOCAL (m_window_flags) -> told him to fix those.
  0b-1 DONE (typed by Marco, his own self-consistent names): advanceUniverseState = handleCollisionEvents ->
  applyShipInputPhysics -> integrateUniversePhysics -> cleanupGravBodies. integrateUniversePhysics (still a PhysicsSystem
  member) = integrateVelocityVerletFirstHalfStep -> recomputeAllForces (addGravityBetweenGravBodies,
  addGravityBetweenShipAndBodies, addGravityBetweenGravBodyPair) -> integrateVelocityVerletSecondHalfStep. Starship:
  applyVelocityVerletFirstHalf/SecondHalf. calculateTotalEnergy pruned.
  0b-2 DONE (Marco, 2026-10-12): Verlet = free functions in Systems/Physics/VelocityVerletIntegration.hpp/.cpp
  (applyVelocityVerletFirst/SecondHalfStep per body, integrateVelocityVerletFirst/SecondHalfStep per universe). Lesson
  from the first attempt: leftover member declarations HID the free functions inside member functions -> link errors.
  Forces stay PhysicsSystem members for now (Marco splits the rest of PhysicsSystem later; REWORK item).
  **PHASE 0 COMPLETE (verified 2026-10-12):** game links, universe_refactor_harness IDENTICAL to the d41ec3d baseline,
  all 12 kept harnesses same as before (shatter_order: same outcomes, random debris counts), real CMake Debug build of
  TransferGame 0 warnings. Engine + Tests untouched (no RunTests/TidyEngine needed). Open nits: VelocityVerletIntegration
  .cpp header comment lacks '.cpp', CustomMathUtilities include belongs in the .cpp (+ PhysicsConstants), no comments in
  its header, Starship's ...FirstHalf vs free ...FirstHalfStep.
  PHASE 1+2 IMPLEMENTED BY CLAUDE 2026-11-24 (Marco: 'implement the whole thing', he reviews as we go; uncommitted):
  - Core/Trajectories.hpp (TracedPoint{position, drawn_at}, TracedPath, TraceStatus Idle/Tracing/Finished, Trajectories
    {status, macro_body_paths, ship_path, clock}); GameState::m_trajectories (Marco typed the member); UIState
    getIsPlanning/setIsPlanning (Marco's copy kept, Claude's duplicate removed); DEPRECATED toggleTrajectoriesRequested;
    Starship::cutThrust; PhysicsConstants TRAJECTORY_TRACE_SPEED 10, MAX_TRACE_TICKS_PER_FRAME 32,
    MAX_TRACE_MILLISECONDS_PER_FRAME 4.0, TRAJECTORY_TICKS_PER_POINT 4, TRAJECTORY_HORIZON_SECONDS 60.
  - Systems/Physics/TrajectoryPrediction.cpp = PhysicsSystem members (option a): updateTrajectories (public, once per
    frame) -> startTrajectories / clearTrajectories / universeChangedSinceTraceStart / advanceIntegratedTrajectories /
    recordTrajectoryPoints; members m_prediction_universe_state, m_trace_start_universe_state, accumulator, ticks.
    Horizon ticks use lround(seconds * INV_PHYSICS_TIME_STEP): a plain cast of 60/(1/120f) gave 7199 (harness caught it).
    TIME BUDGET added after measuring 30 ms/frame at 25 planets + 20000 particles -> now ~4.5 ms avg, 6 worst.
  - Game::Run: time_scale = planning ? 0 : slider; updateTrajectories after the physics loop if runs_simulation.
  - InputSystem::translateGameInputs: P toggles planning; I (only while planning) sets the request.
  - RenderSystem: m_trajectory_line_* buffers (MAX_TRAJECTORY_LINE_VERTICES 131072 in EngineConstants), drawn with the
    DEBUG LINE pipeline (same vertex layout), after the stars and before the bodies. Fade-in by age (trajectoryAlpha:
    0.1 -> 0.7 over 0.75 s), white tip = newest 3 segments while Tracing. First point's 'previous' = the real body's
    previousPosition (ship: center - movementThisTick) so the line starts at the body as DRAWN (alpha interpolation).
    Look constants are static constexpr at the top of the trajectory section in RenderSystem.cpp.
  - Readout: TransferScenes FPS label + planningReadout() label in a UIColumn top-left ('PLANNING: press I to
    integrate' / 'PLANNING: Integrating...' / 'PLANNING').
  - Verified: trajectory_harness (9005 points vs the real game: 0 differ; behaviour checks all as expected),
    universe_refactor IDENTICAL, CMake Release build 0 warnings. NOT yet seen on screen by anyone: Marco to try it.
  - Round 2 (2026-11-24, after Marco tried it: 'solid'): static bodies' lines AMBER (TRAJECTORY_STATIC_COLOR);
    Marco REJECTED drawing relative to a reference body (keep world space). PLACEMENT PREVIEW: while the left button is
    held on the world in planning mode, updateSpawnPreviewTrajectory predicts the planet being placed
    (makeSpawnPreviewBody in PhysicsSystem.cpp = populateGravBodyProperties... with the preview velocity flag, so it IS
    the spawned body); preview universe = planets + ship + that planet, NO debris; restarts when the planet changes;
    no speed limit, MAX_PREVIEW_MILLISECONDS_PER_FRAME 4; drawn white (amber if the spawn will be static), no fade, no
    tip; Trajectories::spawn_preview_path. TIME MARKS every TRAJECTORY_TIME_MARK_SECONDS (5 game s), 5 px across the
    line at any zoom, toggled by UIState::getShowTrajectoryTimeMarks (default on, no key yet: will be a setting).
    Horizon 60 -> 180 s, TRACE_SPEED 10 -> 30, tick cap 32 -> 64. MAX_TRAJECTORY_LINE_VERTICES 327680 with a
    static_assert in PhysicsConstants (worst case 27 lines) so tuning can't silently cut lines. Render rebuilds and
    re-uploads all line vertices every frame: worst case ~9 MB/frame (fine for now; append-only upload if it shows).
    updateTrajectories now also takes const CameraState&. Harness: 27005 trace points + 5401 preview points vs the
    real game, 0 differ.
  - Round 3 (2026-11-24): TRACE_SPEED 30 -> 15 (Marco). APPEND-ONLY UPLOAD: new TrajectoryLine.vert/.frag shaders (added
    to BOTH MakeTransfer scripts; the .msl/.spv only exist after MakeTransfer) compute the fade-in from a per-vertex
    drawn_at + a clock uniform (vertex uniform slot 1, TrajectoryShaderConstants must match the HLSL cbuffer), so a
    traced segment never changes once uploaded. TrajectoryLineVertex {x,y,prevX,prevY,r,g,b,a,drawn_at,fades_in} in
    GPUTypes.hpp. Traced segments: m_trajectory_line_vertex_buffer, each frame uploads only new segments at offset
    uploaded_count with cycle=false (safe: only writes past what earlier frames draw); a new trace (Trajectories::
    trace_id changed) restarts at offset 0 with cycle=true. Per-frame buffer (preview, time marks (zoom-sized), tips)
    rebuilt + uploaded whole. MAX_TRAJECTORY_LINE_VERTICES 294912, MAX_TRAJECTORY_FRAME_VERTICES 16384, both
    static_asserted in PhysicsConstants. Debug lines still use DebugLine shaders. Harness still 0 differ; the GPU
    upload path itself can only be checked on screen.
- **F4 PROFILER (2026-11-24, implemented by Claude at Marco's request, uncommitted):** engine class
  DynamoEngine/Debug/FrameProfiler (addSection -> id, record(id, ms) adds up within a frame, endFrame; every
  PROFILER_WINDOW_FRAMES (300, GameSystemConstants) builds a report: average, '1% low' = average of the slowest 1% of
  frames (at least the slowest one), worst; takeNewReport() true once; millisecondsBetween(perf counters)). 9 googletests
  in Tests/DynamoEngine/Debug/Test_FrameProfiler.cpp (148/148 pass, TidyEngine 40 files clean). Game owns m_profiler +
  ProfilerSections ids: Input, Spawns, Audio, Physics (all ticks summed), Trajectories, Render (CPU side only: the GPU
  runs async), Frame work (before the limiter), Frame total. Game::publishProfilerReport formats lines (monospace
  columns) into UIState::setProfilerLines and printf's them while the overlay is on. F4 toggles
  UIState::getShowProfiler; TransferScenes adds MAX_PROFILER_LINES (12) labels under the FPS/planning readouts.
  The old unused input_time/... locals and commented printf in Game::Run are gone.
- **GHOST CHECKBOX (2026-11-24, Claude):** SpawnSettings::is_ghost -> isCollidable = false (populate... for planets,
  the cluster loop for every particle). Third checkbox in 'Global Instantiation Properties' (that section is now 624
  wide; panel ~1076 of 1280). Ghost = passes through everything incl. the ship, keeps gravity. Built on isCollidable
  because isMacroGhost is NOT checked for particles in handleMacroParticleCollisions (only for the macro side): a
  'ghost to planets only' variant would need that check first. Harness ghost_spawn_. Marco asked whether debris
  alters trajectories: yes for planets in the main trace (it pulls them; no collisions in the prediction), never the
  ship directly (ship gravity loops planets only), not at all in the placement preview.
- **PLANET GHOST (macro ghost) FINISHED (2026-11-24, Marco: 'finish it'):** isMacroGhost was never set before (dead
  since 14827cf). Now SpawnSettings::is_macro_ghost -> isMacroGhost (planets via populate..., every cluster particle).
  Meaning: passes through planets AND the ship (the ship counts as a big body; planets already skipped macro-ghost
  planets vs the ship), still collides with debris; gravity unchanged. Added the missing particle-side check in
  handleMacroParticleCollisions (outer loop skip) and in handleShipCollisions' particle loop. Spawn panel: 'Global'
  section is now TWO rows (Accretable, Bounce / Ghost, Planet Ghost): 4 in a row + the macro section was ~1288 > 1280.
  makeSpawnSection takes unique_ptr<UIElement> (a row or a column of rows); SPAWN_PANEL_CHECKBOX_ROWS = 2 feeds
  SPAWN_PANEL_HEIGHT, so the slider row moved up with it. Both ghost kinds render at 70% alpha (shader, existing).
  Harness ghost_spawn_ covers ghost / planet ghost vs planet, debris and ship.
- **PREVIEW NOW INCLUDES DEBRIS (2026-11-24):** Marco noticed the drag preview ignored particle clouds (it copied only
  planets + ship). Now m_preview_universe_state = the whole universe + the planet being placed: the preview equals the
  I trace of that planet point for point (harness case 5: 0 / 5401 differ; a 2e8 cloud moves its end by 27167 px).
  COST (measured, 4 ms budget): 5 planets 0 particles = full 180 s instantly; 10 + 2400 = 34 s of path per second held
  still, 0.6 s while the mouse moves; 25 + 20000 = 2 s per second held still, ~0 while moving (each move restarts it).
  Options offered to Marco: (a) aggregate debris into <= N clumps above a threshold (exact below), (b) keep showing
  the previous preview while the new one regrows, (c) thread the particle x planet gravity (also speeds up the real
  sim; fits his CPU-threads-then-GPU particle goal).
  -> Marco's DECISION (same day): keep the preview CHEAP, NO debris (reverted to planets + ship + the new planet).
  Harness case 5 now documents it: the cloud moves the preview by 0 px, the I trace by ~27000 px.
- **PREDICTION TIMERS (2026-11-24):** Trajectories::predicted_seconds (= trace ticks done x PHYSICS_TIME_STEP, set in
  advanceIntegratedTrajectories, 0 on start/clear) and preview_predicted_seconds (placement preview, 0 when cleared).
  Readout: 'PLANNING: Integrating... 42.0 / 180 s', 'PLANNING: 180 s predicted' when finished, and a new line under it
  while placing a planet: 'PREVIEW: 12.3 / 180 s' (previewReadout in TransferScenes).
- HARNESSES cleaned 2026-10-12 (Marco asked): 12 kept + updated to the UniverseState API, all print the same as on
  d41ec3d (shatter_order varies run to run anyway: random debris); 7 obsolete files deleted; Notes/drafts/Harnesses/
  README.md = build recipe + what each checks. GOTCHA: .gitignore line 14 ignores Notes/drafts/ -> only 8 harness files
  were ever force-added; the rest (incl. universe_refactor_* and README) DON'T travel between machines. Asked Marco.
- Gotcha for the copy: Starship keeps m_thrust_acceleration from the last real tick (W held when P was pressed) -> the
  prediction must cut thrust (coast), or it thrusts forever.

**Next candidates** (let Marco pick):
- **CAMERA LIMIT CHECK RESOLVED + RENAMES (2026-11-25):** Marco's F3 readout showed the limit WORKS: at zoom 0.01 his
  world (MacBook desktop 1512 x 982 points / MIN_ZOOM = 151,200 x 98,200 around (640, 360)) leaves the view (1280 x
  720 window / 0.01 = 128,000 x 72,000) only x [-10,960, 12,240], y [-12,740, 13,460]; he was parked at (12,240,
  -12,740). Claude's earlier '~32.6k' assumed a 1920 x 1080 display: always compute from HIS screen. What he
  really expected: WORLD SIZE as the master setting (zoom-out limited to where the world fills the window, star field
  sized to cover the world) instead of MIN_ZOOM sizing everything. Agreed to fold that into the camera rework.
  Renames (Marco: 'clamp' terminology unclear): Camera/CameraClamp -> Camera/CameraLimits, clampCameraToWorld ->
  keepViewInsideWorld, clampViewToRectangle -> keepViewInsideRectangle, CameraState::debug_last_clamp_correction ->
  debug_view_push_back (still TEMPORARY with the F3 readout, which now says 'pushed back by'), harness
  camera_clamp_ -> camera_limits_. US ENGLISH from now on (memory transfer-us-english): swept colour/neighbouring/
  normalised from src, HLSL and harnesses (Marco had already done centre -> center).
- **REDUCED TEST WORLD REMOVED (2026-11-25, Marco):** it only moved the world-edge rectangle (camera limit + optional
  walls); physics has no size limit and the star field stays full size, so the collision box covers the use case.
  Gone: USE_SMALL_TEST_WORLD, TEST_WORLD_WIDTH/HEIGHT, the TestVisualScene world save/restore, the grey 'small world'
  outline. Kept: UniverseState world corners (= the star field's extent; the camera clamp uses them) and
  REFLECTIVE_WORLD_WALLS (test scene only, red outline at the star field's edge). Older notes above that mention the
  small test world describe the removed version.
  MARCO OBSERVED (2026-11-25): with the small test world on, the camera was NOT stopped, contrary to Claude's
  code-reading claim ('the view locks to the world's center'). Not investigated (feature removed). The same
  clampCameraToWorld still guards the star-field edge: if the camera ever pans past the star field, look there
  (and check which build/scene was running). Don't describe camera behaviour as fact without running it.
  FOLLOW-UP (same day): clamp moved to Camera/CameraClamp.hpp/.cpp (free function clampCameraToWorld(camera,
  world_min, world_max) returning the correction); camera_clamp_harness proves the MATH is right. TEMPORARY F3
  readout (TransferScenes addSandboxUI: CAMERA center/zoom/view + WORLD edges + 'clamp moved', from
  CameraState::debug_last_clamp_correction) to see it live. Prime suspect for what Marco saw: the star layer's own
  offset (twinkling_star_offset) is only clamped to the star field, so stars keep sliding while the world is locked.
  Remove the readout + debug_last_clamp_correction once confirmed.
  Marco is on a TRACKPAD (no middle button, and the engine only carries VERTICAL scroll = zoom), so he could not
  pan at all. Added: arrow keys pan (InputSystem::updateCamera(game_state, frame_seconds), CAMERA_KEY_PAN_SPEED 900
  UI points/s since 2026-11-26). Future (his call): trackpad-native camera = two-finger swipe pans, pinch zooms; needs horizontal
  wheel (+ pinch) support in SDLInputIntake/InputState and changes zoom for mouse users. Belongs with the camera work.
- **F5 PHYSICS PROFILER (2026-11-25):** engine DynamoEngine::ProfileScope (FrameProfiler.hpp: RAII timer, records into a
  section when its { } block ends; `enabled` false = records nothing; 3 googletests, 151/151, tidy 40 clean).
  PhysicsSystem owns m_profiler (PROFILER_WINDOW_FRAMES) + ProfilerSections, rows: Coll: planets / planet-debris /
  debris grid / debris pairs / ship, Verlet 1st half, Grav: planets / debris / ship, Verlet 2nd half, Walls & box,
  Cleanup, Whole tick, Prediction. One ProfileScope line at the top of each timed function, gated by
  m_profiling_real_tick (true only inside advanceUniverseState) so predictions running the same functions on copies
  don't pollute the rows; Prediction has its own always-on row (updateTrajectories). Per FRAME, summed over ticks.
  Game calls physicsSystem.endProfilerFrame() + publishPhysicsProfilerReport. UIState::ProfilerOverlay {Off, Frame,
  Physics}: F4/F5 switch between them (both wouldn't fit), same key hides. MAX_PROFILER_LINES 17. Shared formatting:
  profilerReportLines() in Game.cpp. First measurement (20 planets + 8000 particles, 2 ticks/frame, O2): whole tick
  3.36 ms/frame = debris pairs 1.45, debris gravity 0.78, planet-debris 0.56, grid build 0.49, the rest ~0.1.
- **STEP 2: PHYSICS ON ITS OWN THREAD (2026-11-27; NOT tested in-game yet).** Claude first built a per-frame
  physics time budget instead (and said why: a thread means the main thread stops sharing the universe and
  PhysicsSystem); Marco: 'why not do the refactor? If it goes wrong we just revert' -> built the thread, the budget
  code removed again. Core/PhysicsThread.hpp/.cpp: owns the live UniverseState + its own PhysicsSystem; run() loop:
  copy the FrameRequest (ship controls, speed, planning, runs_simulation) + pending Commands under m_mutex, apply
  commands (spawn/undo/clear via updateGravBodyInstantiations, or std::function universe edits), accumulate owed game
  time (real x speed, capped at MAX_PHYSICS_BACKLOG_TICKS = 8 (renamed from MAX_PHYSICS_STEPS_PER_FRAME) -> dropped =
  behind), planning: adjustThrottleWhilePlanning in real time, AT MOST ONE TICK PER PASS (new controls/commands read
  between every two ticks), publish, then wait_for on a condition variable until the next tick is due or m_has_news.
  Publish = copy live -> m_staging OUTSIDE the lock, then swap(m_staging, m_published) under it, only if the main thread
  took the previous one (<= 1 copy per frame; ~1 ms at 50k, optimize later with a render-only snapshot). Main thread:
  takeNewestUniverse swaps m_published <-> game_state's universe (3 buffers rotate, no copy under the lock).
  Handoff: ticks, dropped_time, effective speed (measured over 0.5 s windows; = the CURRENT slider unless that window
  dropped time, so no readout flash on slider moves), publish_time, F5 report. Game: sendRequestsToPhysics (scene edits
  from GameState::editUniverse queue, spawn commands, frame request) and takeUniverseFromPhysics; F4 rows 'To physics'
  / 'From physics' replace 'Spawns' / 'Physics' (main-thread hand-off only, ~0.1 ms: Marco asked why they didn't
  match F5's ~3 ms), plus 'Physics thread' = Handoff::busy_ms, the physics thread's work (ticks, commands, publish
  copy) per frame, NOT part of 'Frame work' (runs in parallel); F5 report comes from the handoff; ticks line 'behind N%'; alpha =
  time since publish x effective speed / dt (1 while planning). predictionPhysics = a second PhysicsSystem for the
  trajectory lines and placement preview on the main thread. TestVisual onEnter/onExit queue their box edits
  (GameState::editUniverse); getUniverseStateMutable on GameState is now the presented COPY (comment says so).
  Readout (Marco's wording, always shown): 'FPS: 60    Current Sim Speed: 1.6x' on the FPS line ('stopped (planning)'
  while planning); the slider is now 'Requested Sim Speed' and SNAPS to 0.1 steps shown with 1 decimal (new
  UISlider::setStep(step, decimal_places) + text(); snapToStep divides by steps-per-unit so 0.7 is the double nearest
  0.7; 4 tests, 173/173, TidyEngine clean; radius/mass sliders unchanged, still continuous). Verified with
  Notes/drafts/Harnesses/physics_thread_harness.cpp under ThreadSanitizer: 0 reports; spawn + undo + edit arrive;
  thread result after 302 ticks bit-identical to direct stepping. 169/169, determinism MATCH.
- **DETERMINISM STEP 1 DONE (2026-11-27, Marco: 'go with the original plan'; Notes/PHYSICS_REQUIREMENTS.md is now
  the reference for physics work).** Decisions that day: debris is gameplay, CPU physics identical on Mac+Windows (no
  GPU physics), fixed rubble per level, slow machine = slow motion (sim on its own thread next, so the speed slider
  survives), 50k on a 4-core laptop, extensible materials, Soft Step + coloring as the solver plan, rotating shaped
  bodies later as a separate type. Implemented: DynamoEngine/Math/DeterministicMath (deterministicSin/Cos: Cody-Waite
  3-part pi/2 reduction + Taylor to r^17/r^18; deterministicCbrt: frexp/ldexp + 6 Newton steps) and
  DeterministicRandom.hpp (mixBits = SplitMix64 finalizer, combineHash, doubleBits, hashToUnitInterval/Range);
  13 tests (169/169). UniverseState::random_seed ("TRANSFER") + tick_count (++ in advanceUniverseState).
  bodyRandomKey(universe, body) = seed, tick, macroIdentifier, position bits; substituteWithParticles(+FromImpact)
  take a random_key (fragment speed 0.8-1.1 and impact size 0.5-0.9 hashed per fragment number); directionFromTo hashes
  the two positions; FragmentLayout, Starship, stackedPairDirection, handleAccretion use the deterministic functions;
  randomDouble removed from CustomMathUtilities.hpp. CMakeLists: -ffp-contract=off / /fp:precise for all targets;
  harness README + run_particle_benchmark.sh got the flag. NEW: Core/DeterminismCheck + main.cpp
  `--determinism-check` (6 scenes x 1200 ticks, checksum every 300, OVERALL line, determinism_check.txt; ~1 s).
  Verified: run-to-run identical; Release build == -O2 harness build; without the flag every line differs; seed probe:
  same seed identical, other seed different (1600 fragments). Mac reference: Notes/determinism_check_mac.txt
  (OVERALL 02F18F7B1A87A525). NOT yet run on Windows. Scripts/Apple/CheckDeterminism.sh [nobuild|update] and
  Scripts/Windows/CheckDeterminism.bat [nobuild] build Release, run the check, diff/fc against the reference (update =
  Mac saves a new reference after an intended physics change); documented in Scripts/README.md. The .bat is untested
  (written on the Mac). Side effect: all physics bits changed (flag + hashed randomness),
  so universe_refactor_baseline.txt and contact_identity outputs need regenerating before they're used as guards again.
  GATHER SOLVER REMOVED the same day (Marco: dead code): constants, DebrisSolverSettings, neighbor lists, F5 rows
  'Coll: debris nbrs/solve', UniformParticleGrid::queryNeighborhood all gone (grid files back to f174b8c);
  determinism check still MATCH, 169/169. debris_solver_harness (args now: scene, particles, settle seconds) and
  conservation_probe_harness kept for grading the sequential solver. The findings stay in the entry below.
- **GATHER (JACOBI) DEBRIS SOLVER TRIED: COLLAPSES DEEP PILES (2026-11-27, Marco: 'go ahead and attempt', threads
  later). [REMOVED later that day; this entry keeps the findings.] Built, behind a flag, DEFAULT OFF (USE_GATHER_DEBRIS_SOLVER = false; sequential path bit-identical to
  f174b8c, contact_identity checksums equal). Code: PhysicsSystem::resolveDebrisContactsGather ->
  buildDebrisNeighborLists (flat start/neighbors arrays, skin DEBRIS_NEIGHBOR_SKIN 0.1, contact counts) ->
  DEBRIS_SOLVER_ITERATIONS x runDebrisGatherPass (snapshot pos/vel, each particle sums debrisContactPush from its
  neighbors, computed from the LOWER index's side so both sides match bit for bit, writes only itself);
  UniformParticleGrid::queryNeighborhood (full 3x3); stackedPairDirection = hash of the two indices (no RNG).
  Runtime knobs: physics.debrisSolverSettings() {use_gather, iterations, overlap_relaxation, average_bounces}.
  FINDINGS: momentum exact in every mode. Summing full bounces from all neighbors INJECTS energy (two clusters fired
  together: kinetic x30 in 1 s) -> bounces averaged by max(contact count) (DEBRIS_GATHER_AVERAGE_BOUNCES), overlap
  corrections summed: then energy/spin match sequential in free space (conservation_probe_harness: -15% KE, same
  +0.2% L drift from the position push). But piles under gravity COLLAPSE: wall 4800 (1e10 planet 1400 outside the
  box): sequential 6.97 touching/particle, 0.91 ms; gather x1 47.7 touching, 6.9 ms; x16 33.9, 15 ms; x8 with
  relaxation 1.8: 28.3. Wall 14,400: seq 14.4 touching 4.06 ms; gather 72-81 touching, 34-47 ms. Planet scene
  (spin_harness, a hot gas in both): gather keeps it hotter (845 vs 437 px/s) and spins (+15..+50 deg/s, from the
  larger position pushes: with overlap correction off spin -0.13 but it collapses). Causes: (1) a push travels one
  layer per pass, deep piles need ~layer-count passes; (2) the position push doesn't change velocity, so gravity's
  velocity keeps driving particles in, and averaged bounces remove only part of it; (3) the box/wall clamp isn't
  inside the pass loop. Pivot offered: keep Gauss-Seidel quality and parallelize by COLORING (stripes / 3x3 cells),
  plus Verlet-style position-based debris (velocity implied by the position change) + substeps, the usual recipe of
  100k-particle sims (Marco's video: 'AIs try to create a large C++ physics Simulation', content not readable).
  F5 now: rows that didn't run in the window are hidden (FrameProfiler SectionReport::frames_recorded, Game.cpp
  profilerReportLines; tests added, 156/156), Verlet halves merged into one row, new rows Coll: debris nbrs / solve,
  fitToOverlay shows '(+N rows: see console)' if a report ever overflows MAX_PROFILER_LINES (18).
  Harnesses: debris_solver_harness (planet|wall scenes, seq|gather, prints ms/overlap/touching/speed/spin),
  conservation_probe_harness (momentum/spin/KE in free space).
- **OPPOSITE-SIGN MERGES: 10x RULE EVERYWHERE (2026-11-26, Marco's design call):** asked 'what is broken with
  negative accretion, shouldn't it just subtract?'. Findings: subtraction itself was fine (pickMassPair compares abs, so
  the radius pow never sees a negative ratio: the old REWORK note was wrong); only a near-cancelling merge was a problem
  (velocity = momentum / ~0). Claude first made opposite signs merge at ANY ratio + annihilate under 1% left; Marco:
  'gentle contact should bounce' -> '10x rule everywhere; small debris should be able to eat them down; I care about the
  magnitude of mass'. FINAL: no special case. The ratio already uses |mass|, so a merge needs a 10x size ratio
  whatever the signs, and always leaves >= 90% of the heavier mass: new_mass can't get near 0 or flip sign (comment in
  handleAccretion). Small debris eats a planet down to ~10x the eater's size, then bounces; a planet never reaches
  zero. Marco LIKES the momentum fling of negative merges ('rail gun effect is a cool mechanic'): keep the momentum
  rule. REWORK handleAccretion negative-mass item removed. negative_accretion_harness: -1e5 hits take 1e7 down to 9e5
  (planet speed builds 10 -> 858 px/s from the kicks), hit 92 bounces; -1e6: one merge then bounce; -1e7: bounce.
- **NEGATIVE-MASS BOUNCE FIX + NaN SAFETY NET (2026-11-26, Marco: 'add both'):** Marco: a settled pile at 60 FPS drops
  to 2-3 FPS when a NEGATIVE cluster is added (not a particle-count issue). Cause: handleElasticCollisions' dynamic
  branch divided by m_smaller + m_larger (0 for -m hitting +m: instant inf) and, for unequal opposite-sign masses,
  multiplied the speed every contact (-2m on +m: ~3x per bounce; ~1e59 within 0.5 s). inf -> NaN, NaN spreads by
  contact; every NaN position lands in grid cell (0,0) and NaN fails getCollisionInfo's '>' early reject, so all NaN
  pairs count as contacts: one ~1600-particle cell, n^2 (5 ms/tick at 5600 particles, worse in his bigger pile).
  FIX: m_smaller/m_larger = std::abs(mass) in the bounce (negative mass bounces like normal matter of that size;
  gravity keeps the sign). SAFETY NET: cleanupGravBodies first calls markNonFiniteBodiesForDeletion (static) on both
  lists, deleting bodies with non-finite position/velocity, and prints a std::cerr line on real ticks. The net alone
  does NOT save the runaway case (speeds reach ~1e300 while still finite, the box clamps them into corners and they
  pile up there), so the abs fix is the real fix. Checked: negative_cluster_harness (new, in Harnesses) -1e6 / -2e6 /
  -5e5: 0 NaN, fullest cell 4-7, ~0.55 ms/tick; placement_separation identical except the +1e6/-1e6 line (50.0 ->
  50.2 px, it now bounces); positive masses untouched (abs(x) = x). NOTE: at this point the working tree no longer
  has separateOverlappingDebris / DEBRIS_SEPARATION_FRACTION (the entry below) and particle gravity is commented out
  again: Marco rolled those back himself. Negative-mass ACCRETION: settled right after (entry above).
- **DEBRIS SEPARATES WITH VELOCITY (option 1, 2026-11-26, Marco chose it over PBD 'for now'):**
  handleElasticCollisions: when both bodies are particles, the overlap is cleared by separateOverlappingDebris (static
  in PhysicsSystem.cpp): equal and opposite velocity along the normal, shared by inverse mass, only TOPPING UP to the
  separating speed that clears DEBRIS_SEPARATION_FRACTION (0.2) of the overlap next tick. Planets keep the position
  push (OVERLAP_CORRECTION_FRACTION 0.8, OVERLAP_SLOP 0.01, now file constants). First try used 0.8 for debris too:
  each cluster squeezed by its own gravity BURST apart at ~140 px/s around 0.6 s (top-ups add up along contact chains
  and the speed stays); 0.2 and 0.1 both fine. gravity_spin_harness (particle gravity on): spin 0.0000 for 60 s,
  energy falls then levels off (-2.9e10), no climb, no burst. Piles unchanged: 4800 0.90 ms, 14,400 4.28 ms (was
  4.19), same packing, still settle at 1400; the 590 case still doesn't settle (as with one pass before).
  spin_harness (pile on a STATIC planet, no particle gravity) still drifts 0.004-0.008 deg/s: the static-planet
  contact (one-static branch) still moves the debris radially (position push), which doesn't conserve spin about
  the planet; planets could get the same velocity treatment later. Harness changes (all checks pass): stacked_spawn
  worst overlap 0.11 -> 0.00 px; spawn_flags '1 gentle hit, bounce 0' leaves 6 crumble fragments unabsorbed at the
  end (was 0; planets 1 as expected); static_spawn '1' mass 1001000 -> 1000945 (a few fragments not absorbed yet);
  ship_collision ship v changed again; universe_refactor baseline REGENERATED (third time today).
- **RUNAWAY SPIN CAUSE FOUND (2026-11-26): the debris PUSH-APART.** gravity_spin_harness.cpp (PROBE harness: needs
  a scratch copy of src with particle-particle gravity uncommented and g_probe_push_apart / g_probe_alternate_sweep
  globals; header says how). Two 800-particle clusters at rest in empty space, own gravity only, total spin starts at
  exactly 0. Push-apart OFF: spin stays 0.0000 deg/s for 60 s (both sweep modes). Push-apart ON: -0.37 deg/s
  (alternating) / +0.07..0.19 (one-way) from nothing. Why: handleElasticCollisions' penetration correction moves
  positions without touching velocities; it is mass-weighted along the line of centers, so the center of mass and
  momentum are kept (momentum stays ~1e-8), but L = sum m r x v changes by c n x (v2 - v1) for sliding contacts.
  It also ADDS ENERGY: with push-apart ON total energy falls until ~30 s, then CLIMBS (-3.03e10 -> -2.73e10 by 60 s,
  kinetic x4.8): lifting overlapping particles apart against their own gravity for free. Self-gravity keeps clumps
  overlapping, so the push runs every tick -> runaway. (Push-apart OFF: energy only falls; the cloud's rms radius
  still grows to ~1000, which looks like a dense core flinging out a halo, not energy gain; not verified.)
  Fix options offered (physics design, Marco's call): separate with a velocity impulse along the normal instead of
  moving positions (Baumgarte-style: conserves momentum AND spin exactly), or position-based dynamics (velocity
  derived from the corrected positions).
- **ALTERNATING DEBRIS SWEEP (2026-11-26, Marco; hoped it fixes runaway spin he saw with particle-particle gravity
  re-enabled):** handleParticleParticleCollisions flips m_sweep_debris_backward every call (once per real tick) and
  walks the particles forward / backward; resolveContactsOfParticle(i, backward, ...) (split out) also reverses the
  candidates, so the whole pair order reverses. Piles: packing back to pre-sort (14,400 pile 19.9 candidates / 7.06
  touching), tick 4.69 (before sorting) -> 4.19; 4800 pile 1.00 -> 0.90. Open box 80k: ~7.2-7.4 spread / ~7.9
  packed (one-way sort was 6.88 / 7.49: ~5% slower there, measured 3 runs). Spin: new spin_harness.cpp (static planet
  at the origin, no box, 4800 bounce debris, 60 s; nothing in the scene can torque the pile): one-way 0.006-0.011
  deg/s, alternating 0.003-0.006 deg/s. Halved, but neither is a runaway (about one turn per 10+ hours), so the
  runaway spin with particle-particle gravity likely has another cause; not reproduced yet (that code is commented
  out in addGravityBetweenGravBodies). Candidate cause to check: push-apart moves positions without changing
  velocities, which does not conserve angular momentum (sum of m * dr x v != 0).
- **DEBRIS SORTED BY CELL (2026-11-26, Marco: 'get it ready for the GPU'):** UniformParticleGrid::sortByCell =
  chooseCellSize, each particle's cell, key = cell number counted column by column inside the occupied box
  (column-major), LSD RADIX SORT (counting sort per byte, only as many passes as the largest key has bytes), gather
  into m_reordered_particles + swap. Called at the very START of advanceUniverseState (before
  rememberCollisionBoxSides, which keeps particle indices until after the collisions), every
  DEBRIS_SORT_INTERVAL_TICKS = 8 (PhysicsConstants; 1 = every tick, what a GPU version would do). New F5 row 'Sort
  debris by cell' (17 overlay lines now: MAX_PROFILER_LINES is exactly full). FIRST TRY sorted in the grid's SLOT
  order: no good (hash scatters neighbor cells; 80k tick got slower, 7.59 -> 8.47); column-major fixed it.
  Results, ms/tick (before sort -> sorted every tick -> every 8 ticks): 80k spread 7.59 -> 7.46 -> 6.88; 80k packed
  8.78 -> 8.21 -> 7.49; 8k spread 0.77 -> 0.76 -> 0.735. Settled PILES: no net gain: per contact ~12% faster, but
  the new (spatial sweep) resolution order packs piles ~12% tighter (14,400 pile: 20.3 -> 22.4 candidates, 7.06 ->
  7.89 touching per particle; tick 4.69 -> 4.64; 4800 pile 1.00 -> 1.06). Idea offered: alternate the sweep
  direction each tick (Gauss-Seidel sweep bias). Harnesses: all checks pass; ship_collision '3 debris' ship velocity
  changed (debris met in a different order; momentum conserved); universe_refactor baseline REGENERATED again.
  contact_identity_harness's reference must be re-recorded from the current code before the next speed-up.
- **CONTACT MATH SPEED-UP 'A' (2026-11-26, bit-identical):** getCollisionInfo rejects with squared distances first
  (margin 1e-9 so edge pairs still get the exact sqrt test: same answer always) and records bodyA (CollisionInfo
  ::bodyA); handleDynamicCollision's 4 bounces go through bounceUsingCollisionInfo, which reuses getCollisionInfo's
  distance + normal (negated when body A is the heavier; stacked pairs still get directionFromTo's random
  direction); handleElasticCollisions(smaller, larger, distance, smaller_to_larger) no longer recomputes them
  (was 3 sqrt per contact, now 1). NOT done (would change bits): 1/(m1+m2) precompute; the candidate list visitor
  (needs a template callback, not worth the readability). Verified: new contact_identity_harness.cpp (3 scenes,
  checksums of every body's exact bits every 200 ticks: identical before/after), 154/154, universe_refactor
  baseline identical, all deterministic harnesses identical. Gain smaller than guessed: settled 4800 pile pairs
  1.024 -> 0.913 ms (-11%); open box 80k spread tick 9.45 -> 7.59, 80k packed 11.66 -> 8.78, 8k spread 0.86 -> 0.77.
  Finding: a settled 14,400 pile = 4.69 ms/tick with 20 candidates and 7 TOUCHING per particle: more than the 6
  neighbors equal circles can touch without overlapping, so deep piles are compressed (one 80% push-apart pass per
  tick can't hold the weight) and contact count, not per-contact math, is the cost. Box-corner quirk: particles
  clamped into a corner end up exactly stacked (e.g. (997.45, 997.45)) and take random push directions.
- **SINGLE DEBRIS PASS (2026-11-26, Marco's call + his edit):** handleCollisionEvents = planets, planet-debris,
  SHIP, then ONE debris-debris pass (the second pass after the ship is gone). Comments fixed (ProfilerSections, the
  benchmark text). 154/154 tests; harness checks all still pass with tiny number changes (ghost, ship_collision,
  stacked_spawn: still 0 stacked pairs); universe_refactor baseline REGENERATED (header + README note).
  Timing: 8k spread 1.31 -> 0.86 ms/tick, 80k spread 13.3 -> 9.5, 80k packed 15.4 -> 11.7; settled 4800 pile 2.2 ->
  ~1.1-1.2 ms. SIDE EFFECT FOUND: the second pass was also a second round of separating every pile (it bled energy
  out). With one pass, 4800 bounce debris + 1e10 planet 590 from the wall NO LONGER SETTLES (90 s: median 3.7
  diameters/tick, 955 at the wall); 900 and 1400 still settle (median 0.66 / 0.34). The settling threshold moved
  from 234..590 to 590..900. Reported to Marco; fix offered: debris bounces that lose a little energy.
- **DENSE PILE + 'NEVER SETTLES' REPRODUCED (2026-11-26):** Marco's screenshots: 4800 bounce debris (radius-96
  clusters) + 1e10 static planet left of the TestVisual box. Piled: 2.98 ms/tick (pairs 2.86), 41 FPS at 3x (8
  ticks/frame AT CAP). Not piled: 0.9 ms/tick, stays a gas for minutes. testvisual_load_harness args 5/6 (planet
  center's distance from the wall, bounce debris) + end-of-run 'movement per tick in diameters'. Planet 234 from the
  wall (his 'not settling' shot): pull at the wall 1e10/234^2 = 1.8e5 px/s^2 (1500 px/s gained per tick), median
  particle moves 10.6 DIAMETERS per tick (90%: 19.5) -> collisions checked once per tick can't catch them, they pass
  through each other, elastic walls (WALL_RESTITUTION 1) + elastic bounces keep the energy: a gas that never cools.
  Planet 590 from the wall: piles after ~25 s, then 11.3 candidates / 3.9 touching per particle, pairs 2.1 ms, tick
  2.2 ms, median particle still moves 0.95 diameters per tick (a 'boiling' pile). Not a grid bug: time step vs
  gravity strength. Options offered (Marco's call): substeps when particles move > ~0.5 diameter/tick, restitution
  < 1 for debris, speed cap. Dense pile cost is contact count, not grid: second pass halves it; then threads.
- **2 MS/TICK TARGET (2026-11-26, Marco):** goal = worst-case tick <= 2 ms (or less), 'spectacular' particle counts.
  He sees 60 FPS up to 6 ticks/frame (~3 ms/tick at ~19k), drops past that; worst case he found: 10-12k free debris
  in the TestVisual box + one static planet outside -> dense pile, big slowdown. testvisual_load_harness extended
  (args: cluster radius, max particles, planet mass, seconds; prints cand/p, touch/p, worst tick, end-of-run wall
  pile counts). NOT reproduced yet: 11.2k particles (radius-100 clusters) + 1e10 static planet outside the right
  wall = 1.03 ms/tick, cand/p 0.4, debris NOT piled (mean x 106; 887 within 100 px of the right wall) -- WALL_RESTITUTION
  1.0 and elastic particle bounces keep it a hot gas. Need Marco's exact setup (cluster radius, mass, planet mass/
  position, spawn flags, how long he waits). Plan after that: worst-case benchmark scene with pass/fail vs 2 ms,
  then second pair pass decision, sort particles by cell, pair-loop cleanups, threads.
- **F5 PER TICK + TICKS/FRAME LINE (2026-11-26):** Marco saw frame time jump 16 -> 45 ms going from 11.2k to 12k
  particles in TestVisual and suspected catch-up ticks hid the grid gains. Checked with
  Notes/drafts/Harnesses/testvisual_load_harness.cpp (real spawn path, radius-100 clusters into the 2000 x 2000 box,
  built against HEAD's old grid and the new one): per tick at 20k particles old 5.2-5.7 ms, new 2.5-3.1 ms, so the
  grid is ~2x faster in his scenario too; the jump is frame-level. Changes: advanceUniverseState now calls
  m_profiler.endFrame() itself (endProfilerFrame() removed), so F5 = ms PER TICK ('Physics/tick (ms)'); F5's
  Prediction row removed (F4's Trajectories row times the same work once per frame). Game counts ticks per frame
  over the F4 window: both overlays show 'Ticks/frame: X avg, Y max (cap 8)' (+ 'AT CAP'). Mechanism to explain to
  Marco: 60 FPS = 2 ticks/frame at 1x; a frame over 16.7 ms owes more game time -> more ticks next frame -> longer
  frame; settles if a tick costs < 8.33 ms / speed, else runs at MAX_PHYSICS_STEPS_PER_FRAME 8 every frame; vsync
  (SDL_GPU default) snaps frame times to refresh steps, making it jumpier. Open: what pushed the frame past 16.7 ms
  at 12k (needs his F4 readout); overload policy (slow motion = lower cap and drop the backlog) offered.
- **HASHED PARTICLE GRID (2026-11-26, step 1 of the particle speed-up, Marco: 'go ahead and implement'):**
  UniformParticleGrid rewritten: cell size = 2 x largest radius as before; each cell hashed (two odd 64-bit
  multipliers, XOR, fold high half) into a power-of-two table of >= 2 slots per particle; counting sort (count ->
  running total -> place, stable in particle index order) fills m_slot_start + m_sorted_entries {column, row,
  particle_index}; the query walks the same forward half-stencil (own cell with j > i, (0,+1), (+1,-1..+1)) and skips
  entries from other cells sharing the slot. No std::sort, no binary search. Counters (Marco's sudo run) had shown
  the old grid build was mostly std::sort branch mispredicts (~6-7 per particle per sort) and the query ~340
  instructions/particle of binary search. Verified: Notes/drafts/Harnesses/particle_grid_harness.cpp (vs a copy of
  the old grid, OldUniformParticleGrid.hpp/.cpp, and brute force; 8 scenes incl. slot sharing and clusters 1e9 apart:
  identical pair sets, no duplicates, nothing missed); 154/154 googletests; every deterministic harness identical,
  universe_refactor baseline still EXACT (pair order changed but none of those scenes depends on it). Slot sharing
  measured ~= random hashing: 15-20% of occupied cells share a slot. Benchmark, ms per tick before -> after:
  8k spread 2.19 -> 1.31 (grid 0.54 -> 0.08, pairs 1.34 -> 0.94); 80k spread 34.0 -> 13.3 (grid 6.9 -> 0.96,
  pairs 23.8 -> 9.6); 80k packed 36.1 -> 15.4. Sorting particles in memory by cell would still give 1.3-1.7x on pairs
  (next step, not done). Still open: the second pair pass per tick (Marco's call), squared-distance early-out.
  Side fix: the old grid's packed keys (CELL_COORD_BIAS 1e6 cells) broke past ~6 million px from the origin; the
  hashed grid has no such limit.
- **PARTICLE BENCHMARK HARNESS (2026-11-26, start of the 10x particle work):** Notes/drafts/Harnesses/
  particle_benchmark_harness.cpp + KperfCounters.hpp (kperf private framework via dlopen; counters need sudo; the
  counter path was NOT run by Claude, who can't sudo) + run_particle_benchmark.sh (builds -O2 into bin/, runs with
  sudo unless --no-sudo). Scenes: 8k/80k bounce-only particles (radius 3) in the collision box, spacing 3.0 and 2.2
  radii, 8 light static planets outside. First numbers (no counters, M3 Pro), ms per tick:
  8k spread 2.15 (grid 0.54, pairs 1.31); 80k spread 33.1 (grid 6.9, pairs 22.9); 80k packed 37.6 (grid 7.1, pairs
  27.3). Findings: (1) handleCollisionEvents runs handleParticleParticleCollisions TWICE per tick (before and after
  the ship), so grid + pairs are paid twice. (2) Same density, 10x particles: pair cost per particle 75 -> 144 ns
  (spread) = memory effect (1.1 MB of particles fits L2-ish comfortably, 11 MB barely). (3) Sorting particles by
  cell: pairs 1.57-1.73x faster at 80k, 1.16-1.25x at 8k; grid build 5-7x faster (std::sort on nearly-sorted input).
  (4) Only 6-7% of candidate pairs touch; ~65-90 ns per particle even sorted at 8k for ~1.6-3.3 candidates suggests
  the per-particle query (two binary searches) costs more than the narrow phase -> counters will confirm.
  Marco installed Xcode.app on 2026-11-26 (xcode-select now points at it); its license must be accepted
  (sudo xcodebuild -license) before clang works again; Instruments/xctrace come with it.
- **UI FITS ANY WINDOW SHAPE + PROFILER MOVED (2026-11-26, Marco: F5 overlapped the sliders; Mac fullscreen clipped
  the slider row's edge; 'anywhere absolute px values are used ... should be modified'):** cause: UIRoot scaled by
  window HEIGHT only, so Mac fullscreen 1512 x 982 gave a UI space of ~1108 x 720 and the 1150-wide slider row didn't
  fit. FIX (engine): UIRoot::fitScale = min(w / 1280, h / 720) (UI_REFERENCE_WIDTH added to GlobalConstants); 3 new
  googletests (TallerWindowFitsTheWidth, WiderWindowFitsTheHeight, DesignSizedElementFitsInAnyWindowShape), 154/154,
  tidy 40 clean. Profiler overlay (F4/F5) moved out of the top-left readouts into a TOP-RIGHT column: visor dropdown,
  8-point gap, then MAX_PROFILER_LINES labels PROFILER_TEXT_WIDTH 520 wide (46-char lines of 18 pt Space Mono, which
  is monospaced at 0.612 em = 11 points/char). A static_assert in TransferScenes.cpp checks the column ends above the
  slider row at 1280 x 720 (513 < 543), the tightest the UI space ever gets. World overlays that were fixed SCREEN px
  now grow with the window like the UI (sizes in 'UI points'): spawn velocity arrow (THICKNESS/HEAD_*_POINTS) and
  trajectory time marks (RenderSystem.cpp worldUnitsPerUIPoint(camera) = fitScale / zoom), and arrow-key pan speed
  (CAMERA_KEY_PAN_SPEED now UI points/s). Deliberately NOT scaled: how much WORLD a bigger window shows (normal for a
  world view; part of the camera/world-size rework), star sizes (world space), 1 px debug lines. Known small
  overlap: with F3 AND F5 on, the temporary 'WORLD: ... pushed back by' line runs ~3 characters under the profiler
  column; goes away when that temporary readout is removed.
- **Small-world constants clarification (Marco, 2026-11-25):** USE_SMALL_TEST_WORLD / REFLECTIVE_WORLD_WALLS apply ONLY in
  the TestVisual scene since the collision-box change (the game scene always uses the star-field world). Without
  REFLECTIVE_WORLD_WALLS the small world only limits the camera; bodies can still leave it.
  Marco then found TEST_WORLD_WIDTH/HEIGHT 'did nothing': applied, but invisible (no outline without walls) and at
  zoom 0.35 the window shows ~3657 x 2057 world px, so any smaller world just locks the camera to its center. FIX:
  buildDebugLines outlines the world edge whenever reflective (red) OR smaller than the star field (grey, derived
  from the sizes, no flag).
- **F4 PROFILER, FASTER (2026-11-25):** PROFILER_WINDOW_FRAMES 300 -> 60 (refreshes every second); the first overlay
  line is 'Particles: N   Planets: M', sampled when each report is made. With < 100 frames per window the '1% low'
  equals the slowest frame (= 'worst'); if both fast refresh AND a meaningful 1% low are wanted, FrameProfiler needs a
  rolling window (report every N frames over the last M) -- offered, not done.
- **TEST TOOLS: WORLD EDGES + COLLISION BOX (2026-11-25, Marco's request, to push collision density to a known
  upper bound before the 10x particle work):** UniverseState::world_min_corner/world_max_corner + reflective_walls
  (Game::setUpWorldEdges = the star field's extent, walls off; the camera clamp (InputSystem clampCameraToWorld)
  uses them; the star layer still clamps to the star field) and has_collision_box + collision_box_min/max_corner.
  Walls run at the END of integrateUniversePhysics (bounceOffWorldWalls, bounceOffCollisionBox), so the I trace and
  the placement preview bounce too (the preview copies walls + box: it builds its own universe). WALL_RESTITUTION 1.
  The box is TWO-SIDED: a body whose center was inside stays inside, one from outside bounces off the outside.
  GOTCHA found by the harness: handleElasticCollisions' push-apart sets previousPosition = position, so the box
  can't judge sides by previousPosition after collisions: advanceUniverseState remembers every position BEFORE
  handleCollisionEvents (rememberCollisionBoxSides) and runs bounceOffCollisionBox(BoxSides::BeforeCollisions) right
  after; inside the integrator it uses BoxSides::StartOfStep (previousPosition from the Verlet first half).
  Before the fix a crammed 800-particle cluster leaked up to 141 particles; after: 0. Drawn as lines (always, not
  only F3): world walls red, box cyan. TestVisual is now the TEST SANDBOX: TestVisualScene (TransferScenes.cpp) has
  the game's UI (addSandboxUI, shared with the game scene), a 'TEST SCENE' label, and onEnter puts a
  TEST_COLLISION_BOX_WIDTH x HEIGHT box (2000 x 2000) around the current view and, if USE_SMALL_TEST_WORLD, shrinks
  the world around it (+ REFLECTIVE_WORLD_WALLS); onExit removes the box and restores the saved world. Start menu has
  a 'Test Scene' button; Esc/Resume return to the scene that was paused (UIState::getResumeScene). Known: the two
  scenes have separate slider widgets (knobs can show stale values if both are used; they can't be switched between
  without a restart today). Harness world_walls_. Marco's constant name RUN_GAME_WITH_PLAYER_SHIP has a typo (SHP).
- **DEPRECATED_InputState RETIRED (2026-11-25, branch FixInputState, implemented by Claude, Marco reviewing):**
  24 of its 41 fields were dead. Replaced by Core/PlayerInput.hpp: ShipControls {thrust, rotation: +1/-1/0},
  SpawnProperties {is_cluster, with_velocity, world_position, world_velocity}, PlayerCommands {spawn + spawn_aim,
  undo_newest_spawn, clear_all_bodies, toggle_trajectories} (InputSystem clears them EVERY frame: no more
  resetTransientFlags), PlayerInput {is_aiming_spawn, spawn_aim, commands, ship}. Game owns m_player_input (Marco's
  choice) and passes it: inputSystem.processSystemInputFrame(..., m_player_input, ...),
  physicsSystem.updateGravBodyInstantiations(universe, player_input, spawn_settings), advanceUniverseState(universe,
  m_player_input.ship), updateTrajectories(frame_seconds, universe, player_input, trajectories, ui_state),
  renderFullFrame(game_state, ui_state, player_input, scene). The spawn position/velocity is computed ONCE
  (InputSystem::aimSpawn); physics spawns no longer need the camera; populateSpawnedBodyProperties replaces
  populateGravBodyPropertiesFromDEPRECATED_InputState. Sliders: SpawnSettings::radius/mass, UIState::setSimulationSpeed.
  Starship::applyThrust/applyRotation take const ShipControls& (Starship no longer includes UIState).
  VERIFIED as a pure refactor first: all 12 deterministic harnesses byte-identical to before (random ones: same outcomes).
  One deliberate difference: A+D held together now cancel exactly (before: -0.05 then +0.05 rad).
  Harnesses use Notes/drafts/Harnesses/HarnessInput.hpp (releaseSpawn, pressDelete, aimAtScreen, shipControls).
  GOTCHA seen: Marco's editor formatter reflowed long trailing comments in PlayerCommands INTO the code (fields
  vanished into comments). Comments on their own line above a field are safe.
- **SPAWN PANEL CLEANUP (2026-11-25, Marco):** Planet Ghost checkbox REMOVED (SpawnSettings::is_macro_ghost gone;
  physics still honors GravitationalBody::isMacroGhost, nothing sets it; the ghost_spawn_ harness sets it on bodies).
  Exclusion rules in ONE place: enum class SpawnFlag {Static, Shatterable, Accretable, Bounce, Ghost} +
  SpawnSettings::setFlag(flag, on) / flag(flag) (Core/SpawnSettings.cpp). Turning ON: Ghost -> off Bounce, Accretable,
  Shatterable; Bounce -> off Ghost, Accretable, Shatterable; Accretable/Shatterable -> off Ghost, Bounce; Static
  works with everything. Turning OFF never turns anything on. addSpawnCheckbox(row, label, settings, flag) binds
  through it. Global section back to ONE row (Accretable, Bounce, Ghost; SPAWN_PANEL_CHECKBOX_ROWS = 1; the column-of-
  rows support in makeSpawnSection stays). Harness spawn_settings_ checks the table. Marco renamed SpawnAim ->
  SpawnProperties (PlayerInput.hpp) in his pass; HarnessInput.hpp follows.
- **SHIP OPTIONAL (2026-11-25, Marco: off by default while he irons out other problems):** UniverseState::has_ship
  (default TRUE in the struct, so harnesses keep their ship); Game::StartGame sets it from RUN_GAME_WITH_PLAYER_SHIP
  (GameSystemConstants, currently FALSE). Without the ship nothing touches it: handleShipCollisions,
  addGravityBetweenShipAndBodies, applyShipInputPhysics, adjustThrottleWhilePlanning, both Verlet halves, the ship
  trajectory line, uploadStarship, the F3 hitbox and the THROTTLE readout all check has_ship. The placement preview
  copies has_ship. To bring the ship back: set RUN_GAME_WITH_PLAYER_SHIP = true (a key/setting later if wanted).
- **VELOCITY ARROW SIZE FIXED (2026-11-25, Marco's note):** buildVelocityVectorGeometry(start, end, zoom) sizes the
  shaft (12 px), head length (30 px) and head width (36 px) in SCREEN pixels (divided by zoom; the shader multiplies by
  zoom), and shrinks them proportionally when the drag is shorter than the head. appendPreviewBodies gets the zoom.
  Render-only: harnesses can't see it, Marco to check on screen. REWORK item removed.
- **THROTTLE WHILE PLANNING (2026-11-25, Marco's request):** W/S move the throttle while time is stopped, in REAL
  time: Game::Run calls physicsSystem.adjustThrottleWhilePlanning(universe, player_input.ship, frame_seconds) before
  updateTrajectories; Starship::adjustThrottle(direction, seconds) is shared with applyThrust (one tick of game time,
  bit-identical to before). The throttle counts as an EDIT: universeChangedSinceTraceStart compares it (lines wipe,
  press I again) and the placement preview restarts when it changes. Harness case 6: a planned 50% burn re-traced
  matches the real game for 5401 ship points, 0 differ; the ship doesn't move while it's set. Possible follow-up
  (Marco's call): re-trace automatically after a throttle change instead of waiting for I.
- **THROTTLE (2026-11-25, Marco's design):** Starship::m_throttle 0..1, W raises / S lowers it at
  THROTTLE_CHANGE_PER_SECOND (1.0 = off to full in 1 s of GAME time, per physics tick, so it freezes while planning),
  it HOLDS when released, no reverse thrust (S only reaches 0). Thrust = nose x THRUST_ACCELERATION x throttle.
  Predictions KEEP the current throttle (cutThrust removed: the prediction never runs applyThrust, so
  m_thrust_acceleration stays at the current burn) -> harness 'burn exactness': 5401 ship points vs the real game, 0
  differ. HUD: 'THROTTLE: 45%' under FPS. universe_refactor_baseline.txt REGENERATED for this intended change.
- **FUTURE FEATURES requested by Marco 2026-11-24 (after the planning-mode PR):**
  - **CAMERA (Marco: 'really need to fix the game camera'):** C = recenter the camera on the player's ship, and toggle
    between FREE camera (today's behaviour: middle-mouse pan, wheel zoom) and SOFT-CENTERED camera (follows the ship
    smoothly, e.g. ease toward it each frame instead of snapping). Today the camera is CameraState (offset, zoom,
    render_alpha) moved by InputSystem::updateCamera (zoomAroundCursor, panCamera, clampCameraToStarField). Things to
    decide: does panning in soft mode drop back to free? follow the ship's INTERPOLATED (drawn) position, or jitter;
    how the star-field clamp interacts with following; zoom around the ship vs the cursor while following.
  - **EDIT BODIES IN PLANNING MODE (Marco 2026-11-25):** drag planets to new positions while time is stopped, and
    their VELOCITY ARROWS should stay visible in planning mode and be draggable too (edit a body's velocity after it
    spawned). Fits the existing design: any edit already wipes the traced lines (universeChangedSinceTraceStart) and
    the placement preview machinery could show the edited body's new path live. Needs: picking a body under the
    cursor (world-space hit test), a drag state in PlayerInput (like SpawnProperties), arrows drawn per body (the arrow size
    bug in REWORK should be fixed first), and the edit applied by physics (UniverseState stays physics-owned).
  - **SPAWN FLAG UI (Marco 2026-11-25):** an expandable dropdown for the body property flags instead of independent
    checkboxes, so only meaningful combinations can be chosen. Goes with the REWORK 'Body property flags' item
    (decide the precedence rules first).
  - **SPAWN PRESETS:** two quick-swap slots for the instantiation settings. Press-and-HOLD 1 = save the current settings
    into slot 1, hold 2 = save into slot 2; a short TAP of 1/2 = load that slot, so you flip between two setups
    instantly. What a preset holds: SpawnSettings (Static, Shatterable, Accretable, Bounce, Ghost, Planet Ghost) and
    probably the radius + mass sliders (SpawnSettings::radius / mass since 2026-11-25, so ONE struct) -- Marco to confirm.
    Needs: hold-vs-tap detection (InputState has edges + levels; time a key held, e.g. > 0.5 s = save), sliders must
    SHOW loaded values (UISlider has setValue but no value source; the checkboxes already pull their SpawnSettings
    field every frame), some feedback that a save happened (sound / readout). Number keys must not clash with E
    (rebinding) later.
- Build/scripts follow-ups (the rework itself is DONE, see the done list): the [BUG]/[REFACTOR] items in REWORK's BUILD / TOOLING
  section (Windows setup: silent rmdir failure, unchecked submodules, implicit job count; Mac arch / minimum macOS version;
  Windows VC++ runtime).
  Git gotcha: restaging on Windows drops the Mac scripts' executable bit (100755 -> 100644); fix with
  `git update-index --chmod=+x Scripts/Apple/*.sh` before committing.
- Visor styling: row gaps, a selected-row look, a real "Realistic" coloring, then the actual Charge/Temperature views.
- E: action map, i.e. rebindable keys that replace hard-coded scancodes. DEPRECATED_InputState is already gone (2026-11-25):
  E maps keys onto PlayerCommands / ShipControls. All key bindings live in `translateGameInputs` on purpose, to make E easy.
- G: keyboard focus (Tab will then conflict in menus; the visor hotkey only runs in simulation scenes, on purpose).
- H: world-anchored UI.
- R1–R4: rendering into the engine.
  - R1 GPU upload helper (fixes cycle=false).
  - R2 pipeline helper (about -300 lines).
  - R3 engine UIRenderer.
  - R4 RenderContext.
  - SDLContext lifetime fix (SDL_Quit runs before members die; see REWORK).
- REWORK.txt: ~33 items (physics bugs, render risks, audio, build/tooling).
- From the 2026-10-09/10 session (Marco's open threads, roughly in the order they came up):
  - Right-click passes through the UI (UIRoot handles LEFT only): right-clicking a checkbox spawns a cluster behind it and now
    also shows the cluster preview over the UI. Logged in REWORK; the natural next fix.
  - Damage from Starship::lastImpactEnergy() (needs a threshold: resting on a planet reports ~7e6 per tick).
  - Fuel + Tsiolkovsky (thrust force / current mass, fuel burn per tick; seam = Starship::applyThrust); RCS translation fits,
    RCS rotation would mean real angular physics (reverses 'player-only rotation').
  - Shatter threshold relative to mutual escape speed (heavy planets that touch shatter from gravity alone; Marco accepted
    that for now), resting contact that never settles (ELASTIC_LOSS_FACTOR 0.994).
  - More spawn-panel checkboxes: Ghost + Planet Ghost DONE 2026-11-24 ('Global' is now 2 rows). Next one: one SpawnSettings
    field + one addSpawnCheckbox line; the Macro section still has room for a second row.
  - Despawn far-away debris to return particle budget (idea, Marco's call).
  - Small tidy-ups: UIColumn into its own file; StarshipVertex's unused square_size/color.
- **GPU PARTICLES (goal stated by Marco 2026-10-12, not started):** he wants HUNDREDS OF THOUSANDS of particles AND
  particle-particle gravity. He knows collisions are today's bottleneck. Claude's advice: GPU compute (SDL_GPU compute pipelines +
  shadercross), staged: (1) compute pass + render particles straight from a storage buffer, (2) brute-force N^2 self-gravity
  (tiled, ~20-30k max), (3) collisions = hash grid + GPU sort + GATHER contacts (each particle writes only itself, from a snapshot;
  replaces the serial impulse + break), (4) Barnes-Hut / particle-mesh gravity reusing the sort (100k+), (5) events back to the
  CPU (absorption, ship hits, shatter spawns) via small async readbacks. Open issues: Metal has NO double; IntegrationLines is
  exact only if the tracer runs the same code (later: same shaders on a buffer copy). Keep the CPU physics as the reference for
  harness comparisons.
  2026-10-19: target is now ~500k+. Marco's reference: tsp.grantkot.com ('Tutti Space Program', Grant Kot). Claude inspected it:
  it is NOT WebGPU. It's C++ MPM (material point method: particles interact through a sparse tiled GRID, no pair collisions)
  compiled to WASM, multithreaded on all CPU cores (Taskflow, float2), rendered with WebGL2. ~199k particles, self-gravity
  ('N-body', softening setting, method not visible in the WASM) + radial gravity, occupancy merge/split (adaptive particle
  count), temperature + phase changes. Claude's hidden browser pane showed 2-3 FPS (throttled, ignore). MARCO'S MACHINE
  (2026-10-21): 350k particles, 60 FPS, sim 20 ms. Read from its JS: the sim is DECOUPLED from rendering (_mpm_update_mt
  dispatches 1 step of 1/60 s to worker threads and returns; a callback marks it done; FPS = render rate from the last
  snapshot). 20 ms per 16.7 ms step = the sim runs at ~83% real time without stutter (slows down instead of dropping
  frames). Claude had guessed 'CPU tops out ~200k': WRONG, conceded. New advice: CPU threads first (parallel-for, SoA floats,
  gather collisions, grid gravity, sim thread + snapshot), GPU only if aiming for millions.
  2026-10-21 later: Marco DOES want GPU compute long term (>500k). Agreed shape: write the CPU stage AS KERNELS so the port is
  mechanical: flat float arrays (SoA, fixed capacity = future StructuredBuffers), one function per pass called as
  kernel(i) by parallelFor and writing only index i, no push_back/erase inside passes (alive flags + spawn slots), hash-based
  RNG (id + tick) instead of the shared randomDouble generator, counting-sort grid (count, prefix sum, scatter), no float
  atomics, events appended to an output list. Then port kernel by kernel, checking each against the CPU version. Get a GPU
  debugger first (Xcode Metal capture on the Mac or RenderDoc on Windows). OPEN DESIGN DECISION (Marco's): keep particles as hard discs (current rules, GPU hash
  grid + gather) or switch debris to MPM-style material (fluid/rubble/lava, fracture and temperature for free, bigger rewrite).

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

