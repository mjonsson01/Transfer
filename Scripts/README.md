# Scripts

Every script in here runs from any folder: it first moves to the repo root on its own.
Windows scripts are in `Windows/` and macOS scripts are in `Apple/`. Both sets have the same scripts and take the same arguments (CheckDeterminism's `update` is Mac-only).

| Windows                                 | macOS                                | What it does                                                                          |
|-----------------------------------------|--------------------------------------|----------------------------------------------------------------------------------------|
| `Scripts\Windows\SetupDependencies.bat` | `Scripts/Apple/SetupDependencies.sh` | **Run once after cloning.** Installs SDL3, SDL3_ttf and shadercross into `ThirdParty/` |
| `Scripts\Windows\MakeTransfer.bat`      | `Scripts/Apple/MakeTransfer.sh`      | Lints the engine, compiles the shaders, builds the game into `build/`                  |
| `Scripts\Windows\RunTests.bat`          | `Scripts/Apple/RunTests.sh`          | Builds and runs the unit tests in `build-tests/`                                       |
| `Scripts\Windows\TidyEngine.bat`        | `Scripts/Apple/TidyEngine.sh`        | Runs clang-tidy over the engine (MakeTransfer also runs it)                            |
| `Scripts\Windows\CheckDeterminism.bat`  | `Scripts/Apple/CheckDeterminism.sh`  | Builds the game, runs its physics determinism check, compares with the Mac reference   |

## First-time setup

1. Install the build tools:
   - **Windows:** Visual Studio 2022 or its Build Tools, with the "Desktop development with C++" workload. Also CMake, git and Python 3.
   - **macOS:** the Xcode Command Line Tools (`xcode-select --install`, which includes clang, git and python3), plus CMake (`brew install cmake`).
2. Run `SetupDependencies`. The first run builds shadercross from source, which includes DirectXShaderCompiler, so **expect ~5 minutes**. On macOS it also builds SDL from source.
3. Run `MakeTransfer`.

`SetupDependencies` is safe to run again at any time. Anything already installed at the pinned version is skipped.

## Dependencies (`ThirdParty/`)

`ThirdParty/` is gitignored. Everything in it is produced by `SetupDependencies`, and CMake, TidyEngine and MakeTransfer look for SDL and shadercross **only** there. If something in there looks broken, delete the folder and run `SetupDependencies` again.

| Folder             | Windows                                                | macOS                                    |
|--------------------|--------------------------------------------------------|------------------------------------------|
| `SDL3`, `SDL3_ttf` | SDL's prebuilt Visual C++ packages, checked by SHA-256 | built from source at a pinned git commit |
| `ShaderCross`      | built from source at a pinned git commit               | built from source at a pinned git commit |

Each folder has a `VERSION.txt` that records what's installed.

**Upgrading a dependency:** change its pinned version (and SHA-256 or commit) at the top of **both** `SetupDependencies` scripts, then run them. The comments there explain where to get the values. Use stable SDL releases only: in SDL3's numbering, an odd minor version (3.3.x, 3.5.x) is a prerelease.

## Useful options

- `MakeTransfer [debug|release|clean]`: no argument gives an incremental Release build, `release` gives a full clean rebuild, `debug` gives an incremental Debug build, and `clean` deletes `build/`.
- `SKIP_TIDY=1` (`set SKIP_TIDY=1` on Windows) makes MakeTransfer skip the clang-tidy pass.
- `RunTests -R <name>`: any extra arguments go to ctest (`-R` runs only the tests whose names match).
- `TidyEngine --fix`: any extra arguments go to clang-tidy.
- `CheckDeterminism [nobuild]`: `nobuild` skips the build and runs the existing one. On the Mac, `CheckDeterminism update` saves the output as the new reference (`Notes/determinism_check_mac.txt`).

## Determinism check (Mac vs Windows)

A puzzle must play out the same on every machine, so the physics must give bit-identical results on the Mac and on Windows (Notes/PHYSICS_REQUIREMENTS.md). `TransferGame --determinism-check` plays six fixed physics scenes without opening a window and writes a checksum of every body's exact bits to `determinism_check.txt`. CheckDeterminism builds the game, runs that, and compares it with the reference made on the Mac:

1. **Mac:** `Scripts/Apple/CheckDeterminism.sh` should say MATCH. After a change that alters the physics on purpose, run `Scripts/Apple/CheckDeterminism.sh update` and commit the new `Notes/determinism_check_mac.txt`.
2. **Windows:** pull, then `Scripts\Windows\CheckDeterminism.bat`. MATCH means Windows and the Mac agree bit for bit. If not, the first scene that differs points at the physics code to look at.

Shaders are compiled **only** by `MakeTransfer`. After you change an `.hlsl` file, RunTests or a plain CMake build will still use the old compiled shaders.
