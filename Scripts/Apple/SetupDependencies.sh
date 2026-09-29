#!/bin/bash
# Installs everything the macOS build needs into ThirdParty/ (gitignored). Run it once after cloning,
# and again whenever a pinned version below changes. It's safe to re-run: installed versions are skipped.
#   Scripts/Apple/SetupDependencies.sh
# What it installs (all built from source, so the FIRST run takes a while):
#   ThirdParty/SDL3, ThirdParty/SDL3_ttf  SDL itself (CMakeLists.txt and TidyEngine.sh look for SDL here and
#                                         nowhere else, the same as on Windows)
#   ThirdParty/ShaderCross                the HLSL shader compiler MakeTransfer.sh uses
# Needs: Xcode Command Line Tools (clang, git, python3: `xcode-select --install`) and CMake (`brew install cmake`).
# (Windows equivalent: Scripts\Windows\SetupDependencies.bat)
set -eEuo pipefail
trap 'echo "SetupDependencies failed at line $LINENO: $BASH_COMMAND" >&2' ERR

# Run from the repo root (two folders up from this script) no matter where the script was launched from
cd "$(dirname "$0")/../.."
THIRD_PARTY_DIR="$PWD/ThirdParty"

# --- Pinned versions: the ONLY place they're written down ---
# Everything is pinned to a git commit. A commit hash names its exact contents, including (through the
# submodule pointers stored in it) the exact FreeType/HarfBuzz/DirectXShaderCompiler sources, so it plays the
# role the SHA-256 plays on Windows. Git refuses to check out anything that doesn't match the hash.
# Stable SDL releases only: in SDL3's numbering an ODD minor version (3.3.x, 3.5.x) is a prerelease.
# To upgrade: look up the tag's commit with  git ls-remote https://github.com/libsdl-org/SDL.git refs/tags/release-X.Y.Z
SDL3_VERSION="3.4.16"
SDL3_COMMIT="fa2c02bb6e21974a89ea9824bc53c9932abe5f9c"         # tag release-3.4.16
SDL3_TTF_VERSION="3.2.2"
SDL3_TTF_COMMIT="a1ce3670aec736ecbf0936c43f2f0cc53aa61e5b"     # tag release-3.2.2
SHADERCROSS_COMMIT="1ff05bec573988a98ef9e0260b4da44f512b8367"  # no releases exist; main as of 2026-09-03

# --- Tools ---
for tool in git cmake clang python3; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "$tool not found. See the 'Needs:' line at the top of this script." >&2
        exit 1
    fi
done

# How many compile jobs to run at once. A bare "--parallel" means UNLIMITED with the Makefiles generator (plain
# "make -j"): DirectXShaderCompiler then starts a clang for nearly every file at once and runs the Mac out of memory.
# One job per CPU core, but at most one per 2 GB of RAM (a single DXC/LLVM file can take over 1 GB to compile).
CPU_COUNT="$(sysctl -n hw.ncpu)"
MEMORY_GB=$(( $(sysctl -n hw.memsize) / 1073741824 ))
BUILD_JOBS=$(( MEMORY_GB / 2 ))
if (( BUILD_JOBS > CPU_COUNT )); then BUILD_JOBS=$CPU_COUNT; fi
if (( BUILD_JOBS < 1 )); then BUILD_JOBS=1; fi
echo "Building with $BUILD_JOBS parallel jobs ($CPU_COUNT CPU cores, $MEMORY_GB GB RAM)."

mkdir -p "$THIRD_PARTY_DIR"
SOURCE_ROOT="$THIRD_PARTY_DIR/_src" # temporary: deleted after each install

# fetch_at_commit <git url> <commit> <folder>
#   Downloads exactly one commit (plus its submodules) instead of the whole history: --depth 1 = only that
#   one snapshot.
fetch_at_commit()
{
    local url="$1" commit="$2" folder="$3"
    rm -rf "$folder"
    mkdir -p "$folder"
    git -C "$folder" init --quiet
    git -C "$folder" remote add origin "$url"
    git -C "$folder" fetch --quiet --depth 1 origin "$commit"
    git -C "$folder" checkout --quiet FETCH_HEAD
    git -C "$folder" submodule update --quiet --init --recursive --depth 1
}

# install_from_git <name> <git url> <commit> [extra CMake arguments...]
#   Builds <name> at <commit> in Release and installs it to ThirdParty/<name>/, then deletes the sources.
#   ThirdParty/<name>/VERSION.txt records the commit, so a re-run with the same pin skips all of this.
#   "shift 3" drops the first three arguments; "$@" is then whatever extra CMake arguments are left.
install_from_git()
{
    local name="$1" url="$2" commit="$3"
    shift 3
    local install_dir="$THIRD_PARTY_DIR/$name"
    local source_dir="$SOURCE_ROOT/$name"

    if [[ -f "$install_dir/VERSION.txt" && "$(cat "$install_dir/VERSION.txt")" == "$commit" ]]; then
        echo "$name already installed."
        return 0
    fi

    echo "Downloading $name ($commit)..."
    fetch_at_commit "$url" "$commit" "$source_dir"

    echo "Building $name..."
    cmake -S "$source_dir" -B "$source_dir/build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$install_dir" "$@"
    cmake --build "$source_dir/build" --parallel "$BUILD_JOBS"
    rm -rf "$install_dir"
    cmake --install "$source_dir/build"

    echo "$commit" > "$install_dir/VERSION.txt"
    rm -rf "$source_dir"
    echo "$name installed."
}

# --- SDL3 ---
install_from_git SDL3 https://github.com/libsdl-org/SDL.git "$SDL3_COMMIT" \
    -DSDL_SHARED=ON -DSDL_STATIC=OFF -DSDL_TEST_LIBRARY=OFF -DSDL_TESTS=OFF -DSDL_EXAMPLES=OFF

# --- SDL3_ttf: VENDORED builds its own FreeType/HarfBuzz from the submodules, so nothing comes from Homebrew ---
install_from_git SDL3_ttf https://github.com/libsdl-org/SDL_ttf.git "$SDL3_TTF_COMMIT" \
    -DSDL3_DIR="$THIRD_PARTY_DIR/SDL3/lib/cmake/SDL3" -DSDLTTF_VENDORED=ON -DSDLTTF_SAMPLES=OFF -DBUILD_SHARED_LIBS=ON

# --- ShaderCross: VENDORED builds DirectXShaderCompiler + SPIRV-Cross from source (the slow part) ---
# Upstream doesn't give the installed shadercross a way to find its own libraries, so INSTALL_RPATH tells it to
# look in ../lib (next to its bin/ folder); libSDL3 is then copied in there too, like SDL3.dll on Windows.
install_from_git ShaderCross https://github.com/libsdl-org/SDL_shadercross.git "$SHADERCROSS_COMMIT" \
    -DSDL3_DIR="$THIRD_PARTY_DIR/SDL3/lib/cmake/SDL3" -DSDLSHADERCROSS_VENDORED=ON -DSDLSHADERCROSS_INSTALL=ON \
    -DSDLSHADERCROSS_STATIC=OFF -DCMAKE_INSTALL_RPATH="@executable_path/../lib"
cp -P "$THIRD_PARTY_DIR"/SDL3/lib/libSDL3*.dylib "$THIRD_PARTY_DIR/ShaderCross/lib/"

rm -rf "$SOURCE_ROOT"
echo
echo "Dependencies ready in ThirdParty/. Next: Scripts/Apple/MakeTransfer.sh or Scripts/Apple/RunTests.sh."
