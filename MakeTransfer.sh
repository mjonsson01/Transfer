#!/bin/bash
set -eEuo pipefail

BUILD_DIR="build"
trap 'echo "Build failed at line $LINENO: $BASH_COMMAND" >&2' ERR

# --- Normalize input to lowercase ---
ARG="${1:-}"
ARG_LOWER=$(echo "$ARG" | tr '[:upper:]' '[:lower:]')

# --- Determine build type ---
BUILD_TYPE="Debug"  # default
case "$ARG_LOWER" in
    debug) BUILD_TYPE="Debug" ;;
    release) BUILD_TYPE="Release" ;;
    clean)
        rm -rf "$BUILD_DIR"
        echo "Clean complete."
        exit 0
        ;;
esac

echo "Cleaning and configuring $BUILD_TYPE build..."

# --- Always clean first ---
rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# --- Run CMake and build ---
cmake .. -DCMAKE_BUILD_TYPE=$BUILD_TYPE 
# cmake --build . --config $BUILD_TYPE
cmake --build .

cd ..

# --- Finish message ---
if [[ "$(uname)" == "Darwin" && "$BUILD_TYPE" == "Release" ]]; then
    BUILD_PATH="$BUILD_DIR/TransferGame.app"
    echo "Build complete. App bundle is in $BUILD_PATH"
elif [[ "$BUILD_TYPE" == "Debug" ]]; then
    if [[ "$(uname)" == "Darwin" ]]; then
        BUILD_PATH="$BUILD_DIR/TransferGame"
    else
        BUILD_PATH="$BUILD_DIR/TransferGame.exe"
    fi
    echo "Build complete. Executable is in $BUILD_PATH ($BUILD_TYPE)"
else
    # Windows/Linux Release
    BUILD_PATH="$BUILD_DIR/TransferGame.exe"
    echo "Build complete. Executable is in $BUILD_PATH ($BUILD_TYPE)"
fi

# --- Prompt to run build ---
echo
read -p "Press Enter to run the build, or any other key + Enter to skip: " input
if [[ -z "$input" ]]; then
    echo "Launching..."
    if [[ "$(uname)" == "Darwin" && "$BUILD_TYPE" == "Release" ]]; then
        open "$BUILD_PATH"
    else
        "$BUILD_PATH"
    fi
else
    echo "Skipping launch."
fi