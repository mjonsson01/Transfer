#!/usr/bin/env bash

set -euo pipefail

# Directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

LOCAL_DIR="${SCRIPT_DIR}/LocalShaderCross"
BUILD_ROOT="${SCRIPT_DIR}/../ShaderCross"
REPO_DIR="${BUILD_ROOT}/SDL_shadercross"

echo "Creating LocalShaderCross..."
mkdir -p "${LOCAL_DIR}"

echo "Creating build directory..."
mkdir -p "${BUILD_ROOT}"

if [ ! -d "${REPO_DIR}" ]; then
    echo "Cloning SDL_shadercross..."
    git clone https://github.com/libsdl-org/SDL_shadercross.git "${REPO_DIR}"
fi

cd "${REPO_DIR}"

echo "Initializing submodules..."
git submodule update --init --recursive

echo "Configuring CMake..."
cmake -S . -B build -DSDLSHADERCROSS_VENDORED=ON

echo "Building..."
cmake --build build --config Release

# Locate the executable
if [ -f build/shadercross ]; then
    EXECUTABLE="build/shadercross"
elif [ -f build/Release/shadercross ]; then
    EXECUTABLE="build/Release/shadercross"
else
    echo "Error: shadercross executable not found."
    exit 1
fi

echo "Copying executable..."
cp "${EXECUTABLE}" "${LOCAL_DIR}/"

echo "Done!"
echo "Executable copied to:"
echo "  ${LOCAL_DIR}/shadercross"j