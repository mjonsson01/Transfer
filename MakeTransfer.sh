# #!/bin/bash
# set -eEuo pipefail

# BUILD_DIR="build"
# trap 'echo "Build failed at line $LINENO: $BASH_COMMAND" >&2' ERR

# # --- Normalize input to lowercase ---
# ARG="${1:-}"
# ARG_LOWER=$(echo "$ARG" | tr '[:upper:]' '[:lower:]')

# # --- Determine build type ---
# BUILD_TYPE="Debug"  # default
# case "$ARG_LOWER" in
#     debug) BUILD_TYPE="Debug" ;;
#     release) BUILD_TYPE="Release" ;;
#     clean)
#         rm -rf "$BUILD_DIR"
#         echo "Clean complete."
#         exit 0
#         ;;
# esac


# # =====================================================
# # Compile HLSL -> MSL
# # =====================================================

# echo "Compiling shaders..."

# SHADER_SRC="Transfer/src/HLSL"
# SHADER_OUT="Transfer/Assets/Shaders"
# mkdir -p "$SHADER_OUT"

# SHADERS=(UnifiedGravBody TwinklingStar UIElement VelocityVector)

# for name in "${SHADERS[@]}"; do
#     ./LocalShaderCross/shadercross "$SHADER_SRC/$name.vert.hlsl" -o "$SHADER_OUT/$name.vert.msl"
#     ./LocalShaderCross/shadercross "$SHADER_SRC/$name.frag.hlsl" -o "$SHADER_OUT/$name.frag.msl"
# done

# echo "Shaders compiled successfully."


# echo "Cleaning and configuring $BUILD_TYPE build..."

# # --- Always clean first ---
# rm -rf "$BUILD_DIR"
# mkdir -p "$BUILD_DIR"
# cd "$BUILD_DIR"


# # --- Run CMake and build ---
# cmake .. -DCMAKE_BUILD_TYPE=$BUILD_TYPE 
# # cmake --build . --config $BUILD_TYPE
# cmake --build .

# cd ..

# # --- Finish message ---
# if [[ "$(uname)" == "Darwin" && "$BUILD_TYPE" == "Release" ]]; then
#     BUILD_PATH="$BUILD_DIR/TransferGame.app"
#     echo "Build complete. App bundle is in $BUILD_PATH"
# elif [[ "$BUILD_TYPE" == "Debug" ]]; then
#     if [[ "$(uname)" == "Darwin" ]]; then
#         BUILD_PATH="$BUILD_DIR/TransferGame"
#     else
#         BUILD_PATH="$BUILD_DIR/TransferGame.exe"
#     fi
#     echo "Build complete. Executable is in $BUILD_PATH ($BUILD_TYPE)"
# else
#     # Windows/Linux Release
#     BUILD_PATH="$BUILD_DIR/TransferGame.exe"
#     echo "Build complete. Executable is in $BUILD_PATH ($BUILD_TYPE)"
# fi

# # --- Prompt to run build ---
# echo
# read -p "Press Enter to run the build, or any other key + Enter to skip: " input
# if [[ -z "$input" ]]; then
#     echo "Launching..."
#     if [[ "$(uname)" == "Darwin" && "$BUILD_TYPE" == "Release" ]]; then
#         open "$BUILD_PATH"
#     else
#         "$BUILD_PATH"
#     fi
# else
#     echo "Skipping launch."
# fi
#!/bin/bash
set -eEuo pipefail

BUILD_DIR="build"
trap 'echo "Build failed at line $LINENO: $BASH_COMMAND" >&2' ERR

# --- Normalize input to lowercase ---
ARG="${1:-}"
ARG_LOWER=$(echo "$ARG" | tr '[:upper:]' '[:lower:]')

# --- Determine build type and whether to force a scratch rebuild ---
# Bare invocation: Release, incremental (fast iteration).
# 'release' explicitly: Release, full clean rebuild.
# 'debug': Debug, incremental.
# 'clean': wipe build dir and exit.
BUILD_TYPE="Release"
DO_CLEAN=false

case "$ARG_LOWER" in
    "")
        BUILD_TYPE="Release"
        DO_CLEAN=false
        ;;
    debug)
        BUILD_TYPE="Debug"
        DO_CLEAN=false
        ;;
    release)
        BUILD_TYPE="Release"
        DO_CLEAN=true
        ;;
    clean)
        rm -rf "$BUILD_DIR"
        echo "Clean complete."
        exit 0
        ;;
    *)
        echo "Unknown argument: $ARG" >&2
        echo "Usage: $0 [debug|release|clean]" >&2
        exit 1
        ;;
esac


# =====================================================
# Compile HLSL -> MSL
# =====================================================

echo "Compiling shaders..."

SHADER_SRC="Transfer/src/HLSL"
SHADER_OUT="Transfer/Assets/Shaders"
mkdir -p "$SHADER_OUT"

SHADERS=(UnifiedGravBody TwinklingStar UIElement VelocityVector)

for name in "${SHADERS[@]}"; do
    ./LocalShaderCross/shadercross "$SHADER_SRC/$name.vert.hlsl" -o "$SHADER_OUT/$name.vert.msl"
    ./LocalShaderCross/shadercross "$SHADER_SRC/$name.frag.hlsl" -o "$SHADER_OUT/$name.frag.msl"
done

echo "Shaders compiled successfully."


if [[ "$DO_CLEAN" == true ]]; then
    echo "Full scratch rebuild: configuring $BUILD_TYPE build..."
    rm -rf "$BUILD_DIR"
else
    echo "Incremental build: configuring $BUILD_TYPE build..."
fi

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# --- Run CMake and build ---
cmake .. -DCMAKE_BUILD_TYPE=$BUILD_TYPE
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

# --- Resolve the actual runnable binary (differs from BUILD_PATH when it's a macOS .app bundle) ---
if [[ "$(uname)" == "Darwin" && "$BUILD_TYPE" == "Release" ]]; then
    EXEC_PATH="$BUILD_PATH/Contents/MacOS/TransferGame"
else
    EXEC_PATH="$BUILD_PATH"
fi

# --- Prompt to run build ---
echo
read -p "Press Enter to run the build, 'd' + Enter to run it directly (see printouts), or any other key + Enter to skip: " input

case "$input" in
    "")
        echo "Launching..."
        if [[ "$(uname)" == "Darwin" && "$BUILD_TYPE" == "Release" ]]; then
            open "$BUILD_PATH"
        else
            "$EXEC_PATH"
        fi
        ;;
    d|D)
        echo "Running directly: $EXEC_PATH"
        "$EXEC_PATH"
        ;;
    *)
        echo "Skipping launch."
        ;;
esac