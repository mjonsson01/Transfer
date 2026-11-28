#!/bin/bash
# Builds the game (incremental Release, into build/), runs its determinism check (TransferGame --determinism-check:
# six fixed physics scenes, no window, ~1 s) and compares the output with the saved reference,
# Notes/determinism_check_mac.txt. Identical = the physics plays out bit for bit the same as when the reference was
# made. On Windows, Scripts\Windows\CheckDeterminism.bat compares with the same file: that's the Mac vs Windows check.
#   Scripts/Apple/CheckDeterminism.sh            build, run, compare with the reference
#   Scripts/Apple/CheckDeterminism.sh nobuild    run the existing build, compare
#   Scripts/Apple/CheckDeterminism.sh update     build, run, and SAVE the output as the new reference
#                                                (after an intended physics change; commit the new reference)
# Needs MakeTransfer.sh to have run once (it compiles the shaders the app bundle copies in).
set -eEuo pipefail
trap 'echo "CheckDeterminism failed at line $LINENO: $BASH_COMMAND" >&2' ERR
cd "$(dirname "$0")/../.." # the repo root, two folders up from this script

BUILD_DIR="build"
REFERENCE="Notes/determinism_check_mac.txt"
OUTPUT="determinism_check.txt" # written by the game into the folder it runs in (gitignored)
EXECUTABLE="$BUILD_DIR/TransferGame.app/Contents/MacOS/TransferGame"

ARG_LOWER=$(echo "${1:-}" | tr '[:upper:]' '[:lower:]')
case "$ARG_LOWER" in
    "" | update | nobuild) ;;
    *)
        echo "Unknown argument: $1" >&2
        echo "Usage: $0 [nobuild|update]" >&2
        exit 1
        ;;
esac

# --- Build (the same Release build MakeTransfer makes, without its shader and lint steps or its run prompt) ---
if [[ "$ARG_LOWER" != "nobuild" ]]; then
    cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
    cmake --build "$BUILD_DIR" --target TransferGame
fi
if [[ ! -x "$EXECUTABLE" ]]; then
    echo "No Release build at $EXECUTABLE. Run Scripts/Apple/MakeTransfer.sh first." >&2
    exit 1
fi

# --- Run the check ---
"$EXECUTABLE" --determinism-check

# --- Save as the reference, or compare with it ---
if [[ "$ARG_LOWER" == "update" ]]; then
    cp "$OUTPUT" "$REFERENCE"
    echo
    echo "Saved as the new reference: $REFERENCE (commit it, then run CheckDeterminism on Windows)."
    exit 0
fi

echo
if diff "$REFERENCE" "$OUTPUT"; then
    echo "MATCH: identical to $REFERENCE"
else
    echo "DIFFERENT from $REFERENCE (lines above: < reference, > this run)." >&2
    echo "The first scene that differs points at the code to look at. If the physics changed on purpose, run" >&2
    echo "Scripts/Apple/CheckDeterminism.sh update and commit the new reference." >&2
    exit 1
fi
