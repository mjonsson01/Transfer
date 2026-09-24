#!/bin/bash
# Builds and runs the unit tests in their own build folder (build-tests/), so the game's
# build/ folder and MakeTransfer.sh are never touched.
#   ./RunTests.sh                    run every test
#   ./RunTests.sh -R InputDevices    extra args go straight to ctest (-R filters by test name)
set -eEuo pipefail
trap 'echo "RunTests failed at line $LINENO: $BASH_COMMAND" >&2' ERR
cd "$(dirname "$0")"

BUILD_DIR="build-tests"

# Configure (first run downloads googletest into build-tests/_deps)
cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Debug -DTRANSFER_BUILD_TESTS=ON

if [[ -z "$(find Tests -name '*.cpp' 2>/dev/null)" ]]; then
    echo "googletest is set up, but there are no test sources in Tests/ yet."
    exit 0
fi

cmake --build "$BUILD_DIR" --target TransferTests
ctest --test-dir "$BUILD_DIR" --output-on-failure "$@"