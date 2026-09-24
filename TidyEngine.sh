#!/bin/bash
# Runs clang-tidy over every DynamoEngine file -- headers included, which VS Code's automatic
# analysis skips (it only analyzes .cpp files). Rules come from Transfer/src/DynamoEngine/.clang-tidy.
#   ./TidyEngine.sh                      lint every engine file
#   ./TidyEngine.sh --fix                extra args go straight to clang-tidy (--fix applies suggested
#                                        renames, but only inside engine files -- review the diff, and
#                                        prefer the editor's Rename Symbol for anything used by the game)
#   CLANG_TIDY=/path/to/clang-tidy ./TidyEngine.sh   use a specific clang-tidy
# Exit code: 0 = clean, 1 = findings (or clang-tidy couldn't run).
set -uo pipefail # no -e: keep linting the remaining files after one reports findings
cd "$(dirname "$0")"

ENGINE_DIR="Transfer/src/DynamoEngine"

# --- Locate clang-tidy: $CLANG_TIDY, then PATH, then Homebrew LLVM, then the copy bundled with VS Code's C/C++ extension ---
find_clang_tidy()
{
    if [[ -n "${CLANG_TIDY:-}" ]]; then
        echo "$CLANG_TIDY"
    elif command -v clang-tidy >/dev/null 2>&1; then
        command -v clang-tidy
    elif [[ -x /opt/homebrew/opt/llvm/bin/clang-tidy ]]; then
        echo /opt/homebrew/opt/llvm/bin/clang-tidy
    else
        # Newest installed version of the extension wins
        ls -d "$HOME"/.vscode/extensions/ms-vscode.cpptools-*/LLVM/bin/clang-tidy 2>/dev/null | sort -V | tail -1
    fi
}

CLANG_TIDY_BIN="$(find_clang_tidy)"
if [[ -z "$CLANG_TIDY_BIN" || ! -x "$CLANG_TIDY_BIN" ]]; then
    echo "clang-tidy not found. Install the VS Code C/C++ extension, run 'brew install llvm', or set CLANG_TIDY." >&2
    exit 1
fi
echo "Using $CLANG_TIDY_BIN"

# --- Compiler flags (mirror CMakeLists.txt): -xc++ so headers are parsed as C++, SDL as a system include ---
COMPILE_FLAGS=(-xc++ -std=c++20 "-I$PWD/Transfer/src")
# A standalone clang-tidy doesn't know where Apple keeps the C++ standard library headers (<cmath> etc.);
# point it at the active Xcode / Command Line Tools SDK, the same one the real compiler uses.
if command -v xcrun >/dev/null 2>&1; then
    COMPILE_FLAGS+=(-isysroot "$(xcrun --show-sdk-path)")
fi
# VS Code's bundled clang-tidy ships without clang's own builtin headers (<stdarg.h> etc.), which normally
# live in <llvm>/lib/clang. If they're missing, borrow the system compiler's copy.
CLANG_TIDY_ROOT="$(cd "$(dirname "$CLANG_TIDY_BIN")/.." && pwd)"
if [[ ! -d "$CLANG_TIDY_ROOT/lib/clang" ]] && command -v clang >/dev/null 2>&1; then
    COMPILE_FLAGS+=(-resource-dir "$(clang -print-resource-dir)")
fi
if [[ -d /opt/homebrew/include ]]; then
    COMPILE_FLAGS+=(-isystem /opt/homebrew/include)
fi

# --- Lint each file on its own ---
# --header-filter='^$' reports only the file being linted: every header gets its own pass, so this keeps
# e.g. a Vector2 warning from being repeated once per file that includes Vector2.hpp.
files_checked=0
files_with_findings=0
while IFS= read -r -d '' file; do
    files_checked=$((files_checked + 1))
    output="$("$CLANG_TIDY_BIN" --quiet --header-filter='^$' "$@" "$file" -- "${COMPILE_FLAGS[@]}" 2>&1)"
    if printf '%s\n' "$output" | grep -qE "(warning|error):"; then
        files_with_findings=$((files_with_findings + 1))
        printf '%s\n\n' "$output"
    fi
done < <(find "$ENGINE_DIR" -type f \( -name '*.cpp' -o -name '*.hpp' -o -name '*.hh' -o -name '*.h' \) -print0 | sort -z)

# --- Summary ---
if [[ $files_with_findings -eq 0 ]]; then
    echo "TidyEngine: $files_checked files checked, all clean."
    exit 0
fi
echo "TidyEngine: $files_checked files checked, $files_with_findings with findings."
exit 1
