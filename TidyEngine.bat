@echo off
REM Runs clang-tidy over every DynamoEngine file -- headers included, which VS Code's automatic
REM analysis skips (it only analyzes .cpp files). Rules come from Transfer\src\DynamoEngine\.clang-tidy.
REM   TidyEngine.bat                       lint every engine file
REM   TidyEngine.bat --fix                 extra args go straight to clang-tidy (--fix applies suggested
REM                                        renames, but only inside engine files -- review the diff, and
REM                                        prefer the editor's Rename Symbol for anything used by the game)
REM   set CLANG_TIDY=C:\path\clang-tidy.exe   use a specific clang-tidy
REM Exit code: 0 = clean, 1 = findings (or clang-tidy couldn't run).
setlocal enabledelayedexpansion

REM Run from the repo root no matter where the script was launched from
cd /d "%~dp0"

set "ENGINE_DIR=Transfer\src\DynamoEngine"
set "OUTPUT_FILE=%TEMP%\TidyEngine_output.txt"

REM --- Locate clang-tidy: %CLANG_TIDY%, then PATH, then LLVM's installer, then VS Code's C/C++ extension ---
set "CLANG_TIDY_BIN=%CLANG_TIDY%"
if not defined CLANG_TIDY_BIN for %%P in (clang-tidy.exe) do set "CLANG_TIDY_BIN=%%~$PATH:P"
if not defined CLANG_TIDY_BIN if exist "C:\Program Files\LLVM\bin\clang-tidy.exe" set "CLANG_TIDY_BIN=C:\Program Files\LLVM\bin\clang-tidy.exe"
if not defined CLANG_TIDY_BIN (
    REM Folders are listed alphabetically, so the last match is (approximately) the newest extension version
    for /d %%D in ("%USERPROFILE%\.vscode\extensions\ms-vscode.cpptools-*") do (
        if exist "%%D\LLVM\bin\clang-tidy.exe" set "CLANG_TIDY_BIN=%%D\LLVM\bin\clang-tidy.exe"
    )
)
if not defined CLANG_TIDY_BIN (
    echo clang-tidy not found. Install the VS Code C/C++ extension or LLVM, or set CLANG_TIDY.
    exit /b 1
)
echo Using %CLANG_TIDY_BIN%

REM --- Compiler flags (mirror CMakeLists.txt): -xc++ so headers are parsed as C++, SDL as a system include ---
set "COMPILE_FLAGS=-xc++ -std=c++20 -I"%~dp0Transfer\src" -isystem "C:\SDL3\include" -isystem "%~dp0..\..\SDL3_TTF\include""

REM VS Code's bundled clang-tidy ships without clang's own builtin headers (<stdarg.h> etc.), which normally
REM live in <llvm>\lib\clang. If they're missing and a full clang is on PATH, borrow its copy.
for %%B in ("%CLANG_TIDY_BIN%") do set "CLANG_TIDY_ROOT=%%~dpB.."
if not exist "%CLANG_TIDY_ROOT%\lib\clang" (
    for /f "usebackq delims=" %%R in (`clang -print-resource-dir 2^>nul`) do set "COMPILE_FLAGS=!COMPILE_FLAGS! -resource-dir "%%R""
)

REM --- Lint each file on its own ---
REM --header-filter=^$ reports only the file being linted: every header gets its own pass, so this keeps
REM e.g. a Vector2 warning from being repeated once per file that includes Vector2.hpp.
set /a FILES_CHECKED=0
set /a FILES_WITH_FINDINGS=0
for /r "%ENGINE_DIR%" %%F in (*.cpp *.hpp *.hh *.h) do (
    set /a FILES_CHECKED+=1
    "%CLANG_TIDY_BIN%" --quiet "--header-filter=^$" %* "%%F" -- %COMPILE_FLAGS% > "%OUTPUT_FILE%" 2>&1
    findstr /c:"warning:" /c:"error:" "%OUTPUT_FILE%" >nul
    if !ERRORLEVEL! EQU 0 (
        set /a FILES_WITH_FINDINGS+=1
        type "%OUTPUT_FILE%"
        echo.
    )
)
if exist "%OUTPUT_FILE%" del "%OUTPUT_FILE%"

REM --- Summary ---
if %FILES_WITH_FINDINGS% EQU 0 (
    echo TidyEngine: %FILES_CHECKED% files checked, all clean.
    exit /b 0
)
echo TidyEngine: %FILES_CHECKED% files checked, %FILES_WITH_FINDINGS% with findings.
exit /b 1
