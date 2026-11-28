@echo off
REM Builds the game (incremental Release, into build\), runs its determinism check (TransferGame --determinism-check:
REM six fixed physics scenes, no window, ~1 s) and compares the output with the Mac's saved reference,
REM Notes\determinism_check_mac.txt. Identical = Windows plays the physics out bit for bit the same as the Mac.
REM   Scripts\Windows\CheckDeterminism.bat            build, run, compare with the reference
REM   Scripts\Windows\CheckDeterminism.bat nobuild    run the existing build, compare
REM (The reference is made on the Mac: Scripts/Apple/CheckDeterminism.sh update. Pull it before comparing.)
REM Needs MakeTransfer.bat to have run once (it compiles the shaders the build copies next to the game).
setlocal

REM Run from the repo root (two folders up from this script) no matter where the script was launched from
cd /d "%~dp0..\.."

set BUILD_DIR=build
set CONFIG=Release
set REFERENCE=Notes\determinism_check_mac.txt
REM Written by the game into the folder it runs in (gitignored)
set OUTPUT=determinism_check.txt

if not "%1"=="" if /i not "%1"=="nobuild" (
    echo Unknown argument: %1
    echo Usage: Scripts\Windows\CheckDeterminism.bat [nobuild]
    exit /b 1
)

REM --- Build (the same Release build MakeTransfer makes, without its shader and lint steps or its run prompt) ---
if /i not "%1"=="nobuild" (
    cmake -S . -B %BUILD_DIR% -DCMAKE_BUILD_TYPE=%CONFIG%
    if errorlevel 1 (
        echo CheckDeterminism failed: CMake configure step.
        exit /b 1
    )
    REM --config is required by multi-config generators (Visual Studio); single-config generators ignore it.
    cmake --build %BUILD_DIR% --config %CONFIG% --target TransferGame
    if errorlevel 1 (
        echo CheckDeterminism failed: build step.
        exit /b 1
    )
)

REM --- Locate the executable (Visual Studio puts it in build\Release\, single-config generators in build\) ---
set EXE_PATH=%BUILD_DIR%\%CONFIG%\TransferGame.exe
if not exist "%EXE_PATH%" set EXE_PATH=%BUILD_DIR%\TransferGame.exe
if not exist "%EXE_PATH%" (
    echo No Release build found in %BUILD_DIR%. Run Scripts\Windows\MakeTransfer.bat first.
    exit /b 1
)

REM --- Run the check ---
"%EXE_PATH%" --determinism-check
if errorlevel 1 (
    echo CheckDeterminism failed: the check itself did not finish.
    exit /b 1
)

REM --- Compare with the Mac's reference (fc /L: compare as text; it prints the differing lines) ---
echo.
fc /L "%REFERENCE%" "%OUTPUT%"
if errorlevel 1 (
    echo DIFFERENT from %REFERENCE%: Windows and the Mac don't produce the same physics.
    echo The first scene that differs points at the code to look at.
    exit /b 1
)
echo MATCH: identical to %REFERENCE%
exit /b 0
