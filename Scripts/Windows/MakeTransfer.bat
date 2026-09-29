@echo off
REM Builds the game into build\: lints the engine, compiles the HLSL shaders to SPIR-V, then configures and
REM builds with CMake. Needs Scripts\Windows\SetupDependencies.bat to have been run once.
REM   Scripts\Windows\MakeTransfer.bat           Release, incremental (fast iteration)
REM   Scripts\Windows\MakeTransfer.bat release   Release, full clean rebuild
REM   Scripts\Windows\MakeTransfer.bat debug     Debug, incremental
REM   Scripts\Windows\MakeTransfer.bat clean     delete build\ and stop
REM   set SKIP_TIDY=1                            skip the clang-tidy pass
setlocal enabledelayedexpansion

REM Run from the repo root (two folders up from this script) no matter where the script was launched from
cd /d "%~dp0..\.."

set BUILD_DIR=build
set CONFIG=Release
set DO_CLEAN=0

REM --- Handle Arguments ---
if /i "%1"=="clean" (
    echo Cleaning build directory...
    if exist "%BUILD_DIR%" rmdir /s /q "%BUILD_DIR%"
    exit /b 0
)

REM Bare invocation: Release, incremental (fast iteration).
REM 'release' explicitly: Release, full clean rebuild.
REM 'debug': Debug, incremental.
if /i "%1"=="debug" (
    set CONFIG=Debug
    set DO_CLEAN=0
) else if /i "%1"=="release" (
    set CONFIG=Release
    set DO_CLEAN=1
) else (
    set CONFIG=Release
    set DO_CLEAN=0
)

if %DO_CLEAN% EQU 0 (
    echo Incremental build requested for %CONFIG%.
)

REM --- Conditional Clean ---
if %DO_CLEAN% EQU 1 (
    if exist "%BUILD_DIR%" (
        echo Performing fresh build for %CONFIG%...
        rmdir /s /q "%BUILD_DIR%"
    )
)

if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"

REM =====================================================
REM Lint engine code (non-fatal; set SKIP_TIDY=1 to skip)
REM =====================================================
if not "%SKIP_TIDY%"=="1" (
    echo Running clang-tidy on DynamoEngine...
    call "%~dp0TidyEngine.bat"
    if errorlevel 1 echo WARNING: clang-tidy reported findings ^(see above^). Continuing build.
)

REM =====================================================
REM Compile HLSL -> SPIR-V
REM =====================================================
set "SHADERCROSS=ThirdParty\ShaderCross\bin\shadercross.exe"
set "SHADER_SRC=Transfer\src\HLSL"
set "SHADER_OUT=Transfer\Assets\Shaders"

if not exist "%SHADERCROSS%" (
    echo shadercross not found in ThirdParty\ShaderCross\. Run Scripts\Windows\SetupDependencies.bat first.
    exit /b 1
)
if not exist "%SHADER_OUT%" mkdir "%SHADER_OUT%"

echo Compiling shaders...
REM Every shader has a .vert.hlsl and a .frag.hlsl. "if errorlevel 1" is checked when each line RUNS (unlike
REM %ERRORLEVEL%, which a () block expands once, up front), so the first failing shader stops the build.
REM The failure jumps OUT of the loop to :shaderFailed: an "exit /b 1" inside a for loop can lose its exit code.
for %%S in (UnifiedGravBody TwinklingStar UIElement VelocityVector Starship) do (
    for %%T in (vert frag) do (
        "%SHADERCROSS%" "%SHADER_SRC%\%%S.%%T.hlsl" -o "%SHADER_OUT%\%%S.%%T.spv"
        if errorlevel 1 (
            set "FAILED_SHADER=%%S.%%T.hlsl"
            goto :shaderFailed
        )
    )
)
echo Shaders compiled successfully.

REM --- Configure and Build ---
echo Building TransferGame (%CONFIG%)...
cmake -S . -B %BUILD_DIR% -DCMAKE_BUILD_TYPE=%CONFIG%
if %ERRORLEVEL% NEQ 0 exit /b 1

cmake --build %BUILD_DIR% --config %CONFIG%
if %ERRORLEVEL% NEQ 0 exit /b 1

REM --- Locate Executable ---
set EXE_PATH=%BUILD_DIR%\%CONFIG%\TransferGame.exe
if not exist "%EXE_PATH%" set EXE_PATH=%BUILD_DIR%\TransferGame.exe

echo Build complete!
if exist "%EXE_PATH%" (
    set /p RUN="Press Enter to run, 'd' + Enter to run directly (see printouts), or any other key + Enter to skip: "
    if /i "!RUN!"=="" (
        "%EXE_PATH%"
    ) else if /i "!RUN!"=="d" (
        echo Running directly: %EXE_PATH%
        "%EXE_PATH%"
    ) else (
        echo Skipping launch.
    )
)
exit /b 0

:shaderFailed
echo Shader compilation failed: !FAILED_SHADER!
exit /b 1
