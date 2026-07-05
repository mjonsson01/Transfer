@echo off
setlocal enabledelayedexpansion

set BUILD_DIR=build
set CONFIG=Release
set DO_CLEAN=0

REM --- Handle Arguments ---
if /i "%1"=="clean" (
    echo Cleaning build directory...
    if exist "%BUILD_DIR%" rmdir /s /q "%BUILD_DIR%"
    exit /b 0
)

if /i "%1"=="debug" (
    set CONFIG=Debug
    set DO_CLEAN=1
) else if /i "%1"=="release" (
    set CONFIG=Release
    set DO_CLEAN=1
) else (
    echo No configuration specified. Attempting incremental build...
    set DO_CLEAN=0
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
REM Compile HLSL -> SPIR-V
REM =====================================================

echo Compiling shaders...

.\LocalShaderCross\shadercross.exe ^
    .\Transfer\src\HLSL\UnifiedGravBody.vert.hlsl ^
    -o .\Transfer\Assets\Shaders\UnifiedGravBody.vert.spv

.\LocalShaderCross\shadercross.exe ^
    .\Transfer\src\HLSL\TwinklingStar.vert.hlsl ^
    -o .\Transfer\Assets\Shaders\TwinklingStar.vert.spv

.\LocalShaderCross\shadercross.exe ^
    .\Transfer\src\HLSL\UIElement.vert.hlsl ^
    -o .\Transfer\Assets\Shaders\UIElement.vert.spv

.\LocalShaderCross\shadercross.exe ^
.\Transfer\src\HLSL\UIElement.frag.hlsl ^
 -o .\Transfer\Assets\Shaders\UIElement.frag.spv

if %ERRORLEVEL% NEQ 0 (
    echo Vertex shader compilation failed!
    exit /b 1
)

.\LocalShaderCross\shadercross.exe ^
    .\Transfer\src\HLSL\UnifiedGravBody.frag.hlsl ^
    -o .\Transfer\Assets\Shaders\UnifiedGravBody.frag.spv

.\LocalShaderCross\shadercross.exe ^
    .\Transfer\src\HLSL\TwinklingStar.frag.hlsl ^
    -o .\Transfer\Assets\Shaders\TwinklingStar.frag.spv
if %ERRORLEVEL% NEQ 0 (
    echo Fragment shader compilation failed!
    exit /b 1
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
    set /p RUN="Press Enter to run, or N to skip: "
    if "!RUN!"=="" "%EXE_PATH%"
)