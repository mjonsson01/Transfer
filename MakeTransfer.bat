@REM @echo off
@REM setlocal enabledelayedexpansion

@REM set BUILD_DIR=build
@REM set CONFIG=Release
@REM set DO_CLEAN=0

@REM REM --- Handle Arguments ---
@REM if /i "%1"=="clean" (
@REM     echo Cleaning build directory...
@REM     if exist "%BUILD_DIR%" rmdir /s /q "%BUILD_DIR%"
@REM     exit /b 0
@REM )

@REM if /i "%1"=="debug" (
@REM     set CONFIG=Debug
@REM     set DO_CLEAN=1
@REM ) else if /i "%1"=="release" (
@REM     set CONFIG=Release
@REM     set DO_CLEAN=1
@REM ) else (
@REM     echo No configuration specified. Attempting incremental build...
@REM     set DO_CLEAN=0
@REM )

@REM REM --- Conditional Clean ---
@REM if %DO_CLEAN% EQU 1 (
@REM     if exist "%BUILD_DIR%" (
@REM         echo Performing fresh build for %CONFIG%...
@REM         rmdir /s /q "%BUILD_DIR%"
@REM     )
@REM )

@REM if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"
@REM REM =====================================================
@REM REM Compile HLSL -> SPIR-V
@REM REM =====================================================

@REM echo Compiling shaders...

@REM .\LocalShaderCross\shadercross.exe ^
@REM     .\Transfer\src\HLSL\UnifiedGravBody.vert.hlsl ^
@REM     -o .\Transfer\Assets\Shaders\UnifiedGravBody.vert.spv

@REM .\LocalShaderCross\shadercross.exe ^
@REM     .\Transfer\src\HLSL\TwinklingStar.vert.hlsl ^
@REM     -o .\Transfer\Assets\Shaders\TwinklingStar.vert.spv

@REM .\LocalShaderCross\shadercross.exe ^
@REM     .\Transfer\src\HLSL\UIElement.vert.hlsl ^
@REM     -o .\Transfer\Assets\Shaders\UIElement.vert.spv

@REM .\LocalShaderCross\shadercross.exe ^
@REM     .\Transfer\src\HLSL\VelocityVector.vert.hlsl ^
@REM     -o .\Transfer\Assets\Shaders\VelocityVector.vert.spv

@REM if %ERRORLEVEL% NEQ 0 (
@REM     echo Vertex shader compilation failed!
@REM     exit /b 1
@REM )

@REM .\LocalShaderCross\shadercross.exe ^
@REM     .\Transfer\src\HLSL\UnifiedGravBody.frag.hlsl ^
@REM     -o .\Transfer\Assets\Shaders\UnifiedGravBody.frag.spv

@REM .\LocalShaderCross\shadercross.exe ^
@REM     .\Transfer\src\HLSL\TwinklingStar.frag.hlsl ^
@REM     -o .\Transfer\Assets\Shaders\TwinklingStar.frag.spv

@REM .\LocalShaderCross\shadercross.exe ^
@REM     .\Transfer\src\HLSL\UIElement.frag.hlsl ^
@REM      -o .\Transfer\Assets\Shaders\UIElement.frag.spv

@REM .\LocalShaderCross\shadercross.exe ^
@REM     .\Transfer\src\HLSL\VelocityVector.frag.hlsl ^
@REM     -o .\Transfer\Assets\Shaders\VelocityVector.frag.spv


@REM if %ERRORLEVEL% NEQ 0 (
@REM     echo Fragment shader compilation failed!
@REM     exit /b 1
@REM )

@REM echo Shaders compiled successfully.
@REM REM --- Configure and Build ---
@REM echo Building TransferGame (%CONFIG%)...
@REM cmake -S . -B %BUILD_DIR% -DCMAKE_BUILD_TYPE=%CONFIG%
@REM if %ERRORLEVEL% NEQ 0 exit /b 1

@REM cmake --build %BUILD_DIR% --config %CONFIG%
@REM if %ERRORLEVEL% NEQ 0 exit /b 1

@REM REM --- Locate Executable ---
@REM set EXE_PATH=%BUILD_DIR%\%CONFIG%\TransferGame.exe
@REM if not exist "%EXE_PATH%" set EXE_PATH=%BUILD_DIR%\TransferGame.exe

@REM echo Build complete!
@REM if exist "%EXE_PATH%" (
@REM     set /p RUN="Press Enter to run, or N to skip: "
@REM     if "!RUN!"=="" "%EXE_PATH%"
@REM )


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
    .\Transfer\src\HLSL\VelocityVector.vert.hlsl ^
    -o .\Transfer\Assets\Shaders\VelocityVector.vert.spv

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

.\LocalShaderCross\shadercross.exe ^
    .\Transfer\src\HLSL\UIElement.frag.hlsl ^
     -o .\Transfer\Assets\Shaders\UIElement.frag.spv

.\LocalShaderCross\shadercross.exe ^
    .\Transfer\src\HLSL\VelocityVector.frag.hlsl ^
    -o .\Transfer\Assets\Shaders\VelocityVector.frag.spv


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