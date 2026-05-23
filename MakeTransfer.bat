@REM @echo off
@REM setlocal enabledelayedexpansion

@REM set BUILD_DIR=build
@REM set CONFIG=Release

@REM REM --- Normalize argument ---
@REM set ARG=%1
@REM if "%ARG%"=="" set ARG=release

@REM REM --- Handle clean ---
@REM if /i "%ARG%"=="clean" (
@REM     echo Cleaning build directory...
@REM     rmdir /s /q "%BUILD_DIR%"
@REM     if %ERRORLEVEL% NEQ 0 (
@REM         echo Failed to clean build directory.
@REM         exit /b 1
@REM     )
@REM     echo Clean complete.
@REM     exit /b 0
@REM )

@REM REM --- Determine build type ---
@REM if /i "%ARG%"=="debug" set CONFIG=Debug
@REM if /i "%ARG%"=="release" set CONFIG=Release

@REM REM --- Auto-clean previous build ---
@REM if exist "%BUILD_DIR%" (
@REM     echo Auto-cleaning previous build...
@REM     rmdir /s /q "%BUILD_DIR%"
@REM     if %ERRORLEVEL% NEQ 0 (
@REM         echo Failed to clean build directory.
@REM         exit /b 1
@REM     )
@REM )

@REM REM --- Create build dir ---
@REM mkdir "%BUILD_DIR%"
@REM cd "%BUILD_DIR%"

@REM echo Configuring project (%CONFIG%)...
@REM cmake .. -DCMAKE_BUILD_TYPE=%CONFIG%
@REM if %ERRORLEVEL% NEQ 0 (
@REM     echo Configuration failed!
@REM     exit /b 1
@REM )

@REM echo Building TransferGame...
@REM cmake --build . --config %CONFIG%
@REM if %ERRORLEVEL% NEQ 0 (
@REM     echo Build failed!
@REM     exit /b 1
@REM )

@REM cd ..

@REM REM --- Determine executable path ---
@REM set EXE_PATH=%BUILD_DIR%\%CONFIG%\TransferGame.exe
@REM if not exist "%EXE_PATH%" set EXE_PATH=%BUILD_DIR%\TransferGame.exe
@REM if not exist "%EXE_PATH%" set EXE_PATH=%BUILD_DIR%\TransferGame\TransferGame.exe

@REM echo Build complete! Executable is at %EXE_PATH%

@REM REM --- Optionally launch ---
@REM set /p RUN="Press Enter to run the build, or type N + Enter to skip: "
@REM if "!RUN!"=="" (
@REM     echo Launching executable...
@REM     "%EXE_PATH%"
@REM ) else (
@REM     echo Skipping launch.
@REM )

@REM exit /b 0

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
    .\Transfer\src\HLSL\Body.vert.hlsl ^
    -o .\Transfer\Assets\Shaders\Body.vert.spv

if %ERRORLEVEL% NEQ 0 (
    echo Vertex shader compilation failed!
    exit /b 1
)

.\LocalShaderCross\shadercross.exe ^
    .\Transfer\src\HLSL\Body.frag.hlsl ^
    -o .\Transfer\Assets\Shaders\Body.frag.spv

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