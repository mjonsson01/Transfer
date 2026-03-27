@REM @REM @echo off
@REM @REM setlocal
@REM @REM set BUILD_DIR=build

@REM @REM if "%1"=="clean" (
@REM @REM     echo Cleaning build directory...
@REM @REM     rmdir /s /q "%BUILD_DIR%"
@REM @REM     if %ERRORLEVEL% NEQ 0 (
@REM @REM         echo Failed to clean build directory.
@REM @REM         exit /b 1
@REM @REM     )
@REM @REM     exit /b 0
@REM @REM )

@REM @REM if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"
@REM @REM cd "%BUILD_DIR%"

@REM @REM echo Configuring project...
@REM @REM cmake .. 
@REM @REM if %ERRORLEVEL% NEQ 0 (
@REM @REM     echo Configuration failed!
@REM @REM     exit /b 1
@REM @REM )

@REM @REM echo Building TransferGame...
@REM @REM cmake --build . --config Release
@REM @REM if %ERRORLEVEL% NEQ 0 (
@REM @REM     echo Build failed!
@REM @REM     exit /b 1
@REM @REM )

@REM @REM cd ..
@REM @REM echo Build complete! Executable is in %BUILD_DIR%\Release\
@REM @REM exit /b 0


@REM @echo off
@REM setlocal enabledelayedexpansion

@REM set BUILD_DIR=build
@REM set CONFIG=Release

@REM REM --- Parse command line argument ---
@REM if "%1"=="debug" set CONFIG=Debug
@REM if "%1"=="release" set CONFIG=Release
@REM if "%1"=="clean" (
@REM     echo Cleaning build directory...
@REM     rmdir /s /q "%BUILD_DIR%"
@REM     if %ERRORLEVEL% NEQ 0 (
@REM         echo Failed to clean build directory.
@REM         exit /b 1
@REM     )
@REM     exit /b 0
@REM )

@REM REM --- Create build dir ---
@REM if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"
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


@REM @echo off
@REM setlocal enabledelayedexpansion

@REM set BUILD_DIR=build
@REM set CONFIG=Release

@REM REM --- Parse command line argument ---
@REM if /i "%1"=="debug" set CONFIG=Debug
@REM if /i "%1"=="release" set CONFIG=Release
@REM if /i "%1"=="clean" (
@REM     echo Cleaning build directory...
@REM     rmdir /s /q "%BUILD_DIR%"
@REM     if %ERRORLEVEL% NEQ 0 (
@REM         echo Failed to clean build directory.
@REM         exit /b 1
@REM     )
@REM     echo Clean complete.
@REM     exit /b 0
@REM )

@REM REM --- Always auto-clean first ---
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

REM --- Normalize argument ---
set ARG=%1
if "%ARG%"=="" set ARG=release

REM --- Handle clean ---
if /i "%ARG%"=="clean" (
    echo Cleaning build directory...
    rmdir /s /q "%BUILD_DIR%"
    if %ERRORLEVEL% NEQ 0 (
        echo Failed to clean build directory.
        exit /b 1
    )
    echo Clean complete.
    exit /b 0
)

REM --- Determine build type ---
if /i "%ARG%"=="debug" set CONFIG=Debug
if /i "%ARG%"=="release" set CONFIG=Release

REM --- Auto-clean previous build ---
if exist "%BUILD_DIR%" (
    echo Auto-cleaning previous build...
    rmdir /s /q "%BUILD_DIR%"
    if %ERRORLEVEL% NEQ 0 (
        echo Failed to clean build directory.
        exit /b 1
    )
)

REM --- Create build dir ---
mkdir "%BUILD_DIR%"
cd "%BUILD_DIR%"

echo Configuring project (%CONFIG%)...
cmake .. -DCMAKE_BUILD_TYPE=%CONFIG%
if %ERRORLEVEL% NEQ 0 (
    echo Configuration failed!
    exit /b 1
)

echo Building TransferGame...
cmake --build . --config %CONFIG%
if %ERRORLEVEL% NEQ 0 (
    echo Build failed!
    exit /b 1
)

cd ..

REM --- Determine executable path ---
set EXE_PATH=%BUILD_DIR%\%CONFIG%\TransferGame.exe
if not exist "%EXE_PATH%" set EXE_PATH=%BUILD_DIR%\TransferGame.exe
if not exist "%EXE_PATH%" set EXE_PATH=%BUILD_DIR%\TransferGame\TransferGame.exe

echo Build complete! Executable is at %EXE_PATH%

REM --- Optionally launch ---
set /p RUN="Press Enter to run the build, or type N + Enter to skip: "
if "!RUN!"=="" (
    echo Launching executable...
    "%EXE_PATH%"
) else (
    echo Skipping launch.
)

exit /b 0