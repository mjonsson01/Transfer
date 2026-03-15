@echo off
setlocal
set BUILD_DIR=build

if "%1"=="clean" (
    echo Cleaning build directory...
    rmdir /s /q "%BUILD_DIR%"
    if %ERRORLEVEL% NEQ 0 (
        echo Failed to clean build directory.
        exit /b 1
    )
    exit /b 0
)

if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"
cd "%BUILD_DIR%"

echo Configuring project...
cmake .. 
if %ERRORLEVEL% NEQ 0 (
    echo Configuration failed!
    exit /b 1
)

echo Building TransferGame...
cmake --build . --config Release
if %ERRORLEVEL% NEQ 0 (
    echo Build failed!
    exit /b 1
)

cd ..
echo Build complete! Executable is in %BUILD_DIR%\Release\
exit /b 0