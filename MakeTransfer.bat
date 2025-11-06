@echo off
setlocal
set BUILD_DIR=build

if "%1"=="clean" (
    echo Cleaning build directory...
    rmdir /s /q "%BUILD_DIR%"
    exit /b 0
)

if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"
cd "%BUILD_DIR%"

echo Configuring project...
cmake .. -DCMAKE_BUILD_TYPE=Release

echo Building TransferGame...
cmake --build . --config Release

cd ..
echo Build complete! Executable is in build\Release\
