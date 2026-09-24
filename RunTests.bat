@echo off
REM Builds and runs the unit tests in their own build folder (build-tests\), so the game's
REM build\ folder and MakeTransfer.bat are never touched.
REM   RunTests.bat                    run every test
REM   RunTests.bat -R InputDevices    extra args go straight to ctest (-R filters by test name)
setlocal

REM Run from the repo root no matter where the script was launched from
cd /d "%~dp0"

set BUILD_DIR=build-tests
set CONFIG=Debug

REM --- Configure (first run downloads googletest into build-tests\_deps) ---
cmake -S . -B %BUILD_DIR% -DCMAKE_BUILD_TYPE=%CONFIG% -DTRANSFER_BUILD_TESTS=ON
if %ERRORLEVEL% NEQ 0 (
    echo RunTests failed: CMake configure step.
    exit /b 1
)

REM --- Bail out early if there are no tests yet (dir returns an error when nothing matches) ---
dir /s /b "Tests\*.cpp" >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo googletest is set up, but there are no test sources in Tests\ yet.
    exit /b 0
)

REM --- Build the test executable ---
REM --config is required by multi-config generators (Visual Studio); single-config generators ignore it.
cmake --build %BUILD_DIR% --target TransferTests --config %CONFIG%
if %ERRORLEVEL% NEQ 0 (
    echo RunTest failed: build step.
    exit /b 1
)

REM --- Run the tests ---
REM -C picks the configuration for multi-config generators; %* forwards any extra arguments to ctest.
ctest --test-dir %BUILD_DIR% -C %CONFIG% --output-on-failure %*
if %ERRORLEVEL% NEQ 0 (
    echo RunTests failed: one or more tests failed.
    exit /b 1
)

echo All tests passed.
exit /b 0
