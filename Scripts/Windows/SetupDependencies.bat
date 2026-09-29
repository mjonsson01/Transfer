@echo off
REM Installs everything the Windows build needs into ThirdParty\ (gitignored). Run it once after cloning,
REM and again whenever a pinned version below changes. It's safe to re-run: installed versions are skipped.
REM   Scripts\Windows\SetupDependencies.bat
REM What it installs:
REM   ThirdParty\SDL3, ThirdParty\SDL3_ttf  SDL's prebuilt Visual C++ dev packages (CMakeLists.txt and
REM                                         TidyEngine.bat look for SDL here and nowhere else)
REM   ThirdParty\ShaderCross                the HLSL shader compiler MakeTransfer.bat uses. It has no prebuilt
REM                                         releases, so it's built from source: the FIRST run takes a while.
REM Needs: Visual Studio 2022 (or its Build Tools) with C++, CMake, git, Python 3. curl/tar come with Windows.
REM (macOS equivalent: Scripts/Apple/SetupDependencies.sh)
setlocal

REM Run from the repo root (two folders up from this script) no matter where the script was launched from
cd /d "%~dp0..\.."
set "REPO_ROOT=%CD%"
set "THIRD_PARTY_DIR=ThirdParty"

REM --- Pinned versions: the ONLY place they're written down ---
REM Stable SDL releases only: in SDL3's numbering an ODD minor version (3.3.x, 3.5.x) is a prerelease.
REM The SHA-256 of each .zip makes sure we unpack exactly the file that was pinned, not a corrupted
REM or tampered download. Get it from the release page's asset list, or run
REM   powershell -Command "(Get-FileHash <file>.zip).Hash"   on a download you trust.
set "SDL3_VERSION=3.4.16"
set "SDL3_SHA256=1A784CB2A5C64D56FE7A62090FE9D242D9865F235E4EA9678F1A6BA4E693E7DE"
set "SDL3_TTF_VERSION=3.2.2"
set "SDL3_TTF_SHA256=67805C5BABFC49CA0C56882DC9B8CABBCDD1E6F9EDDE10DDAC91DDB38F3AFB8C"
REM SDL_shadercross has no releases, so it's pinned to a git commit. A commit hash names its exact contents
REM (and, through the submodule pointers stored in it, the exact DirectXShaderCompiler/SPIRV sources too).
set "SHADERCROSS_COMMIT=1ff05bec573988a98ef9e0260b4da44f512b8367"

REM --- Tools ---
REM A missing tool jumps OUT of the loop to :missingTool: an "exit /b 1" inside a for loop can lose its exit code.
for %%T in (curl tar git cmake) do (
    where %%T >nul 2>&1
    if errorlevel 1 (
        set "MISSING_TOOL=%%T"
        goto :missingTool
    )
)

if not exist "%THIRD_PARTY_DIR%" mkdir "%THIRD_PARTY_DIR%"

call :installPackage SDL3 SDL %SDL3_VERSION% %SDL3_SHA256%
if errorlevel 1 exit /b 1
call :installPackage SDL3_ttf SDL_ttf %SDL3_TTF_VERSION% %SDL3_TTF_SHA256%
if errorlevel 1 exit /b 1
call :installShaderCross
if errorlevel 1 exit /b 1

echo.
echo Dependencies ready in %THIRD_PARTY_DIR%\. Next: Scripts\Windows\MakeTransfer.bat or Scripts\Windows\RunTests.bat.
exit /b 0

:missingTool
echo %MISSING_TOOL% not found. See the "Needs:" line at the top of this script.
exit /b 1


REM =====================================================================================================
REM installPackage <package> <github repo> <version> <sha256>
REM   Downloads <package>-devel-<version>-VC.zip from github.com/libsdl-org/<repo>'s release page, checks its
REM   SHA-256, and unpacks it to ThirdParty\<package>\. The folder name has no version in it, so CMake and
REM   TidyEngine.bat never need to know the version; ThirdParty\<package>\VERSION.txt records it instead.
REM   Called with "call :installPackage ..."; "exit /b" returns to the caller (not out of the script).
REM =====================================================================================================
:installPackage
set "PACKAGE=%~1"
set "REPO=%~2"
set "VERSION=%~3"
set "EXPECTED_SHA256=%~4"
set "INSTALL_DIR=%THIRD_PARTY_DIR%\%PACKAGE%"
set "ZIP_NAME=%PACKAGE%-devel-%VERSION%-VC.zip"
set "ZIP_PATH=%THIRD_PARTY_DIR%\%ZIP_NAME%"
set "URL=https://github.com/libsdl-org/%REPO%/releases/download/release-%VERSION%/%ZIP_NAME%"

REM Already installed at this exact version? Nothing to do.
set "INSTALLED_VERSION="
if exist "%INSTALL_DIR%\VERSION.txt" set /p INSTALLED_VERSION=<"%INSTALL_DIR%\VERSION.txt"
if "%INSTALLED_VERSION%"=="%VERSION%" (
    echo %PACKAGE% %VERSION% already installed.
    exit /b 0
)

echo Downloading %ZIP_NAME% ...
REM --fail: an HTTP error (e.g. 404 for a mistyped version) is an error, instead of saving the error page
REM --location: follow GitHub's redirect to the actual file
curl --fail --location --silent --show-error --output "%ZIP_PATH%" "%URL%"
if errorlevel 1 (
    echo Download failed: %URL%
    exit /b 1
)

set "ACTUAL_SHA256="
for /f "usebackq delims=" %%H in (`powershell -NoProfile -Command "(Get-FileHash -Algorithm SHA256 '%ZIP_PATH%').Hash"`) do set "ACTUAL_SHA256=%%H"
if /i not "%ACTUAL_SHA256%"=="%EXPECTED_SHA256%" (
    echo SHA-256 mismatch for %ZIP_NAME%
    echo   expected %EXPECTED_SHA256%
    echo   got      %ACTUAL_SHA256%
    del "%ZIP_PATH%"
    exit /b 1
)

REM Replace any older version. The zip holds one folder named <package>-<version>; unpack, then rename it.
if exist "%INSTALL_DIR%" rmdir /s /q "%INSTALL_DIR%"
if exist "%THIRD_PARTY_DIR%\%PACKAGE%-%VERSION%" rmdir /s /q "%THIRD_PARTY_DIR%\%PACKAGE%-%VERSION%"
tar -xf "%ZIP_PATH%" -C "%THIRD_PARTY_DIR%"
if errorlevel 1 (
    echo Unpacking %ZIP_NAME% failed.
    exit /b 1
)
ren "%THIRD_PARTY_DIR%\%PACKAGE%-%VERSION%" "%PACKAGE%"
if errorlevel 1 (
    echo Expected %ZIP_NAME% to contain a folder named %PACKAGE%-%VERSION%.
    exit /b 1
)
>"%INSTALL_DIR%\VERSION.txt" echo %VERSION%
del "%ZIP_PATH%"

echo %PACKAGE% %VERSION% installed.
exit /b 0


REM =====================================================================================================
REM installShaderCross
REM   Fetches SDL_shadercross at SHADERCROSS_COMMIT (plus its submodules), builds it with its "vendored"
REM   dependencies (DirectXShaderCompiler + SPIRV-Cross, built from source too), and installs it to
REM   ThirdParty\ShaderCross\. bin\ then holds shadercross.exe next to every DLL it loads. The multi-GB source
REM   and build folders are deleted afterwards; VERSION.txt records the commit so a re-run can skip all this.
REM
REM   Why a drive letter: DirectXShaderCompiler's build nests so deep that under a normal repo path some of its
REM   files pass Windows' 260-character path limit, and MSBuild fails. "subst T: <folder>" makes T:\ show that
REM   folder's contents, so every path gets ~60 characters shorter. The mapping is removed again whether the
REM   build works or not. (If the script is interrupted, remove it by hand: subst T: /D)
REM =====================================================================================================
:installShaderCross
set "INSTALL_DIR=%THIRD_PARTY_DIR%\ShaderCross"
set "SRC_ROOT=%REPO_ROOT%\%THIRD_PARTY_DIR%\_src"

set "INSTALLED_VERSION="
if exist "%INSTALL_DIR%\VERSION.txt" set /p INSTALLED_VERSION=<"%INSTALL_DIR%\VERSION.txt"
if "%INSTALLED_VERSION%"=="%SHADERCROSS_COMMIT%" (
    echo ShaderCross %SHADERCROSS_COMMIT% already installed.
    exit /b 0
)

REM Pick the first free drive letter from T: to Z:
set "BUILD_DRIVE="
for %%D in (T U V W X Y Z) do (
    if not defined BUILD_DRIVE if not exist %%D:\ set "BUILD_DRIVE=%%D:"
)
if not defined BUILD_DRIVE (
    echo No free drive letter between T: and Z: for the temporary build drive.
    exit /b 1
)
if not exist "%SRC_ROOT%" mkdir "%SRC_ROOT%"
subst %BUILD_DRIVE% "%SRC_ROOT%"
if errorlevel 1 (
    echo Could not map %BUILD_DRIVE% to %SRC_ROOT%.
    exit /b 1
)
echo Building on temporary drive %BUILD_DRIVE% ^(= %SRC_ROOT%^)

call :buildShaderCross %BUILD_DRIVE%\sc
set "BUILD_RESULT=%ERRORLEVEL%"

REM Delete the sources while the short drive still exists (some paths are too long to delete any other way)
if "%BUILD_RESULT%"=="0" rmdir /s /q %BUILD_DRIVE%\sc
subst %BUILD_DRIVE% /D
if not "%BUILD_RESULT%"=="0" exit /b 1
rmdir /s /q "%SRC_ROOT%"

echo ShaderCross installed.
exit /b 0


REM buildShaderCross <source folder on the temporary drive>
:buildShaderCross
set "SOURCE_DIR=%~1"

echo Downloading SDL_shadercross %SHADERCROSS_COMMIT% and its submodules ^(several hundred MB^) ...
REM A failed earlier run may have left its sources behind: start clean
if exist "%SOURCE_DIR%" rmdir /s /q "%SOURCE_DIR%"
mkdir "%SOURCE_DIR%"
REM Fetch exactly the pinned commit instead of cloning the whole history. --depth 1 = only that one snapshot.
REM core.longpaths lets git itself write the deep DirectXShaderCompiler paths.
git -C "%SOURCE_DIR%" -c core.longpaths=true init --quiet
git -C "%SOURCE_DIR%" remote add origin https://github.com/libsdl-org/SDL_shadercross.git
git -C "%SOURCE_DIR%" -c core.longpaths=true fetch --quiet --depth 1 origin %SHADERCROSS_COMMIT%
if errorlevel 1 (
    echo Fetching SDL_shadercross failed.
    exit /b 1
)
git -C "%SOURCE_DIR%" -c core.longpaths=true checkout --quiet FETCH_HEAD
git -C "%SOURCE_DIR%" -c core.longpaths=true submodule update --quiet --init --recursive --depth 1
if errorlevel 1 (
    echo Fetching SDL_shadercross's submodules failed.
    exit /b 1
)

echo Building shadercross ^(first time: expect 20+ minutes, most of it DirectXShaderCompiler^) ...
REM SDL3_DIR: build against our own ThirdParty SDL3, so no other SDL on the machine gets involved
cmake -S "%SOURCE_DIR%" -B "%SOURCE_DIR%\b" -DSDL3_DIR="%REPO_ROOT%\%THIRD_PARTY_DIR%\SDL3\cmake" -DSDLSHADERCROSS_VENDORED=ON -DSDLSHADERCROSS_INSTALL=ON -DSDLSHADERCROSS_STATIC=OFF
if errorlevel 1 (
    echo Configuring shadercross failed.
    exit /b 1
)
cmake --build "%SOURCE_DIR%\b" --config Release --parallel
if errorlevel 1 (
    echo Building shadercross failed.
    exit /b 1
)
if exist "%INSTALL_DIR%" rmdir /s /q "%INSTALL_DIR%"
cmake --install "%SOURCE_DIR%\b" --config Release --prefix "%REPO_ROOT%\%INSTALL_DIR%"
if errorlevel 1 (
    echo Installing shadercross failed.
    exit /b 1
)
REM Windows looks for DLLs next to the .exe: shadercross also needs SDL3.dll, which isn't part of its install
copy /y "%THIRD_PARTY_DIR%\SDL3\lib\x64\SDL3.dll" "%INSTALL_DIR%\bin\" >nul
>"%INSTALL_DIR%\VERSION.txt" echo %SHADERCROSS_COMMIT%
exit /b 0
