cmake_minimum_required(VERSION 3.20)
project(TransferGame LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ----------------------------------
# Collect source files
# ----------------------------------
file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/Transfer/src/*.cpp"
    "${CMAKE_SOURCE_DIR}/Transfer/src/*.hpp"
)

# ----------------------------------
# Platform-specific configuration
# ----------------------------------
if(WIN32)
    message(STATUS "Configuring for Windows")

    # SDL paths
    set(SDL3_DIR "C:/SDL3")
    set(SDL3_TTF_DIR "${CMAKE_SOURCE_DIR}/../../SDL3_TTF")

    set(SDL3_INCLUDE_DIRS
        "${SDL3_DIR}/include"
        "${SDL3_TTF_DIR}/include"
    )

    set(SDL3_LIBS
        "${SDL3_DIR}/lib/SDL3.lib"
        "${SDL3_TTF_DIR}/lib/x64/SDL3_ttf.lib"
    )

    # DLLs to copy
    set(SDL3_DLLS
        "${SDL3_DIR}/bin/SDL3.dll"
        "${SDL3_TTF_DIR}/lib/x64/SDL3_ttf.dll"
    )
    list(APPEND SRC_FILES "${CMAKE_SOURCE_DIR}/Transfer/Assets/Transfer.rc")

elseif(APPLE)
    message(STATUS "Configuring for macOS")

    # Find SDL via Homebrew
    set(SDL3_DIR "/opt/homebrew/lib/cmake/SDL3")
    find_package(SDL3 REQUIRED CONFIG)

    set(SDL3_TTF_DIR "/opt/homebrew/lib/cmake/SDL3_ttf")
    find_package(SDL3_ttf REQUIRED CONFIG)

    set(SDL3_LIBS
        SDL3::SDL3
        SDL3_ttf::SDL3_ttf
    )

    # macOS bundle info
    set(MACOSX_BUNDLE_ICON_FILE Transfer.icns)
    set(MACOSX_BUNDLE_INFO_PLIST "${CMAKE_SOURCE_DIR}/Transfer/Info.plist")
endif()

# ----------------------------------
# DynamoEngine (static library) -- reusable engine code; must never include game code
# ----------------------------------
file(GLOB_RECURSE DYNAMO_ENGINE_SRC_FILES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/Transfer/src/DynamoEngine/*.cpp"
    "${CMAKE_SOURCE_DIR}/Transfer/src/DynamoEngine/*.hpp"
)
list(FILTER SRC_FILES EXCLUDE REGEX "/DynamoEngine/") # the game links the engine instead of compiling it itself

add_library(DynamoEngine STATIC ${DYNAMO_ENGINE_SRC_FILES})
target_include_directories(DynamoEngine PUBLIC "${CMAKE_SOURCE_DIR}/Transfer/src" ${SDL3_INCLUDE_DIRS})
target_link_libraries(DynamoEngine PUBLIC ${SDL3_LIBS})
target_compile_options(DynamoEngine PRIVATE
    $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-Wall -Wextra -Wimplicit-fallthrough>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

# ----------------------------------
# Executable definition
# ----------------------------------
if(APPLE AND CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Building Release as .app bundle")
    add_executable(${PROJECT_NAME} MACOSX_BUNDLE ${SRC_FILES})
    set_target_properties(${PROJECT_NAME} PROPERTIES
        MACOSX_BUNDLE TRUE
        MACOSX_BUNDLE_ICON_FILE Transfer.icns
    )
else()
    add_executable(${PROJECT_NAME} ${SRC_FILES})
endif()

# ----------------------------------
# Include directories
# ----------------------------------
target_include_directories(${PROJECT_NAME}
    PRIVATE
        "${CMAKE_SOURCE_DIR}/Transfer/src"
        ${SDL3_INCLUDE_DIRS}
)

# ----------------------------------
# Link libraries
# ----------------------------------
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        DynamoEngine
        ${SDL3_LIBS}
)

# ----------------------------------
# Post-build steps
# ----------------------------------
if(WIN32)
    # Copy DLLs
    foreach(dll IN LISTS SDL3_DLLS)
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${dll}" $<TARGET_FILE_DIR:${PROJECT_NAME}>
        )
    endforeach()

    # Copy ICO icon next to executable ONLY
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_SOURCE_DIR}/Transfer/Assets/Visual/Transfer.ico"
            "$<TARGET_FILE_DIR:${PROJECT_NAME}>/Transfer.ico"
        COMMENT "Copying app icon (ICO) next to executable for Windows"
    )

    # Copy other assets
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_SOURCE_DIR}/Transfer/Assets"
            "$<TARGET_FILE_DIR:${PROJECT_NAME}>/Resources"
        COMMENT "Copying Resources next to executable"
    )
elseif(APPLE)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        # Copy icon into bundle
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
                "${CMAKE_SOURCE_DIR}/Transfer/Assets/Visual/Transfer.icns"
                "$<TARGET_BUNDLE_DIR:${PROJECT_NAME}>/Contents/Resources/Transfer.icns"
            COMMENT "Copying app icon into bundle Resources folder"
        )

        # Copy all assets into app bundle
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_directory
                "${CMAKE_SOURCE_DIR}/Transfer/Assets"
                "$<TARGET_BUNDLE_DIR:${PROJECT_NAME}>/Contents/Resources/Assets/"
            COMMENT "Copying game resources into app bundle"
        )
    else()
        # Debug: simplified asset copy
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_directory
                "${CMAKE_SOURCE_DIR}/Transfer/Assets"
                "$<TARGET_FILE_DIR:${PROJECT_NAME}>/Resources/Assets/"
            COMMENT "Copying game resources into Resources/Assets for Debug"
        )
    endif()
endif()

# ----------------------------------
# Custom targets
# ----------------------------------
add_custom_target(clean_all
    COMMAND ${CMAKE_COMMAND} -E echo "Cleaning build and binary files ..."
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} -E echo "Clean complete."
    COMMENT "Remove all build artifacts and CMake cache."
)

add_custom_target(buildTransfer
    COMMAND ${CMAKE_COMMAND} -S "${CMAKE_SOURCE_DIR}" -B "${CMAKE_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} --build "${CMAKE_BINARY_DIR}" --target ${PROJECT_NAME}
    COMMENT "Building TransferGame target."
)


# ----------------------------------
# Unit tests (googletest) -- off for game builds; RunTests.sh turns them on
# ----------------------------------
option(TRANSFER_BUILD_TESTS "Build the TransferTests unit test executable" OFF)

if(TRANSFER_BUILD_TESTS)
    if(POLICY CMP0135)
        cmake_policy(SET CMP0135 NEW) # use extraction time as file timestamps (silences a FetchContent warning)
    endif()

    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.tar.gz
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE) # Windows: same C runtime as the game
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)         # we never install googletest
    FetchContent_MakeAvailable(googletest)

    enable_testing()

    # Tests + every engine source. Engine code must build without any game code -- if this target
    # fails to link, something in DynamoEngine/ reached into the game.
    file(GLOB_RECURSE TEST_SRC_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/Tests/*.cpp")
    file(GLOB_RECURSE ENGINE_SRC_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/Transfer/src/DynamoEngine/*.cpp")

    if(TEST_SRC_FILES)
        add_executable(TransferTests ${TEST_SRC_FILES})
        target_include_directories(TransferTests PRIVATE
            "${CMAKE_SOURCE_DIR}/Transfer/src"  # engine / game headers:  #include "DynamoEngine/..."
            "${CMAKE_SOURCE_DIR}/Tests"         # test utilities:         #include "TestingUtilities/..."
            ${SDL3_INCLUDE_DIRS}
        )
        target_link_libraries(TransferTests PRIVATE GTest::gtest_main DynamoEngine)
        target_link_libraries(TransferTests PRIVATE GTest::gtest_main ${SDL3_LIBS})

        include(GoogleTest)
        gtest_discover_tests(TransferTests)
    else()
        message(STATUS "TransferTests: no test sources in Tests yet -- skipping the test target")
    endif()
endif()

# ----------------------------------
# VS Code IntelliSense config -- generates .vscode/c_cpp_properties.json on every configure
# ----------------------------------
# Every build script (MakeTransfer.sh/.bat, RunTest.sh) runs a configure, so they all refresh this file.
# It points IntelliSense at each build folder's compile_commands.json (game build/ + tests build-tests/),
# so game sources and test sources both get their exact compile flags. .vscode/ is gitignored and
# machine-specific: hand edits to this file will be overwritten. Turn off with -DTRANSFER_GENERATE_VSCODE_CONFIG=OFF.
option(TRANSFER_GENERATE_VSCODE_CONFIG "Write .vscode/c_cpp_properties.json on configure" ON)

if(TRANSFER_GENERATE_VSCODE_CONFIG)
    # compile_commands.json from the build folder being configured right now (CMake writes it at the end of
    # this run, so it may not exist yet) plus any other known build folder that has been built before.
    set(VSCODE_COMPILE_COMMANDS "")
    foreach(build_dir IN ITEMS "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/build" "${CMAKE_SOURCE_DIR}/build-tests")
        if(build_dir STREQUAL CMAKE_BINARY_DIR OR EXISTS "${build_dir}/compile_commands.json")
            file(RELATIVE_PATH relative_build_dir "${CMAKE_SOURCE_DIR}" "${build_dir}")
            list(APPEND VSCODE_COMPILE_COMMANDS "\"\${workspaceFolder}/${relative_build_dir}/compile_commands.json\"")
        endif()
    endforeach()
    list(REMOVE_DUPLICATES VSCODE_COMPILE_COMMANDS)
    list(JOIN VSCODE_COMPILE_COMMANDS ",\n                " VSCODE_COMPILE_COMMANDS_JSON)

    # Fallback include paths, used for files that aren't in any compile_commands.json (e.g. a header opened on its own)
    set(VSCODE_INCLUDE_PATHS "\"\${workspaceFolder}/Transfer/src\"" "\"\${workspaceFolder}/Tests\"")
    if(SDL3_INCLUDE_DIRS)
        foreach(include_dir IN LISTS SDL3_INCLUDE_DIRS) # Windows: explicit SDL folders
            list(APPEND VSCODE_INCLUDE_PATHS "\"${include_dir}\"")
        endforeach()
    elseif(APPLE)
        list(APPEND VSCODE_INCLUDE_PATHS "\"/opt/homebrew/include\"") # Homebrew SDL3 / SDL3_ttf
    endif()
    if(EXISTS "${CMAKE_SOURCE_DIR}/build-tests/_deps/googletest-src/googletest/include")
        list(APPEND VSCODE_INCLUDE_PATHS "\"\${workspaceFolder}/build-tests/_deps/googletest-src/googletest/include\"")
    endif()
    list(JOIN VSCODE_INCLUDE_PATHS ",\n                " VSCODE_INCLUDE_PATHS_JSON)

    # IntelliSense mode = <platform>-<compiler>-<architecture>
    string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" VSCODE_ARCH)
    if(VSCODE_ARCH MATCHES "arm64|aarch64")
        set(VSCODE_ARCH "arm64")
    else()
        set(VSCODE_ARCH "x64")
    endif()

    if(MSVC)
        set(VSCODE_CONFIG_NAME "Win32")
        set(VSCODE_INTELLISENSE_MODE "windows-msvc-${VSCODE_ARCH}")
    elseif(APPLE)
        set(VSCODE_CONFIG_NAME "Mac")
        set(VSCODE_INTELLISENSE_MODE "macos-clang-${VSCODE_ARCH}")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(VSCODE_CONFIG_NAME "Linux")
        set(VSCODE_INTELLISENSE_MODE "linux-clang-${VSCODE_ARCH}")
    else()
        set(VSCODE_CONFIG_NAME "Linux")
        set(VSCODE_INTELLISENSE_MODE "linux-gcc-${VSCODE_ARCH}")
    endif()

    # [=[ ... ]=] is a literal block (no ${} expansion); @ONLY fills in just the @VARIABLES@.
    # file(CONFIGURE) only rewrites the file when its contents actually change.
    file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/.vscode")
    file(CONFIGURE OUTPUT "${CMAKE_SOURCE_DIR}/.vscode/c_cpp_properties.json" CONTENT [=[
{
    "configurations": [
        {
            "name": "@VSCODE_CONFIG_NAME@",
            "compileCommands": [
                @VSCODE_COMPILE_COMMANDS_JSON@
            ],
            "includePath": [
                @VSCODE_INCLUDE_PATHS_JSON@
            ],
            "compilerPath": "@CMAKE_CXX_COMPILER@",
            "cStandard": "c17",
            "cppStandard": "c++@CMAKE_CXX_STANDARD@",
            "intelliSenseMode": "@VSCODE_INTELLISENSE_MODE@"
        }
    ],
    "version": 4
}
]=] @ONLY)
    message(STATUS "Updated .vscode/c_cpp_properties.json")
endif()