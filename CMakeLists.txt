cmake_minimum_required(VERSION 3.20)
project(TransferGame LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# set(CMAKE_BUILD_TYPE Debug)
# ----------------------------------
# Platform-specific SDL setup
# ----------------------------------
if (WIN32)
    # --- Windows paths --- NEEDS TO BE UPDATED.
    set(SDL3_DIR      "${CMAKE_SOURCE_DIR}/../../SDL3")
    set(SDL3_TTF_DIR  "${CMAKE_SOURCE_DIR}/../../SDL3_TTF")

    set(SDL3_INCLUDE_DIRS
        "${SDL3_DIR}/include"
        "${SDL3_TTF_DIR}/include"
    )

    set(SDL3_LIBS
        "${SDL3_DIR}/lib/x64/SDL3.lib"
        "${SDL3_TTF_DIR}/lib/x64/SDL3_ttf.lib"
    )

elseif (APPLE)
    # --- macOS Homebrew installation paths ---
    message(STATUS "Configuring for macOS (Homebrew SDL3 and SDL3_ttf expected)")

    # Try to find SDL3 and SDL3_ttf automatically
    # find_package(SDL3 REQUIRED CONFIG)
    # find_package(SDL3_ttf REQUIRED CONFIG)

    # If not found automatically, you can override manually like:
    set(SDL3_DIR "/opt/homebrew/lib/cmake/SDL3")
    find_package(SDL3 REQUIRED CONFIG)

    set(SDL3_TTF_DIR "/opt/homebrew/lib/cmake/SDL3_ttf")
    find_package(SDL3_ttf REQUIRED CONFIG)
    # Note: Homebrew installs headers and libraries under /opt/homebrew
    
    set(SDL3_LIBS
        SDL3::SDL3
        SDL3_ttf::SDL3_ttf
    )
else()
    message(FATAL_ERROR "Unsupported platform. Only Windows and macOS are supported right now.")
endif()

# ----------------------------------
# Collect source files
# ----------------------------------
file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/Transfer/src/*.cpp"
    "${CMAKE_SOURCE_DIR}/Transfer/src/*.h"
)

# ----------------------------------
# Executable
# ----------------------------------
# add_executable(${PROJECT_NAME} ${SRC_FILES})
add_executable(${PROJECT_NAME} ${SRC_FILES} Transfer.rc)

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
        ${SDL3_LIBS}
)

# ----------------------------------
# Copy DLLs to output directory (Windows only)
# ----------------------------------
if (WIN32)
    foreach(dll IN ITEMS
        "${SDL3_DIR}/lib/x64/SDL3.dll"
        "${SDL3_TTF_DIR}/lib/x64/SDL3_ttf.dll"
    )
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${dll}" $<TARGET_FILE_DIR:${PROJECT_NAME}>
        )
    endforeach()
endif()

# --- Copy Assets folder next to the executable ---
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/Transfer/Assets"
        "$<TARGET_FILE_DIR:${PROJECT_NAME}>/Assets"
)

# --- Custom targets ---
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
