cmake_minimum_required(VERSION 3.20)
project(TransferGame LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ----------------------------------
# Paths to SDL dependencies
# (these could be turned into cache vars if you want to override them)
# ----------------------------------
set(SDL3_DIR      "${CMAKE_SOURCE_DIR}/../../SDL3")
set(SDL3_TTF_DIR  "${CMAKE_SOURCE_DIR}/../../SDL3_TTF")

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
add_executable(${PROJECT_NAME} ${SRC_FILES})

# ----------------------------------
# Include directories
# ----------------------------------
target_include_directories(${PROJECT_NAME}
    PRIVATE
        "${CMAKE_SOURCE_DIR}/Transfer/src"
        "${SDL3_DIR}/include"
        "${SDL3_TTF_DIR}/include"
)

# ----------------------------------
# Link libraries
# ----------------------------------
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        "${SDL3_DIR}/lib/x64/SDL3.lib"
        "${SDL3_TTF_DIR}/lib/x64/SDL3_ttf.lib"
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


# --- Copy Assets folder next to the executable ---
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/Transfer/Assets"
        "$<TARGET_FILE_DIR:${PROJECT_NAME}>/Assets"
)