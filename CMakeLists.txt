# cmake_minimum_required(VERSION 3.20)
# project(TransferGame LANGUAGES CXX)

# set(CMAKE_CXX_STANDARD 20)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)
# set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# # ----------------------------------
# # Platform-specific SDL setup
# # ----------------------------------
# if (WIN32)
#     set(SDL3_DIR      "${CMAKE_SOURCE_DIR}/../../SDL3")
#     set(SDL3_TTF_DIR  "${CMAKE_SOURCE_DIR}/../../SDL3_TTF")

#     set(SDL3_INCLUDE_DIRS
#         "${SDL3_DIR}/include"
#         "${SDL3_TTF_DIR}/include"
#     )

#     set(SDL3_LIBS
#         "${SDL3_DIR}/lib/x64/SDL3.lib"
#         "${SDL3_TTF_DIR}/lib/x64/SDL3_ttf.lib"
#     )

# elseif (APPLE)
#     message(STATUS "Configuring for macOS (Homebrew SDL3 and SDL3_ttf expected)")

#     set(SDL3_DIR "/opt/homebrew/lib/cmake/SDL3")
#     find_package(SDL3 REQUIRED CONFIG)

#     set(SDL3_TTF_DIR "/opt/homebrew/lib/cmake/SDL3_ttf")
#     find_package(SDL3_ttf REQUIRED CONFIG)

#     set(SDL3_LIBS
#         SDL3::SDL3
#         SDL3_ttf::SDL3_ttf
#     )

#     # macOS bundle info
#     set(MACOSX_BUNDLE_ICON_FILE Transfer.icns)
#     set(MACOSX_BUNDLE_INFO_PLIST "${CMAKE_SOURCE_DIR}/Transfer/Info.plist")
# endif()

# # ----------------------------------
# # Collect source files
# # ----------------------------------
# file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS
#     "${CMAKE_SOURCE_DIR}/Transfer/src/*.cpp"
#     "${CMAKE_SOURCE_DIR}/Transfer/src/*.h"
# )

# # ----------------------------------
# # Executable and macOS bundle setup
# # ----------------------------------
# if(APPLE)
#     if(CMAKE_BUILD_TYPE STREQUAL "Release")
#         message(STATUS "Building Release as .app bundle")

#         add_executable(${PROJECT_NAME} MACOSX_BUNDLE ${SRC_FILES})

#         set_target_properties(${PROJECT_NAME} PROPERTIES
#             MACOSX_BUNDLE TRUE
#             MACOSX_BUNDLE_ICON_FILE Transfer.icns
#         )

#         add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
#             COMMAND ${CMAKE_COMMAND} -E copy
#                 "${CMAKE_SOURCE_DIR}/Transfer/Assets/Visual/Transfer.icns"
#                 "$<TARGET_BUNDLE_DIR:${PROJECT_NAME}>/Contents/Resources/Transfer.icns"
#             COMMENT "Copying app icon into bundle Resources folder"
#         )

#         add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
#             COMMAND ${CMAKE_COMMAND} -E copy_directory
#                 "${CMAKE_SOURCE_DIR}/Transfer/Assets"
#                 "$<TARGET_BUNDLE_DIR:${PROJECT_NAME}>/Contents/Resources/Assets/"
#             COMMENT "Copying game resources into app bundle"
#         )

#     else()
#         message(STATUS "Building Debug as normal executable")
#         add_executable(${PROJECT_NAME} ${SRC_FILES})

#         # --- Simplified Debug copy: only copy to Resources/Assets ---
#         add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
#             COMMAND ${CMAKE_COMMAND} -E copy_directory
#                 "${CMAKE_SOURCE_DIR}/Transfer/Assets"
#                 "$<TARGET_FILE_DIR:${PROJECT_NAME}>/Resources/Assets/"
#             COMMENT "Copying game resources into Resources/Assets for Debug"
#         )
#     endif()
# else()
#     # Windows/Linux: normal executable
#     add_executable(${PROJECT_NAME} ${SRC_FILES})
# endif()

# # ----------------------------------
# # Include directories
# # ----------------------------------
# target_include_directories(${PROJECT_NAME}
#     PRIVATE
#         "${CMAKE_SOURCE_DIR}/Transfer/src"
#         ${SDL3_INCLUDE_DIRS}
# )

# # ----------------------------------
# # Link libraries
# # ----------------------------------
# target_link_libraries(${PROJECT_NAME}
#     PRIVATE
#         ${SDL3_LIBS}
# )

# # ----------------------------------
# # Copy DLLs to output directory (Windows only)
# # ----------------------------------
# if (WIN32)
#     foreach(dll IN ITEMS
#         "${SDL3_DIR}/lib/x64/SDL3.dll"
#         "${SDL3_TTF_DIR}/lib/x64/SDL3_ttf.dll"
#     )
#         add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
#             COMMAND ${CMAKE_COMMAND} -E copy_if_different
#                 "${dll}" $<TARGET_FILE_DIR:${PROJECT_NAME}>
#         )
#     endforeach()

#     add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E copy_if_different
#             "${CMAKE_SOURCE_DIR}/Transfer/Assets/Visual/Transfer.ico"
#             "$<TARGET_FILE_DIR:${PROJECT_NAME}>/Transfer.ico"
#         COMMENT "Copying app icon (ICO) next to executable for Windows"
#     )

#     # Copy resources next to exe
#     add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E copy_directory
#             "${CMAKE_SOURCE_DIR}/Transfer/Assets"
#             "$<TARGET_FILE_DIR:${PROJECT_NAME}>/Resources"
#         COMMENT "Copying Resources next to executable"
#     )
# endif()

# # ----------------------------------
# # Custom targets
# # ----------------------------------
# add_custom_target(clean_all
#     COMMAND ${CMAKE_COMMAND} -E echo "Cleaning build and binary files ..."
#     COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}"
#     COMMAND ${CMAKE_COMMAND} -E echo "Clean complete."
#     COMMENT "Remove all build artifacts and CMake cache."
# )

# add_custom_target(buildTransfer
#     COMMAND ${CMAKE_COMMAND} -S "${CMAKE_SOURCE_DIR}" -B "${CMAKE_BINARY_DIR}"
#     COMMAND ${CMAKE_COMMAND} --build "${CMAKE_BINARY_DIR}" --target ${PROJECT_NAME}
#     COMMENT "Building TransferGame target."
# )

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
    "${CMAKE_SOURCE_DIR}/Transfer/src/*.h"
)

# ----------------------------------
# Platform-specific configuration
# ----------------------------------
if(WIN32)
    message(STATUS "Configuring for Windows")

    # SDL paths
    set(SDL3_DIR "${CMAKE_SOURCE_DIR}/../../SDL3")
    set(SDL3_TTF_DIR "${CMAKE_SOURCE_DIR}/../../SDL3_TTF")

    set(SDL3_INCLUDE_DIRS
        "${SDL3_DIR}/include"
        "${SDL3_TTF_DIR}/include"
    )

    set(SDL3_LIBS
        "${SDL3_DIR}/lib/x64/SDL3.lib"
        "${SDL3_TTF_DIR}/lib/x64/SDL3_ttf.lib"
    )

    # DLLs to copy
    set(SDL3_DLLS
        "${SDL3_DIR}/lib/x64/SDL3.dll"
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