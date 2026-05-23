// File: Transfer/src/Utilities/System/SystemPathUtility.h

#pragma once

#include <SDL3/SDL.h>
#include <SDL3/SDL_filesystem.h>
#include <string>

namespace Utilities
{

inline std::string GetBasePath()
{
    const char* basePathC = SDL_GetBasePath();
    if (basePathC)
    {
        return std::string(basePathC); // SDL owns the memory
    }
    return "./"; // fallback
}

// Returns full path to a resource folder/file
inline std::string GetResourcePath(const std::string& relativePath)
{
    std::string base = GetBasePath();
#ifdef __APPLE__
    // If we're inside a .app bundle, base will end with Contents/Resources/
    if (base.find("Contents/Resources/") != std::string::npos)
    {
        // Release bundle: base already points to Resources
        return base + "Assets/" + relativePath;
    }
    else
    {
        // Debug build: executable outside bundle
        return base + "Resources/Assets/" + relativePath;
    }
#else
    // Windows/Linux: relative to exe folder need to implement.
    return base + "Resources/" + relativePath;
#endif
}

}; // namespace Utilities