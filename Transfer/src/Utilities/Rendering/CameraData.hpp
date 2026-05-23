// File: Transfer/src/Utilities/Rendering/CameraData.hpp

#pragma once

struct CameraConstants
{
    float screenWidth;
    float screenHeight;
    float zoom;
    float offsetX;
    float offsetY;
    uint32_t viewMode; // Real View, Charge View, Mass View, Magnetism View, Temp View, etc.
    float _padding[2]; // fills in last 8 bytes to make a 32 bit struct.
};