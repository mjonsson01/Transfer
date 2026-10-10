// File: Transfer/src/DynamoEngine/UI/UISound.hpp

#pragma once

// Standard Library Imports
#include <cstdint>

namespace DynamoEngine
{
enum class UISound : uint8_t
{
    SpecialClick,  // a special button was clicked
    Tick,          // a slider moved past one of its tick marks
    StandardClick, // a standard button was clicked
};
} // namespace DynamoEngine
