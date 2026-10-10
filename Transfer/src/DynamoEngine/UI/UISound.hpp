// File: Transfer/src/DynamoEngine/UI/UISound.hpp

#pragma once

// Standard Library Imports
#include <cstdint>

namespace DynamoEngine
{
enum class UISound : uint8_t
{
    Click,    // a button was clicked
    Tick,     // a slider moved past one of its tick marks
    Checkbox, // a checkbox was toggled
};
} // namespace DynamoEngine
