// File: Transfer/src/Core/CameraState.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/Math/Vector2.hpp"
#include "Utilities/Constants/GameSystemConstants.hpp"

// Standard Library Imports
#include <cstdint>

enum class VisorView : uint32_t
{
    Realistic = 0,
    Mass = 1,
    Charge = 2,
    Temperature = 3
};

struct CameraState
{
    double zoom = STARTUP_ZOOM_VALUE;
    DynamoEngine::Vector2D offset = {0.0, 0.0}; // pan offset
    DynamoEngine::Vector2D twinkling_star_offset = {0.0, 0.0};

    float window_width = static_cast<float>(SCREEN_WIDTH);
    float window_height = static_cast<float>(SCREEN_HEIGHT);

    float max_display_width = static_cast<float>(SCREEN_WIDTH);
    float max_display_height = static_cast<float>(SCREEN_HEIGHT);

    float render_alpha;
    VisorView visor_view = VisorView::Realistic;
};

// The view Tab switches to: Realistic -> Mass -> Charge -> Temperature -> back to Realistic
inline VisorView nextVisorView(VisorView view)
{
    switch (view)
    {
    case VisorView::Realistic:
        return VisorView::Mass;
    case VisorView::Mass:
        return VisorView::Charge;
    case VisorView::Charge:
        return VisorView::Temperature;
    case VisorView::Temperature:
        return VisorView::Realistic;
    }
    return VisorView::Realistic; // unreachable: every view is handled above
}