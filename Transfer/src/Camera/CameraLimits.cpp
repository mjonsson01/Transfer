// File: Transfer/src/Camera/CameraLimits.cpp
#include "Camera/CameraLimits.hpp"

// Custom Imports
#include "Constants/GameSystemConstants.hpp"

// Standard Library Imports
#include <algorithm>

// Keeps the view's center where the view stays inside a rectangle (center +- half size, world units). When the view
// is bigger than the rectangle, the view is centered on it.
static void keepViewInsideRectangle(DynamoEngine::Vector2D& offset_to_limit, DynamoEngine::Vector2D rectangle_center,
                                    double rectangle_half_width, double rectangle_half_height, double view_half_width,
                                    double view_half_height)
{
    double slack_x = std::max(0.0, rectangle_half_width - view_half_width);
    double slack_y = std::max(0.0, rectangle_half_height - view_half_height);

    DynamoEngine::Vector2D view_center_world = {view_half_width - offset_to_limit.x_val,
                                                view_half_height - offset_to_limit.y_val};
    view_center_world.x_val =
        std::clamp(view_center_world.x_val, rectangle_center.x_val - slack_x, rectangle_center.x_val + slack_x);
    view_center_world.y_val =
        std::clamp(view_center_world.y_val, rectangle_center.y_val - slack_y, rectangle_center.y_val + slack_y);

    offset_to_limit.x_val = view_half_width - view_center_world.x_val;
    offset_to_limit.y_val = view_half_height - view_center_world.y_val;
}

DynamoEngine::Vector2D keepViewInsideWorld(CameraState& camera_state, const DynamoEngine::Vector2D& world_min_corner,
                                           const DynamoEngine::Vector2D& world_max_corner)
{
    double view_half_width = (camera_state.window_width / 2.0) / camera_state.zoom;
    double view_half_height = (camera_state.window_height / 2.0) / camera_state.zoom;

    const DynamoEngine::Vector2D offset_before = camera_state.offset;
    const DynamoEngine::Vector2D world_center = (world_min_corner + world_max_corner) / 2.0;
    const DynamoEngine::Vector2D world_half_size = (world_max_corner - world_min_corner) / 2.0;
    keepViewInsideRectangle(camera_state.offset, world_center, world_half_size.x_val, world_half_size.y_val,
                         view_half_width, view_half_height);

    double star_field_half_width = camera_state.max_display_width / (2.0 * MIN_ZOOM);
    double star_field_half_height = camera_state.max_display_height / (2.0 * MIN_ZOOM);
    DynamoEngine::Vector2D star_field_center = {SCREEN_WIDTH / 2.0, SCREEN_HEIGHT / 2.0};
    keepViewInsideRectangle(camera_state.twinkling_star_offset, star_field_center, star_field_half_width,
                         star_field_half_height, view_half_width, view_half_height);

    return camera_state.offset - offset_before;
}
