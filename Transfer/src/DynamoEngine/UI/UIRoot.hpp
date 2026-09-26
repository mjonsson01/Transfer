// File: Transfer/src/DynamoEngine/UI/UIRoot.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/Math/Vector2.hpp"
#include "DynamoEngine/UI/UIElement.hpp"
#include "DynamoEngine/UI/UIGeometryBuilder.hpp"

// Standard Library Imports
#include <vector>

namespace DynamoEngine
{
// The top of a UI tree (one per scene). It covers the whole window and is responsible for the things
// that only make sense for the tree as a whole: UI scale, draw order, and "what is under the mouse?".
class UIRoot : public UIElement
{
  public:
    UIRoot();
    // --- Window and scale --- //

    // Call at startup and whenever the window is resized (size in screen points, e.g. 1280 x 720)
    void setWindowSize(Vector2F window_size) { m_window_size = window_size; }

    // The player's UI scale setting from the options menu (1.0 = normal, 1.5 = 50% bigger)
    void setPlayerScale(float player_scale);

    // Final scale = (window height / 720) * player scale. The UI is designed at 720p and grows with the window.
    float uiScale() const;

    // The size of the space the UI is laid out in: window size / uiScale()
    Vector2F uiSpaceSize() const;

    // Converts a mouse position (screen points) into UI space, so it can be compared with element rects
    Vector2F screenToUISpace(Vector2F screen_point) const;

    // --- Every frame --- //

    // Lays out the whole tree for the current window and scale, then calls update() on every visible element
    void updateElements(float delta_seconds);

    // Draws every visible element, back to front: layer first, then parents before children,
    // then siblings by z-index
    void drawElements(UIGeometryBuilder& builder) const;

    // --- Queries --- //

    // The element drawn on top at `ui_point` (the one a click would reach first), or nullptr if none
    UIElement* topmostElementAt(Vector2F ui_point) const;

    // Every visible element, back to front (the order drawElements uses). The root itself is not included.
    std::vector<UIElement*> elementsInDrawOrder() const;

  private:
    Vector2F m_window_size = {1280.0f, 720.0f};
    float m_player_scale = 1.0f;
};
} // namespace DynamoEngine