// File: Transfer/src/DynamoEngine/UI/UIRoot.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/Input/InputState.hpp"
#include "DynamoEngine/Math/Vector2.hpp"
#include "DynamoEngine/UI/UIElement.hpp"
#include "DynamoEngine/UI/UIGeometryBuilder.hpp"
#include "DynamoEngine/UI/UIInputResult.hpp"

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

    // --- Mouse input --- //

    // Delivers this frame's mouse input to the elements (hover, press, drag, release) and reports what the UI took.
    // Call once per frame, after updateElements() so every element is in its current place.
    UIInputResult processInput(const InputState& input);

    // The element drawn on top at `ui_point` (the one a click would reach first), or nullptr if none
    UIElement* topmostElementAt(Vector2F ui_point) const;

    // Every visible element, back to front (the order drawElements uses). The root itself is not included.
    std::vector<UIElement*> elementsInDrawOrder() const;

  private:
    Vector2F m_window_size = {1280.0f, 720.0f};
    float m_player_scale = 1.0f;
    void updateHover(Vector2F ui_mouse_position);
    void pressTopmostElementThatWantsIt(Vector2F ui_mouse_position);
    // Non-owning: these point INTO the tree. Safe because elements are only destroyed along with the root.
    // (If a removeChild() is ever added, it must clear these when it removes one of them.)
    UIElement* m_hovered_element = nullptr;  // the element under the cursor, or nullptr
    UIElement* m_captured_element = nullptr; // the element that took the current press, or nullptr
};
} // namespace DynamoEngine