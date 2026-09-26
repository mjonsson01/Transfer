// File: Transfer/src/DynamoEngine/UI/UIRoot.cpp

#include "DynamoEngine/UI/UIRoot.hpp"

// Custom Imports
#include "DynamoEngine/Constants/GlobalConstants.hpp"

// Standard Library Imports
#include <algorithm>
#include <memory>

namespace
{
using DynamoEngine::UIElement;

// Sorting rules, written as named functions so the sort calls below read as sentences

bool hasLowerZIndex(const UIElement* first, const UIElement* second) { return first->zIndex() < second->zIndex(); }

bool isInLowerLayer(const UIElement* first, const UIElement* second) { return first->layer() < second->layer(); }

// An element's children, ordered by z-index. stable_sort: equal z-indexes keep the order they were added in.
std::vector<UIElement*> childrenByZIndex(const UIElement& element)
{
    std::vector<UIElement*> children;
    for (const std::unique_ptr<UIElement>& child : element.children())
    {
        children.push_back(child.get());
    }
    std::stable_sort(children.begin(), children.end(), hasLowerZIndex);
    return children;
}

// Adds `element` and its visible descendants to `ordered_elements`: the element first, then each of its
// children (by z-index) followed by that child's own children. A hidden element hides its whole subtree.
void collectInTreeOrder(UIElement& element, std::vector<UIElement*>& ordered_elements)
{
    if (!element.isVisible())
    {
        return;
    }
    ordered_elements.push_back(&element);

    for (UIElement* child : childrenByZIndex(element))
    {
        collectInTreeOrder(*child, ordered_elements);
    }
}
} // namespace

namespace DynamoEngine
{
UIRoot::UIRoot()
{
    setPlacement({.align = UIAlign::Fill}); // the root always covers the whole UI space
}

void UIRoot::setPlayerScale(float player_scale)
{
    if (player_scale > 0.0f) // ignore nonsense values instead of collapsing the UI to nothing
    {
        m_player_scale = player_scale;
    }
}

float UIRoot::uiScale() const
{
    if (m_window_size.y_val <= 0.0f) // minimized window: keep a usable scale instead of dividing by zero later
    {
        return m_player_scale;
    }
    return (m_window_size.y_val / UI_REFERENCE_HEIGHT) * m_player_scale;
}

Vector2F UIRoot::uiSpaceSize() const { return m_window_size / uiScale(); }

Vector2F UIRoot::screenToUISpace(Vector2F screen_point) const { return screen_point / uiScale(); }

void UIRoot::updateElements(float delta_seconds)
{
    // Layout runs every frame so elements added, moved, or resized since last frame are always in the right place
    const Vector2F space = uiSpaceSize();
    updateLayout(SDL_FRect{0.0f, 0.0f, space.x_val, space.y_val});

    for (UIElement* element : elementsInDrawOrder())
    {
        element->update(delta_seconds);
    }
}

void UIRoot::drawElements(UIGeometryBuilder& builder) const
{
    for (const UIElement* element : elementsInDrawOrder())
    {
        element->draw(builder);
    }
}
UIInputResult UIRoot::processInput(const InputState& input)
{
    const Vector2F mouse_position = screenToUISpace(input.mousePosition());
    UIInputResult result;

    // 1. Hover: tell elements when the cursor moves onto or off them
    updateHover(mouse_position);

    // 2. Press: offer the click to elements under the cursor, top to bottom, until one takes it
    if (input.wasMouseButtonPressed(MouseButton::Left) && m_captured_element == nullptr)
    {
        pressTopmostElementThatWantsIt(mouse_position);
    }

    // 3. While a press is captured, the element that took it gets every drag and the release
    if (m_captured_element != nullptr)
    {
        result.pointer_captured = true; // stays true on the release frame too, so the game never sees that release

        if (input.isMouseButtonDown(MouseButton::Left))
        {
            m_captured_element->onMouseDragged(mouse_position);
        }
        else // released this frame (or the window lost focus, which also lets go of the button)
        {
            const bool released_inside = m_captured_element->containsPoint(mouse_position);
            m_captured_element->onMouseReleased(mouse_position, released_inside);
            m_captured_element = nullptr;
        }
    }

    result.pointer_over_ui = (m_hovered_element != nullptr);
    return result;
}

void UIRoot::updateHover(Vector2F ui_mouse_position)
{
    UIElement* element_under_mouse = topmostElementAt(ui_mouse_position);
    if (element_under_mouse == m_hovered_element)
    {
        return; // still over the same element (or still over nothing)
    }

    if (m_hovered_element != nullptr)
    {
        m_hovered_element->onMouseExited();
    }
    if (element_under_mouse != nullptr)
    {
        element_under_mouse->onMouseEntered();
    }
    m_hovered_element = element_under_mouse;
}

void UIRoot::pressTopmostElementThatWantsIt(Vector2F ui_mouse_position)
{
    const std::vector<UIElement*> draw_order = elementsInDrawOrder();

    // Front to back: the top element gets the first chance; if it says no (returns false), the next one down does
    for (auto element = draw_order.rbegin(); element != draw_order.rend(); ++element)
    {
        if ((*element)->containsPoint(ui_mouse_position) && (*element)->onMousePressed(ui_mouse_position))
        {
            m_captured_element = *element;
            return;
        }
    }
}

UIElement* UIRoot::topmostElementAt(Vector2F ui_point) const
{
    const std::vector<UIElement*> draw_order = elementsInDrawOrder();

    // Walk the draw order backwards (front to back): the first element containing the point is the top one
    for (auto element = draw_order.rbegin(); element != draw_order.rend(); ++element)
    {
        if ((*element)->containsPoint(ui_point))
        {
            return *element;
        }
    }
    return nullptr;
}

std::vector<UIElement*> UIRoot::elementsInDrawOrder() const
{
    // 1. Tree order: parents before children, siblings by z-index
    std::vector<UIElement*> draw_order;
    for (UIElement* child : childrenByZIndex(*this))
    {
        collectInTreeOrder(*child, draw_order);
    }

    // 2. Group by layer. stable_sort keeps the tree order inside each layer.
    std::stable_sort(draw_order.begin(), draw_order.end(), isInLowerLayer);
    return draw_order;
}
} // namespace DynamoEngine