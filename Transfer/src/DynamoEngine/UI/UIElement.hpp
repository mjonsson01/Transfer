// File: Transfer/src/DynamoEngine/UI/UIElement.hpp

#pragma once

// SDL Imports
#include <SDL3/SDL_rect.h>

// Custom Imports
#include "DynamoEngine/Math/Vector2.hpp"
#include "DynamoEngine/UI/UIGeometryBuilder.hpp"
#include "DynamoEngine/UI/UIPlacement.hpp"

// Standard Library Imports
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

namespace DynamoEngine
{
// Draw-order bands: everything in a higher layer draws above (and is clicked before) lower layers
enum class UILayer : uint8_t
{
    HUD,     // always-on overlay: sliders, FPS counter
    Menu,    // menus and panels
    Overlay, // open dropdown lists, popups
};

// Base class for everything in the UI. Each element owns its children.
// Subclasses override draw() and whichever mouse events they care about.
class UIElement
{
  public:
    virtual ~UIElement() = default;

    // --- Children --- //
    UIElement& addChild(std::unique_ptr<UIElement> child); // takes ownership
    const std::vector<std::unique_ptr<UIElement>>& children() const { return m_children; }
    UIElement* parent() const { return m_parent; }

    // --- Position and size --- //
    void setPlacement(const UIPlacement& placement) { m_placement = placement; }
    virtual void updateLayout(const SDL_FRect& parent_rect); // sets rect(), then updates the children
    const SDL_FRect& rect() const { return m_rect; }
    virtual bool containsPoint(Vector2F point) const;

    // --- Visibility and draw order --- //
    void setVisible(bool is_visible) { m_is_visible = is_visible; }
    bool isVisible() const { return m_is_visible; }
    void setLayer(UILayer layer) { m_layer = layer; }
    UILayer layer() const;                               // if never set: the parent's layer
    void setZIndex(int z_index) { m_z_index = z_index; } // order among siblings, higher on top
    int zIndex() const { return m_z_index; }

    // --- Override these --- //
    virtual void update(float delta_seconds) {}
    virtual void draw(UIGeometryBuilder& builder) const {}                 // this element only, not children
    virtual bool onMousePressed(Vector2F mouse_position) { return false; } // true = "this click is mine"
    virtual void onMouseDragged(Vector2F mouse_position) {}
    virtual void onMouseReleased(Vector2F mouse_position, bool released_inside) {}
    virtual void onMouseEntered() {}
    virtual void onMouseExited() {}

  protected:
    SDL_FRect m_rect = {0.0f, 0.0f, 0.0f, 0.0f};
    UIPlacement m_placement;

  private:
    UIElement* m_parent = nullptr;                      // the parent owns us
    std::vector<std::unique_ptr<UIElement>> m_children; // we own these
    std::optional<UILayer> m_layer;                     // empty = use the parent's layer
    int m_z_index = 0;
    bool m_is_visible = true;
};
} // namespace DynamoEngine