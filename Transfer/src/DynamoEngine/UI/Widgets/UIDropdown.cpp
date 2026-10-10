// File: Transfer/src/DynamoEngine/UI/Widgets/UIDropdown.cpp

#include "DynamoEngine/UI/Widgets/UIDropdown.hpp"

// Standard Library Imports
#include <cassert>
#include <utility>

namespace
{
using DynamoEngine::Vector2F;

// Same look as UIButton
constexpr SDL_Color BUTTON_COLOR = {128, 128, 128, 255};         // gray
constexpr SDL_Color BUTTON_HOVERED_COLOR = {108, 108, 108, 255}; // slightly darker
constexpr SDL_Color BUTTON_HELD_COLOR = {88, 88, 88, 255};       // darker still while held down
constexpr SDL_Color BUTTON_TEXT_COLOR = {255, 255, 255, 255};    // white

// Same edge rule as UIElement::containsPoint: left and top edges are inside, right and bottom are not
bool isInside(const SDL_FRect& rect, Vector2F point)
{
    const bool inside_x = point.x_val >= rect.x && point.x_val < rect.x + rect.w;
    const bool inside_y = point.y_val >= rect.y && point.y_val < rect.y + rect.h;
    return inside_x && inside_y;
}
} // namespace

namespace DynamoEngine
{
UIDropdown::UIDropdown(std::string title, std::vector<std::string> option_labels)
    : m_title(std::move(title)), m_option_labels(std::move(option_labels))
{
    assert(!m_option_labels.empty() && "UIDropdown needs at least one option");
    setLayer(UILayer::Overlay);
}

void UIDropdown::selectOption(int option_index)
{
    const bool is_valid_index = option_index >= 0 && option_index < static_cast<int>(m_option_labels.size());
    assert(is_valid_index && "selectOption: there is no option with this index");
    if (is_valid_index)
    {
        m_selected_option = option_index;
    }
}

void UIDropdown::update(float delta_seconds)
{
    if (m_selection_source)
    {
        selectOption(m_selection_source());
    }
}

// Geometry, button is m_rect, option rows stack directly below with the same size

bool UIDropdown::isOnButton(Vector2F point) const
{
    return UIElement::containsPoint(point); // normal inside rect test (not the overridden class version)
}

SDL_FRect UIDropdown::optionRect(int option_index) const
{
    const float row_y = m_rect.y + m_rect.h * static_cast<float>(option_index + 1); // row 0 sits right under the button
    return {m_rect.x, row_y, m_rect.w, m_rect.h};
}

int UIDropdown::optionAt(Vector2F point) const
{
    if (!m_is_open)
    {
        return NO_OPTION; // closed, no row options for selection
    }
    for (int option_index = 0; option_index < static_cast<int>(m_option_labels.size()); ++option_index)
    {
        if (isInside(optionRect(option_index), point))
        {
            return option_index;
        }
    }
    return NO_OPTION;
}

bool UIDropdown::containsPoint(Vector2F point) const
{
    if (m_is_open)
    {
        return true; // while its open, every press is routed here
    }
    return isOnButton(point);
}

bool UIDropdown::onMousePressed(Vector2F mouse_position)
{
    if (!m_is_open)
    {
        // Closed: containsPoint only lets this press through because its on the button, so open the list
        m_is_open = true;
        m_is_button_held = true;
        requestSound(UISound::StandardClick);
        return true;
    }

    // Open: this press closes the list, whatever it landed on
    const int pressed_option = optionAt(mouse_position); // before closing
    m_is_open = false;
    m_hovered_option = NO_OPTION;

    if (isOnButton(mouse_position))
    {
        m_is_button_held = true;
        requestSound(UISound::StandardClick);
    }
    else if (pressed_option != NO_OPTION)
    {
        requestSound(UISound::StandardClick);
        const bool is_new_choice = (pressed_option != m_selected_option);
        m_selected_option = pressed_option;
        if (is_new_choice && m_on_option_chosen)
        {
            m_on_option_chosen(pressed_option);
        }
    }
    // Game never acts on the press that closed the menu (outside the press)
    return true;
}

void UIDropdown::onMouseReleased(Vector2F mouse_position, bool released_inside) { m_is_button_held = false; }

void UIDropdown::onMouseHover(Vector2F mouse_position)
{
    // While open we contain every point so check the actual buttons and rows
    m_is_button_hovered = isOnButton(mouse_position);
    m_hovered_option = optionAt(mouse_position);
}

void UIDropdown::onMouseExited()
{
    m_is_button_hovered = false;
    m_hovered_option = NO_OPTION;
}

void UIDropdown::draw(UIGeometryBuilder& builder) const
{
    SDL_Color button_color = BUTTON_COLOR;
    if (m_is_button_held)
    {
        button_color = BUTTON_HELD_COLOR;
    }
    else if (m_is_button_hovered)
    {
        button_color = BUTTON_HOVERED_COLOR;
    }
    builder.addRect(m_rect, button_color);
    builder.addTextCentered(m_title + ": " + m_option_labels[m_selected_option], m_rect, BUTTON_TEXT_COLOR);

    if (!m_is_open)
    {
        return;
    }

    // Draw open list
    for (int option_index = 0; option_index < static_cast<int>(m_option_labels.size()); ++option_index)
    {
        const SDL_FRect row = optionRect(option_index);
        const SDL_Color row_color = (option_index == m_hovered_option) ? BUTTON_HOVERED_COLOR : BUTTON_COLOR;
        builder.addRect(row, row_color);
        builder.addTextCentered(m_option_labels[option_index], row, BUTTON_TEXT_COLOR);
    }
}
} // namespace DynamoEngine