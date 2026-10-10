// File: Transfer/src/DynamoEngine/UI/Widgets/UIButton.cpp

#include "DynamoEngine/UI/Widgets/UIButton.hpp"

namespace
{
constexpr SDL_Color BUTTON_COLOR = {128, 128, 128, 255};         // gray
constexpr SDL_Color BUTTON_HOVERED_COLOR = {108, 108, 108, 255}; // slightly darker
constexpr SDL_Color BUTTON_HELD_COLOR = {88, 88, 88, 255};       // darker still while held down
constexpr SDL_Color BUTTON_TEXT_COLOR = {255, 255, 255, 255};    // white
} // namespace

namespace DynamoEngine
{
void UIButton::draw(UIGeometryBuilder& builder) const
{
    SDL_Color background = BUTTON_COLOR;
    if (m_is_held_down)
    {
        background = BUTTON_HELD_COLOR;
    }
    else if (m_is_hovered)
    {
        background = BUTTON_HOVERED_COLOR;
    }

    builder.addRect(m_rect, background);
    builder.addTextCentered(m_text, m_rect, BUTTON_TEXT_COLOR);
}

bool UIButton::onMousePressed(Vector2F mouse_position)
{
    m_is_held_down = true;
    return true; // buttons always take the click
}

void UIButton::onMouseReleased(Vector2F mouse_position, bool released_inside)
{
    m_is_held_down = false;

    if (released_inside)
    {
        requestSound(UISound::SpecialClick);
        if (m_on_click) // an empty UIAction is "false": nothing to do
        {
            m_on_click();
        }
    }
}
} // namespace DynamoEngine
