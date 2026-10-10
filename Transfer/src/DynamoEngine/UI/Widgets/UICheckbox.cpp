// File: Transfer/src/DynamoEngine/UI/Widgets/UICheckbox.cpp

#include "DynamoEngine/UI/Widgets/UICheckbox.hpp"

namespace
{
constexpr SDL_Color CHECKBOX_COLOR = {128, 128, 128, 255};         // the same gray as a button
constexpr SDL_Color CHECKBOX_HOVERED_COLOR = {108, 108, 108, 255}; // slightly darker
constexpr SDL_Color CHECKBOX_HELD_COLOR = {88, 88, 88, 255};       // darker still while held down
constexpr SDL_Color BOX_COLOR = {50, 50, 50, 255};                 // the empty square
constexpr SDL_Color CHECK_COLOR = {255, 255, 255, 255};            // the white fill while checked
constexpr SDL_Color LABEL_COLOR = {255, 255, 255, 255};            // white

constexpr float BOX_SIZE_FRACTION = 0.6f;    // the square's side, as a share of the row's height
constexpr float CHECK_INSET_FRACTION = 0.2f; // the gap around the white fill, as a share of the square's side
} // namespace

namespace DynamoEngine
{
UICheckbox::UICheckbox(std::string label, bool starts_checked) : m_label(std::move(label)), m_is_checked(starts_checked)
{
}

void UICheckbox::update(float delta_seconds)
{
    if (m_checked_source)
    {
        m_is_checked = m_checked_source();
    }
}

SDL_FRect UICheckbox::boxRect() const
{
    // A square, vertically centred, with the same gap to the left edge as above and below it
    const float side = m_rect.h * BOX_SIZE_FRACTION;
    const float gap = (m_rect.h - side) / 2.0f;
    return SDL_FRect{m_rect.x + gap, m_rect.y + gap, side, side};
}

void UICheckbox::draw(UIGeometryBuilder& builder) const
{
    SDL_Color background = CHECKBOX_COLOR;
    if (m_is_held_down)
    {
        background = CHECKBOX_HELD_COLOR;
    }
    else if (m_is_hovered)
    {
        background = CHECKBOX_HOVERED_COLOR;
    }
    builder.addRect(m_rect, background);

    const SDL_FRect box = boxRect();
    builder.addRect(box, BOX_COLOR);
    if (m_is_checked)
    {
        const float inset = box.w * CHECK_INSET_FRACTION;
        builder.addRect(SDL_FRect{box.x + inset, box.y + inset, box.w - (2.0f * inset), box.h - (2.0f * inset)},
                        CHECK_COLOR);
    }

    // The label starts one gap to the right of the box, vertically centred in the row
    const float gap = box.x - m_rect.x;
    const Vector2F label_top_left = {box.x + box.w + gap, m_rect.y + ((m_rect.h - builder.fontHeight()) / 2.0f)};
    builder.addText(m_label, label_top_left, LABEL_COLOR);
}

bool UICheckbox::onMousePressed(Vector2F mouse_position)
{
    m_is_held_down = true;
    return true; // checkboxes always take the click, like buttons
}

void UICheckbox::onMouseReleased(Vector2F mouse_position, bool released_inside)
{
    m_is_held_down = false;
    if (!released_inside)
    {
        return; // the player backed out of the click
    }

    m_is_checked = !m_is_checked;
    requestSound(UISound::StandardClick);
    if (m_on_toggled) // an empty action is "false": nobody is listening
    {
        m_on_toggled(m_is_checked);
    }
}
} // namespace DynamoEngine