// File: Transfer/src/DynamoEngine/UI/Widgets/UIDropdown.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/UI/UIElement.hpp"

// Standard Library Imports
#include <functional>
#include <string>
#include <vector>

namespace DynamoEngine
{
// A collapsible set of buttons showing the current choice. Pressing it pops open the list of sub options as buttons.
// A click on a button should choose it and collapse the list. Clicking anywhere else should close the menu
// Draws in the Overlay layer

class UIDropdown : public UIElement
{
  public:
    using OptionChosenAction = std::function<void(int option_index)>;
    using SelectionSource = std::function<int()>;

    UIDropdown(std::string title, std::vector<std::string> option_labels); // option 0 should be default at construction

    // if the option differs from the current one, the option changes, otherwise nothing changes
    void setOnOptionChosen(OptionChosenAction action) { m_on_option_chosen = std::move(action); }

    // Where the selection comes from every frame
    void setSelectionSource(SelectionSource source) { m_selection_source = std::move(source); }

    void selectOption(int option_index);                     // selects without calling the chosen action
    int selectedOption() const { return m_selected_option; } // getter for the selected option
    bool isOpen() const { return m_is_open; }

    // UIElement overrides
    bool containsPoint(Vector2F point) const override; // while open: EVERY point, so a press outside can close it
    void update(float delta_seconds) override;
    void draw(UIGeometryBuilder& builder) const override;
    bool onMousePressed(Vector2F mouse_position) override;
    void onMouseReleased(Vector2F mouse_position, bool released_inside) override;
    void onMouseHover(Vector2F mouse_position) override;
    void onMouseExited() override;

  private:
    static constexpr int NO_OPTION = -1;

    bool isOnButton(Vector2F point) const;        // if on the top level button that pops open the dropdown
    SDL_FRect optionRect(int option_index) const; // option idx row stacked below button
    int optionAt(Vector2F point) const;           // Option row under point or NO_OPTION;

    std::string m_title;
    std::vector<std::string> m_option_labels;
    OptionChosenAction m_on_option_chosen; // empty = nothing is listening
    SelectionSource m_selection_source;    // empty = only selectOption() and player change it

    int m_selected_option = 0;
    int m_hovered_option = NO_OPTION;
    bool m_is_open = false;
    bool m_is_button_hovered = false;
    bool m_is_button_held = false;
};
} // namespace DynamoEngine