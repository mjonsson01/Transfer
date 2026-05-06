// File: Transfer/src/Systems/UISystem.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Core/GameState.h"
#include "Core/InputState.h"
#include "Core/UIState.h"
#include "Entities/UIElements/Buttons/PlayGameButton.h"
#include "Entities/UIElements/Overlay/FPSCounter.h"
#include "Entities/UIElements/Sliders/MassSlider.h"
#include "Entities/UIElements/Sliders/RadiusSlider.h"
#include "Entities/UIElements/Sliders/SimulationSpeedSlider.h"
#include "Entities/UIElements/UIElement.h"

// Standard Library Imports
#include <vector>

// Owns Logic of UI Components and stores the UIElements while dispatching
// information to the GameState and UIState
class UISystem
{
  public:
    UISystem();
    ~UISystem();
    void CleanUp();
    void UpdateUIElements(GameState& gameState, UIState& UIState);
    void updateGameUIElements(GameState& gameState, UIState& UIState);
    std::vector<UIElement*>& getGameUIElementsMutable() { return allGameUIElements; }
    const std::vector<UIElement*>& getGameUIElements() const { return allGameUIElements; }
    void updatePauseUIElements(GameState& gameState, UIState& UIState);
    std::vector<UIElement*>& getPauseUIElementsMutable() { return allGameUIElements; }
    const std::vector<UIElement*>& getPauseUIElements() const { return allGameUIElements; }
    UIElementType findElementWeAreIn(InputState& inputsReceived);
    void routeSliderInput(UIElementType elementToUpdate, InputState& inputState);
    void routeButtonClick(UIElementType elementToUpdate, InputState& inputState);

  private:
    bool isSlider(UIElementType typeToCheck);
    bool isButton(UIElementType typeToCheck);
    std::vector<UIElement*> allGameUIElements;
    std::vector<UIElement*> allPauseUIElements;
    std::vector<UIElement*> allStartGameUIElements;
    std::vector<UIElement*> allLevelEditorUIElements;
    UIElementType activeElement = UIElementType::NONE;
};
