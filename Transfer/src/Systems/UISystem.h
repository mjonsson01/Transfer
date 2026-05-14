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
#include "Scenes/GameScene/GameScene.h"
#include "Scenes/PauseScene/PauseScene.h"
#include "Scenes/Scene.h"
// Standard Library Imports
#include <iostream>
#include <unordered_map>

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
    std::unordered_map<UIElementIdentifier, UIElement*>& getAllUIElements() { return allUIElements; }
    Scene* getScene(SceneIdentifier sceneID) { return allScenes[sceneID]; }
    void updatePauseUIElements(GameState& gameState, UIState& UIState);

    UIElementIdentifier findElementWeAreIn(InputState& inputsReceived);
    void routeSliderInput(UIElementIdentifier elementToUpdate, InputState& inputState);
    void routeButtonClick(UIElementIdentifier elementToUpdate, InputState& inputState);
    void populateScenes();

  private:
    bool isSlider(UIElementIdentifier typeToCheck);
    bool isButton(UIElementIdentifier typeToCheck);
    std::unordered_map<SceneIdentifier, Scene*> allScenes;
    std::unordered_map<UIElementIdentifier, UIElement*> allUIElements;
};
