// File: Transfer/src/Systems/UISystem.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "Core/GameState.hpp"
#include "Core/InputState.hpp"
#include "Core/UIState.hpp"
#include "Entities/UIElements/UIElement.hpp"
#include "Scenes/GameScene/GameScene.hpp"
#include "Scenes/PauseScene/PauseScene.hpp"
#include "Scenes/Scene.hpp"
#include "Scenes/SceneIdentifierEnum.hpp"
#include "Scenes/StartMenuScene/StartMenuScene.hpp"
#include "Scenes/TestVisualScene/TestVisualScene.hpp"

// Standard Library Imports
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
    Scene* getScene(SceneIdentifier sceneID) { return allScenes[sceneID]; }
    void updateUISystemCurrentSceneID(UIState& UIState) { currentSceneID = UIState.getCurrentSceneID(); }
    void updateGameUIElements(GameState& gameState, UIState& UIState);
    void updateMenuUIElements(GameState& gameState, UIState& UIState);
    UIElementIdentifier findElementWeAreIn(InputState& inputsReceived);
    void routeSliderInput(UIElementIdentifier elementToUpdate, InputState& inputState);
    void routeButtonClick(UIElementIdentifier elementToUpdate, UIState& UIState);
    void populateScenes();

  private:
    SceneIdentifier currentSceneID = SceneIdentifier::NONE;
    bool isSlider(UIElementIdentifier typeToCheck);
    bool isButton(UIElementIdentifier typeToCheck);
    std::unordered_map<SceneIdentifier, Scene*> allScenes;
    UIElementIdentifier activeElementID = UIElementIdentifier::NONE;
};
