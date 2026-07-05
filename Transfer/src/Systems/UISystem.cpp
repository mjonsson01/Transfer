// File: Transfer/src/Systems/UISystem.cpp

#include "UISystem.hpp"

#include <iostream>

// Constructor
UISystem::UISystem() : allScenes()
{
    allScenes.insert({SceneIdentifier::GAME_SCENE, nullptr});
    allScenes.insert({SceneIdentifier::PAUSE_SCENE, nullptr});
    allScenes.insert({SceneIdentifier::START_MENU_SCENE, nullptr});
    allScenes.insert({SceneIdentifier::TEST_VISUAL_SCENE, nullptr});

    populateScenes();
}
// Destructor
UISystem::~UISystem()
{
    // Clean up any allocated resources here
}

void UISystem::UpdateUIElements(GameState& gameState, UIState& UIState)
{
    updateUISystemCurrentSceneID(UIState); // get most up to date current scene

    // Clunky fix but whatever
    std::unordered_map<UIElementIdentifier, UIElement*> allUIElements = allScenes[currentSceneID]->getSceneElements();
    for (auto& [UI_element_ID, UI_element_ptr] : allUIElements)
    {
        if (UI_element_ID == UIElementIdentifier::FPS_COUNTER_INDEX)
        {
            UI_element_ptr->updateMe(UIState);
            break;
        }
    }

    if (currentSceneID == SceneIdentifier::GAME_SCENE)
    {
        updateGameUIElements(gameState, UIState);
    }
    else
    {
        updateMenuUIElements(gameState, UIState);
    }
}

void UISystem::CleanUp()
{
    for (auto& [scene_ID, scene_ptr] : allScenes)
    {
        if (scene_ptr)
        {
            scene_ptr->CleanUpSceneElements();
            delete scene_ptr;
            scene_ptr = nullptr;
        }
    }
    allScenes.clear();
}

void UISystem::populateScenes()
{
    GameScene* game_scene = new GameScene();
    PauseScene* pause_scene = new PauseScene();
    StartMenuScene* start_menu_scene = new StartMenuScene();
    TestVisualScene* test_visual_scene = new TestVisualScene();
    for (auto& [scene_ID, scene_ptr] : allScenes)
    {
        switch (scene_ID)
        {
        case SceneIdentifier::GAME_SCENE:
            scene_ptr = game_scene;
            scene_ptr->populateMe();
            break;
        case SceneIdentifier::PAUSE_SCENE:
            scene_ptr = pause_scene;
            scene_ptr->populateMe();
            break;
        case SceneIdentifier::START_MENU_SCENE:
            scene_ptr = start_menu_scene;
            scene_ptr->populateMe();
            break;
        case SceneIdentifier::TEST_VISUAL_SCENE:
            scene_ptr = test_visual_scene;
            scene_ptr->populateMe();
            break;
        default:
            break;
        }
    }
}

void UISystem::updateGameUIElements(GameState& gameState, UIState& UIState)
{
    bool consumed = false;
    InputState& inputs_received = UIState.getMutableInputState();

    if (inputs_received.leftMouseButtonJustPressed)
    {
        activeElementID = findElementWeAreIn(inputs_received);
    }
    if (activeElementID != UIElementIdentifier::NONE)
    {
        consumed = true;
        if (inputs_received.isClickingLeftMouseButton)
        {
            if (isSlider(activeElementID))
            {
                routeSliderInput(activeElementID, UIState);
            }
        }
        if (inputs_received.leftMouseButtonJustReleased)
        {
            UIElementIdentifier releasedOver = findElementWeAreIn(inputs_received);

            if (releasedOver == activeElementID)
            {
                if (isButton(activeElementID))
                {
                    routeButtonClick(activeElementID, UIState);
                }
            }
            activeElementID = UIElementIdentifier::NONE;
        }
    }
    inputs_received.UIInputConsumed = consumed;
    inputs_received.isPreviewingMacro = inputs_received.isClickingLeftMouseButton && !consumed;
    inputs_received.isPreviewingWithInitialVelocity =
        inputs_received.isPreviewingMacro && inputs_received.isPressingShift;
}

void UISystem::updateMenuUIElements(GameState& gameState, UIState& UIState)
{
    bool consumed = false;
    InputState& inputs_received = UIState.getMutableInputState();

    if (inputs_received.leftMouseButtonJustPressed)
    {
        activeElementID = findElementWeAreIn(inputs_received);
    }
    if (activeElementID != UIElementIdentifier::NONE)
    {
        consumed = true;
        if (inputs_received.isClickingLeftMouseButton)
        {
            if (isSlider(activeElementID))
            {
                routeSliderInput(activeElementID, UIState);
            }
        }
        if (inputs_received.leftMouseButtonJustReleased)
        {
            UIElementIdentifier releasedOver = findElementWeAreIn(inputs_received);

            if (releasedOver == activeElementID)
            {
                if (isButton(activeElementID))
                {
                    routeButtonClick(activeElementID, UIState);
                }
            }
            activeElementID = UIElementIdentifier::NONE;
        }
    }
    inputs_received.UIInputConsumed = consumed;
    inputs_received.isPreviewingMacro = false;
    inputs_received.isPreviewingWithInitialVelocity = false;
}

void UISystem::routeSliderInput(UIElementIdentifier sliderTypeToUpdate, UIState& UIState)
{
    std::unordered_map<UIElementIdentifier, UIElement*> allUIElements = allScenes[currentSceneID]->getSceneElements();
    UIElement* elementToUpdate = allUIElements[sliderTypeToUpdate];
    auto& input_state = UIState.getMutableInputState();
    // this could be further abstracted into passing the full input state so
    // that the UIelement decides what it updates, but I want to make the ui
    // element as stupid as possible

    // Also prob refactor as switch. there is duplicate code here so TODO: Cleanup
    if (sliderTypeToUpdate == UIElementIdentifier::MASS_SLIDER_INDEX)
    {
        double mass_val_to_be_calculated_and_injected = 0.0; // initialize
        elementToUpdate->slideMe(input_state.mouseCurrPosition, mass_val_to_be_calculated_and_injected, UIState);
        input_state.selectedMass = mass_val_to_be_calculated_and_injected;
    }

    else if (sliderTypeToUpdate == UIElementIdentifier::RADIUS_SLIDER_INDEX)
    {
        double radius_val_to_be_calculated_and_injected = 0.0; // initialize
        elementToUpdate->slideMe(input_state.mouseCurrPosition, radius_val_to_be_calculated_and_injected, UIState);
        input_state.selectedRadius = radius_val_to_be_calculated_and_injected;
    }
    else if (sliderTypeToUpdate == UIElementIdentifier::SIMULATION_SPEED_SLIDER_INDEX)
    {
        double simulation_speed_val_to_be_calculated_and_injected = 1.0; // initialize to 1
        elementToUpdate->slideMe(input_state.mouseCurrPosition, simulation_speed_val_to_be_calculated_and_injected,
                                 UIState);
        input_state.selectedSimSpeedScale = simulation_speed_val_to_be_calculated_and_injected;
    }
    else
    {
        // no known identifier, just exit.
        return;
    }
}

void UISystem::routeButtonClick(UIElementIdentifier buttonToUpdate, UIState& UIState)
{
    InputState& input_state = UIState.getMutableInputState();
    std::unordered_map<UIElementIdentifier, UIElement*> allUIElements = allScenes[currentSceneID]->getSceneElements();
    UIElement* elementToUpdate = allUIElements[buttonToUpdate];
    if (buttonToUpdate == UIElementIdentifier::PLAY_GAME_BUTTON_INDEX)
    {
        elementToUpdate->clickMe(input_state.mouseCurrPosition, UIState);
    }
    else if (buttonToUpdate == UIElementIdentifier::RESUME_BUTTON_INDEX)
    {
        elementToUpdate->clickMe(input_state.mouseCurrPosition, UIState);
    }
    else
    {
        return;
    }
    // Will add other ui elements from in game here? Or maybe should just route to scene elements?
}
UIElementIdentifier UISystem::findElementWeAreIn(InputState& inputsReceived)
{

    UIElementIdentifier contacted_element = UIElementIdentifier::NONE; // default result
    Vector2D curr_pos = inputsReceived.mouseCurrPosition;
    std::unordered_map<UIElementIdentifier, UIElement*> allUIElements = allScenes[currentSceneID]->getSceneElements();
    for (auto& [UI_element_ID, UI_element_ptr] : allUIElements)
    {
        contacted_element = UI_element_ptr->checkAndReturnIfHit(curr_pos);
        if (contacted_element != UIElementIdentifier::NONE)
        {
            break; // stop at the first hit
        }
    }
    return contacted_element; // single return
}

bool UISystem::isSlider(UIElementIdentifier typeToCheck)
{
    bool conclusion = false;
    if (typeToCheck == UIElementIdentifier::MASS_SLIDER_INDEX ||
        typeToCheck == UIElementIdentifier::RADIUS_SLIDER_INDEX ||
        typeToCheck == UIElementIdentifier::SIMULATION_SPEED_SLIDER_INDEX)
        conclusion = true;
    return conclusion;
}
bool UISystem::isButton(UIElementIdentifier typeToCheck)
{
    bool conclusion = false;
    if (typeToCheck == UIElementIdentifier::PLAY_GAME_BUTTON_INDEX ||
        typeToCheck == UIElementIdentifier::RESUME_BUTTON_INDEX)
        conclusion = true;
    return conclusion;
}