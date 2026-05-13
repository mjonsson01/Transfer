// File: Transfer/src/Systems/UISystem.cpp

#include "UISystem.h"

// Constructor
UISystem::UISystem() : allScenes()
{
    // FPSCounter* fps_counter = new FPSCounter();
    // allUIElements.insert({fps_counter->getUIElementIdentifier(), fps_counter});
    // MassSlider* mass_slider = new MassSlider();
    // allUIElements.insert({mass_slider->getUIElementIdentifier(), mass_slider});
    // RadiusSlider* radius_slider = new RadiusSlider();
    // allUIElements.insert({radius_slider->getUIElementIdentifier(), radius_slider});
    // SimulationSpeedSlider* simulation_speed_slider = new SimulationSpeedSlider();
    // allUIElements.insert({simulation_speed_slider->getUIElementIdentifier(), simulation_speed_slider});
    // PlayGameButton* play_button = new PlayGameButton();
    // allGameUIElements.push_back(play_button);
    allScenes.insert({START_MENU_SCENE, nullptr});
    allScenes.insert({GAME_SCENE, nullptr});
    allScenes.insert({PAUSE_SCENE, nullptr});
    allScenes.insert({TEST_VISUAL_SCENE, nullptr});
}
// Destructor
UISystem::~UISystem()
{
    // Clean up any allocated resources here
}

void UISystem::UpdateUIElements(GameState& gameState, UIState& UIState)
{
    if (UIState.getGameScene()) // TODO convert to switch of game scene?
    {
        if (UIState.getPauseMenuActive())
        {
            updatePauseUIElements(gameState, UIState);
        }
        else
        {
            updateGameUIElements(gameState, UIState);
        }
    }
    else
    {
        return;
    }
}
void UISystem::populateScenes() {}
void UISystem::updateGameUIElements(GameState& gameState, UIState& UIState)
{
    bool consumed = false;
    InputState& inputsReceived = UIState.getMutableInputState();

    // Grab elements on just pressed
    if (inputsReceived.leftMouseButtonJustPressed)
    {
        activeElement = findElementWeAreIn(inputsReceived);
    }
    // immediately consume if clicking on an active element
    if (activeElement != UIElementIdentifier::NONE)
    {
        consumed = true;
    }
    // continuously update elements that are draggable
    if (inputsReceived.isClickingLeftMouseButton && inputsReceived.isDragging)
    {
        if (isSlider(activeElement) && activeElement != FPS_COUNTER_INDEX)
        {
            routeSliderInput(activeElement, inputsReceived);
            consumed = true;
        }
    }
    else
    {
        // If we just released the mouse
        if (inputsReceived.leftMouseButtonJustReleased)
        {
            UIElementIdentifier releasedOverElement = findElementWeAreIn(inputsReceived);
            if (releasedOverElement != UIElementIdentifier::NONE)
            {
                consumed = true;

                if (isButton(activeElement))
                {
                    routeButtonClick(activeElement, inputsReceived);
                }
                activeElement = UIElementIdentifier::NONE; // Reset for next interaction
            }
        }
    }
    inputsReceived.UIInputConsumed = consumed;
    inputsReceived.isPreviewingMacro = inputsReceived.isClickingLeftMouseButton && !consumed;
    inputsReceived.isPreviewingWithInitialVelocity = inputsReceived.isPreviewingMacro && inputsReceived.isPressingShift;
}

void UISystem::updatePauseUIElements(GameState& gameState, UIState& UIState) { return; }

void UISystem::CleanUp()
{
    for (auto& [id, scene_ptr] : allScenes)
    {
        if (scene_ptr)
        {
            scene_ptr->CleanUpSceneElements();
            delete scene_pair.second;
            scene_ptr = nullptr
        }
    }
    allScenes.clear();
}

void UISystem::routeSliderInput(UIElementIdentifier sliderTypeToUpdate, InputState& inputState)
{
    UIElement* elementToUpdate = allUIElements.find(sliderTypeToUpdate)->second;

    // this could be further abstracted into passing the full input state so
    // that the UIelement decides what it updates, but I want to make the ui
    // element as stupid as possible

    // Also prob refactor as switch. there is duplicate code here so TODO: Cleanup
    if (sliderTypeToUpdate == UIElementIdentifier::MASS_SLIDER_INDEX)
    {
        double mass_val_to_be_calculated_and_injected = 0.0; // initialize
        elementToUpdate->slideMe(inputState.mouseCurrPosition, mass_val_to_be_calculated_and_injected);
        inputState.selectedMass = mass_val_to_be_calculated_and_injected;
    }
    else if (sliderTypeToUpdate == UIElementIdentifier::RADIUS_SLIDER_INDEX)
    {
        double radius_val_to_be_calculated_and_injected = 0.0; // initialize
        elementToUpdate->slideMe(inputState.mouseCurrPosition, radius_val_to_be_calculated_and_injected);
        inputState.selectedRadius = radius_val_to_be_calculated_and_injected;
    }
    else if (sliderTypeToUpdate == UIElementIdentifier::SIMULATION_SPEED_SLIDER_INDEX)
    {
        double simulation_speed_val_to_be_calculated_and_injected = 1.0; // initialize to 1
        elementToUpdate->slideMe(inputState.mouseCurrPosition, simulation_speed_val_to_be_calculated_and_injected);
        inputState.selectedSimSpeedScale = simulation_speed_val_to_be_calculated_and_injected;
    }
    else
    {
        // no known identifier, just exit.
        return;
    }
}

void UISystem::routeButtonClick(UIElementIdentifier buttonToUpdate, InputState& inputState)
{
    UIElement* elementToUpdate = allUIElements.find(buttonToUpdate)->second;
    if (buttonToUpdate == UIElementIdentifier::PLAY_GAME_BUTTON_INDEX)
    {
        elementToUpdate->clickMe(inputState.mouseCurrPosition);
    }
    // Will add other ui elements from in game here? Or maybe should just route to scene elements?
}
UIElementIdentifier UISystem::findElementWeAreIn(InputState& inputsReceived)
{
    // std::cout<<"findElementWeAreInCalled"<<std::endl;
    UIElementIdentifier contacted_element = UIElementIdentifier::NONE; // default result
    Vector2D curr_pos = inputsReceived.mouseCurrPosition;
    for (auto& pair : allUIElements)
    {
        contacted_element = pair.second->checkAndReturnIfHit(curr_pos);
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
    if (typeToCheck == UIElementIdentifier::PLAY_GAME_BUTTON_INDEX)
        conclusion = true;
    return conclusion;
}