// File: Transfer/src/Systems/UISystem.cpp

#include "UISystem.h"

// Constructor
UISystem::UISystem() : allGameUIElements(), allPauseUIElements()
{
    // Initialize any UI system state here
    FPSCounter* fps_counter = new FPSCounter();
    allGameUIElements.push_back(fps_counter);
    MassSlider* mass_slider = new MassSlider();
    allGameUIElements.push_back(mass_slider);
    RadiusSlider* radius_slider = new RadiusSlider();
    allGameUIElements.push_back(radius_slider);
    SimulationSpeedSlider* simulation_speed_slider = new SimulationSpeedSlider();
    allGameUIElements.push_back(simulation_speed_slider);

    // fully implemented, but needs to go to a start menu scene or something
    // PlayGameButton* play_button = new PlayGameButton();
    // allGameUIElements.push_back(play_button);
}
// Destructor
UISystem::~UISystem()
{
    // Clean up any allocated resources here
}

void UISystem::UpdateUIElements(GameState& gameState, UIState& UIState)
{
    if (UIState.getGameScene())
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
    if (activeElement != UIElementType::NONE)
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
        if (inputsReceived.leftMouseButtonJustReleased) // You may need to pass this flag into InputState
        {
            UIElementType releasedOverElement = findElementWeAreIn(inputsReceived);
            if (releasedOverElement != UIElementType::NONE)
            {
                consumed = true;

                if (isButton(activeElement))
                {
                    routeButtonClick(activeElement, inputsReceived);
                }
                activeElement = UIElementType::NONE; // Reset for next interaction
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
    for (auto& element : allGameUIElements)
    {
        delete element;
    }
    for (auto& element : allStartGameUIElements)
    {
        delete element;
    }
    allGameUIElements.clear();
    allStartGameUIElements.clear();
}

void UISystem::routeSliderInput(UIElementType sliderTypeToUpdate, InputState& inputState)
{
    UIElement* elementToUpdate = allGameUIElements[sliderTypeToUpdate];

    // this could be further abstracted into passing the full input state so
    // that the UIelement decides what it updates, but I want to make the ui
    // element as stupid as possible

    // Also prob refactor as switch. there is duplicate code here so TODO: Cleanup
    if (sliderTypeToUpdate == UIElementType::MASS_SLIDER_INDEX)
    {
        double mass_val_to_be_calculated_and_injected = 0.0; // initialize
        elementToUpdate->slideMe(inputState.mouseCurrPosition, mass_val_to_be_calculated_and_injected);
        inputState.selectedMass = mass_val_to_be_calculated_and_injected;
    }
    else if (sliderTypeToUpdate == UIElementType::RADIUS_SLIDER_INDEX)
    {
        double radius_val_to_be_calculated_and_injected = 0.0; // initialize
        elementToUpdate->slideMe(inputState.mouseCurrPosition, radius_val_to_be_calculated_and_injected);
        inputState.selectedRadius = radius_val_to_be_calculated_and_injected;
    }
    else if (sliderTypeToUpdate == UIElementType::SIMULATION_SPEED_SLIDER_INDEX)
    {
        double simulation_speed_val_to_be_calculated_and_injected = 1.0; // initialize to 1
        elementToUpdate->slideMe(inputState.mouseCurrPosition, simulation_speed_val_to_be_calculated_and_injected);
        inputState.selectedSimSpeedScale = simulation_speed_val_to_be_calculated_and_injected;
    }
    else
    {
        return;
    }
}

void UISystem::routeButtonClick(UIElementType buttonToUpdate, InputState& inputState)
{
    UIElement* elementToUpdate = allGameUIElements[buttonToUpdate];
    if (buttonToUpdate == UIElementType::PLAY_GAME_BUTTON_INDEX)
    {
        elementToUpdate->clickMe(inputState.mouseCurrPosition);
    }
    // Will add other ui elements from in game here? Or maybe should just route to scene elements?
}
UIElementType UISystem::findElementWeAreIn(InputState& inputsReceived)
{
    // std::cout<<"findElementWeAreInCalled"<<std::endl;
    UIElementType contacted_element = UIElementType::NONE; // default result
    Vector2D curr_pos = inputsReceived.mouseCurrPosition;
    for (auto& element : allGameUIElements)
    {
        contacted_element = element->checkAndReturnIfHit(curr_pos);
        if (contacted_element != UIElementType::NONE)
        {
            break; // stop at the first hit
        }
    }
    return contacted_element; // single return
}

bool UISystem::isSlider(UIElementType typeToCheck)
{
    bool conclusion = false;
    if (typeToCheck == UIElementType::MASS_SLIDER_INDEX || typeToCheck == UIElementType::RADIUS_SLIDER_INDEX || typeToCheck == UIElementType::SIMULATION_SPEED_SLIDER_INDEX)
        conclusion = true;
    return conclusion;
}
bool UISystem::isButton(UIElementType typeToCheck)
{
    bool conclusion = false;
    if (typeToCheck == UIElementType::PLAY_GAME_BUTTON_INDEX)
        conclusion = true;
    return conclusion;
}