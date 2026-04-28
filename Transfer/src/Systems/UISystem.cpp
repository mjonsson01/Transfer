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
    PlayGameButton* play_button = new PlayGameButton();
    allGameUIElements.push_back(play_button);
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

    if (inputsReceived.leftMouseButtonJustPressed && inputsReceived.isClickingLeftMouseButton)
    {
        activeElement = findElementWeAreIn(inputsReceived);
    }
    // If you are clicking and its a preview macro, disable the previewing macro and consume?

    // If releasing on an element, still count as consumed
    if (!inputsReceived.isClickingLeftMouseButton && inputsReceived.leftMouseButtonJustReleased)
    {
        activeElement = findElementWeAreIn(inputsReceived);
        if (activeElement != UIElementType::NONE)
        {
            consumed = true;
        }
        activeElement = UIElementType::NONE; // release any hold on the element
    }
    if (activeElement != UIElementType::NONE && activeElement != FPS_COUNTER_INDEX)
    {
        updateSpecificElementAndPropagateUpwards(activeElement, inputsReceived);
        consumed = true;
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
    allGameUIElements.clear();
}

void UISystem::updateSpecificElementAndPropagateUpwards(UIElementType elementTypeToUpdate, InputState& inputState)
{
    UIElement* elementToUpdate = allGameUIElements[elementTypeToUpdate];

    // this could be further abstracted into passing the full input state so
    // that the UIelement decides what it updates, but I want to make the ui
    // element as stupid as possible
    if (elementTypeToUpdate == UIElementType::MASS_SLIDER_INDEX)
    {
        double massValToBeCalculatedAndInjected = 0.0; // initialize
        elementToUpdate->updateMe(inputState.mouseCurrPosition, massValToBeCalculatedAndInjected);
        inputState.selectedMass = massValToBeCalculatedAndInjected;
    }
    if (elementTypeToUpdate == UIElementType::RADIUS_SLIDER_INDEX)
    {
        double radiusValToBeCalculatedAndInjected = 0.0; // initialize
        elementToUpdate->updateMe(inputState.mouseCurrPosition, radiusValToBeCalculatedAndInjected);
        inputState.selectedRadius = radiusValToBeCalculatedAndInjected;
    }
    if (elementTypeToUpdate == UIElementType::PLAY_GAME_BUTTON_INDEX)
    {   
        std::cout<<"hit play game button"<<std::endl;
        double placeholder = 0.0;
        elementToUpdate->updateMe(inputState.mouseCurrPosition, placeholder);
    }
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
