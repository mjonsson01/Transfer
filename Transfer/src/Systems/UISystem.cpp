// File: Transfer/src/Systems/UISystem.cpp

#include "UISystem.h"

// Constructor
UISystem::UISystem() : allUIElements()
{
    // Initialize any UI system state here
    FPSCounter* fps_counter = new FPSCounter();
    allUIElements.push_back(fps_counter);
    MassSlider* mass_slider = new MassSlider();
    allUIElements.push_back(mass_slider);
    RadiusSlider* radius_slider = new RadiusSlider();
    allUIElements.push_back(radius_slider);
}
// Destructor
UISystem::~UISystem()
{
    // Clean up any allocated resources here
}

void UISystem::UpdateUIElements(GameState& gameState, UIState& UIState)
{
    bool consumed = false;
    InputState& inputsReceived = UIState.getMutableInputState();

    if (inputsReceived.leftMouseButtonJustPressed && inputsReceived.isClickingLeftMouseButton)
    {
        activeElement = isPositionInUIElementHotZone(inputsReceived);
    }
    // If you are clicking and its a preview macro, disable the previewing macro and consume?

    // If releasing on an element, still count as consumed
    if (!inputsReceived.isClickingLeftMouseButton || inputsReceived.leftMouseButtonJustReleased)
    {
        activeElement = isPositionInUIElementHotZone(inputsReceived);
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

void UISystem::CleanUp()
{
    for (auto& element : allUIElements)
    {
        delete element;
    }
    allUIElements.clear();
}

void UISystem::updateSpecificElementAndPropagateUpwards(UIElementType elementTypeToUpdate, InputState& inputState)
{
    UIElement* elementToUpdate = allUIElements[elementTypeToUpdate];

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
    // if (elementTypeToUpdate == UIElementType::RADIUS_SLIDER_INDEX)
    // {
    // }
}

UIElementType UISystem::isPositionInUIElementHotZone(InputState& inputsReceived)
{
    UIElementType hotzoneType = UIElementType::NONE; // default result
    Vector2D curr_pos = inputsReceived.mouseCurrPosition;
    for (auto& element : allUIElements)
    {
        hotzoneType = element->isInDeadZone(curr_pos);
        if (hotzoneType != UIElementType::NONE)
        {
            break; // stop at the first hit
        }
    }
    return hotzoneType; // single return
}

// void UISystem::ProcessUIFrame(GameState& gameState, UIState& UIState)
// {

// }

// void UISystem::InitializeUIElements(UIState& UIState)
// {
//     // Create and add the FPS counter
//     FPSCounter* fpsCounter = new FPSCounter();
//     UIState.addUIElement(fpsCounter);
//     // Create and add the Mass Slider
//     MassSlider* massSlider = new MassSlider();
//     UIState.addUIElement(massSlider);
//     // Add other UI elements as needed
// }

// void UISystem::DeleteUIElements(UIState& UIState)
// {
//     // Get all UI elements and delete them
//     std::vector<UIElement*> ui_elements = UIState.getUIElements();
//     for (auto& element : ui_elements) {
//         delete element;
//     }
//     UIState.clearUIElements();
// }

// void UISystem::RenderUIElements(SDL_Renderer* renderer, UIState& UIState,
// TTF_Font* UIFont)
// {
//     std::vector<UIElement*> ui_elements = UIState.getUIElements();
//     if (!UIState.getUIElementsVisible()) return; // Skip rendering if UI
//     elements are hidden

//     for (auto& element : ui_elements) {
//         element->renderElement(renderer, UIState, UIFont); // Pass UIFont if
//         needed
//     }
// }
