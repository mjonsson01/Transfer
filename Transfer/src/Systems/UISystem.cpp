// File: Transfer/src/Systems/UISystem.cpp

#include "UISystem.h"
#include "Entities/UIElement.h"
#include "Entities/UIElements/FPSCounter.h"
#include "Entities/UIElements/MassSlider.h"

// Constructor
UISystem::UISystem()
{
    // Initialize any UI system state here
}
// Destructor
UISystem::~UISystem()
{
    // Clean up any allocated resources here
}

void UISystem::ProcessUIFrame(GameState& state, UIState& UIState)
{
    
}

void UISystem::InitializeUIElements(UIState& UIState)
{
    // Create and add the FPS counter
    FPSCounter* fpsCounter = new FPSCounter();
    UIState.addUIElement(fpsCounter);
    // Create and add the Mass Slider
    MassSlider* massSlider = new MassSlider();
    UIState.addUIElement(massSlider);
    // Add other UI elements as needed
}

void UISystem::DeleteUIElements(UIState& UIState)
{
    // Get all UI elements and delete them
    std::vector<UIElement*> ui_elements = UIState.getUIElements();
    for (auto& element : ui_elements) {
        delete element;
    }
    UIState.clearUIElements();
}

void UISystem::RenderUIElements(SDL_Renderer* renderer, UIState& UIState, TTF_Font* UIFont)
{
    std::vector<UIElement*> ui_elements = UIState.getUIElements();
    for (auto& element : ui_elements) {
        element->renderElement(renderer, UIState, UIFont); // Pass UIFont if needed
    }
}   