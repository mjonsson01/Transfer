#include "UISystem.h"

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

void UISystem::RenderUIElements(SDL_Renderer* renderer, UIState& UIState)
{
    for (auto& element : elements) {
        element->renderElement(renderer);
    }
}   