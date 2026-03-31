// File: Transfer/src/Core/UIState.cpp

#include "UIState.h"

// UIState::UIState()
//     : fps(0.0f), showFPSCounter(true), UIElements() //instantiate the ui elements
// {
//     massKnobRect = {50.0f, 3*SCREEN_HEIGHT/4-10.0f, 20.0f, 40.0f}; // default position and size for the mass knob, will need to update based on resolution and slider position
// }

// UIState::~UIState()
// {
//     // Cleanup if necessary
// }

UIState::UIState()
    : framesPerSecond(TARGET_FPS)
{

}

UIState::~UIState()
{

}