// File: Transfer/src/Entities/UIElements/Overlay/FPSCounter.cpp

#include "Entities/UIElements/Overlay/FPSCounter.hpp"

FPSCounter::FPSCounter()
{
    updateLayout(SCREEN_WIDTH, SCREEN_HEIGHT);
    setVisibility(true);
    UIElementID = UIElementIdentifier::FPS_COUNTER_INDEX;
}

void FPSCounter::updateMe(UIState& uiState)
{
    float fps_local = uiState.getFPS();
    fps = static_cast<int>(fps_local);
}

void FPSCounter::buildGeometry(DynamoEngine::UIGeometryBuilder& builder)
{
    builder.addText("FPS: " + getDisplayText(), {getX(), getY()});
}
void FPSCounter::updateLayout(float window_width, float window_height)
{
    setPosition(window_width / SCREEN_WIDTH * 10, window_height / SCREEN_HEIGHT * 10);
}