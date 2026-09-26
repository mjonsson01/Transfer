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
void FPSCounter::updateLayout(float windowWidth, float windowHeight)
{
    setPosition(windowWidth / SCREEN_WIDTH * 10, windowHeight / SCREEN_HEIGHT * 10);
}