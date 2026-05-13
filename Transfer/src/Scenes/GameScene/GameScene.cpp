// File: Transfer/src/Scenes/GameScene/GameScene.cpp

#include "GameScene.h"

GameScene::GameScene() : Scene(SceneIdentifier::GAME_SCENE)
{
    sceneUIElements.insert({UIElementIdentifier::FPS_COUNTER_INDEX, nullptr});
}

void GameScene::populateMe() {}