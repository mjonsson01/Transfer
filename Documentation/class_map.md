# Codebase Structure Report

## Header: `Core\Game.h`
### Class: `Game`
- `StartGame()`
- `EndGame()`
- `Run()`
- `ProcessInput()`
- `UpdatePhysicsFrame()`
- `RenderFrame()`
- `PlayAudio()`
- `UpdateFPS()`
- `LimitFrameRate()`

---
## Header: `Core\GameState.h`
### Class: `GameState`
- `IsPlaying()`
- `SetPlaying()`
- `getIsShuttingDownAudioSystem()`
- `setIsShuttingDownAudioSystem()`
- `getParticles()`
- `getParticlesMutable()`
- `getMacroBodies()`
- `getMacroBodiesMutable()`
- `getAlpha()`
- `setAlpha()`
- `getTimeScaleFactor()`
- `getToggleSlow()`
- `invertToggleSlow()`
- `getToggleFast()`
- `invertToggleFast()`
- `getPlayMusic()`
- `invertPlayMusic()`
- `incrementMaxIDInstantiated()`
- `getMaxIDInstantiated()`

---
## Header: `Core\InputState.h`
### Class: `InputState`
- `resetTransientFlags()`
- `togglePhysicsPause()`
- `clearAllBodies()`

---
## Header: `Core\UIState.h`
### Class: `UIState`
- `getMutableInputState()`
- `getInputState()`
- `getFPS()`
- `setFPS()`
- `getAllUIVisibility()`
- `invertUIElementsVisibility()`
- `getGameScene()`
- `setGameScene()`
- `getLevelEditorScene()`
- `setLevelEditorScene()`
- `getLevelSelectScene()`
- `setLevelSelectScene()`
- `getPauseMenuActive()`
- `setPauseMenuActive()`
- `getStartMenuActive()`
- `setStartMenuActive()`
- `getRenderDebug()`
- `setRenderDebug()`

---
## Header: `Entities\Physics\GravitationalBody.h`
### Class: `GravitationalBody`
- *No methods found*

---
## Header: `Entities\UIElements\UIElement.h`
### Class: `UIElement`
- `renderMe()`
- `slideMe()`
- `clickMe()`
- `setPosition()`
- `getX()`
- `getY()`
- `setVisibility()`
- `checkAndReturnIfHit()`
- `getUIElementIdentifier()`

---
## Header: `Entities\UIElements\Buttons\Button.h`
### Class: `Button`
- `renderMe()`
- `getDisplayText()`
- `clickMe()`
- `getButtonState()`

---
## Header: `Entities\UIElements\Buttons\PlayGameButton.h`
### Class: `PlayGameButton`
- *No methods found*

---
## Header: `Entities\UIElements\Checkboxes\Checkbox.h`
### Class: `Checkbox`
- *No methods found*

---
## Header: `Entities\UIElements\Overlay\FPSCounter.h`
### Class: `FPSCounter`
- `renderMe()`

---
## Header: `Entities\UIElements\Sliders\MassSlider.h`
### Class: `MassSlider`
- `getDisplayText()`

---
## Header: `Entities\UIElements\Sliders\RadiusSlider.h`
### Class: `RadiusSlider`
- `getDisplayText()`

---
## Header: `Entities\UIElements\Sliders\Slider.h`
### Class: `Slider`
- `renderMe()`
- `getDisplayText()`
- `slideMe()`
- `getSliderValue()`

---
## Header: `Entities\VisualElements\TwinklingStars.h`
### Class: `TwinklingStar`
- *No methods found*

---
## Header: `Scenes\Scene.h`
### Class: `Scene`
- *No methods found*

---
## Header: `Systems\AudioSystem.h`
### Class: `AudioSystem`
- `ProcessSystemAudioFrame()`
- `LoadTrack()`
- `AddAllMusicToPlaylist()`
- `CleanUp()`

---
## Header: `Systems\InputSystem.h`
### Class: `InputSystem`
- `ProcessSystemInputFrame()`
- `CleanUp()`
- `routeSDL_EventInputInGame()`
- `routeSDL_EventInputInMenu()`
- `translateAndPassTransferInputsOff()`
- `translateAndPassMenuInputsOff()`

---
## Header: `Systems\PhysicsSystem.h`
### Class: `SpawnLimiter`
- `canSpawn()`
- `reset()`

### Class: `PhysicsSystem`
- `UpdateSystemFrame()`
- `CleanUp()`
- `updateGravBodyInstantiations()`
- `updateGravityForSystem()`
- `calculateGravity()`
- `integrateForwardsPhase1()`
- `integrateForwardsPhase2()`
- `handleCollisions()`
- `handleElasticCollisions()`
- `handleDynamicExplosionCollision()`
- `handleAccretion()`
- `createMacroBody()`
- `createParticle()`
- `createParticleCluster()`
- `calculateTotalEnergy()`
- `substituteWithParticles()`
- `populateCollisionProxyFromMacroBody()`
- `cleanupParticles()`
- `cleanupMacroBodies()`

---
## Header: `Systems\RenderSystem.h`
### Class: `RenderSystem`
- `RenderFullFrame()`
- `CleanUp()`
- `getRenderer()`
- `getUIFontRegular()`
- `getUIFontTitle()`
- `renderPreviewBodies()`
- `renderBodies()`
- `renderDragLine()`
- `renderUIElements()`
- `getColorForProperty()`
- `buildCircleTextureCache()`
- `clearCachedCircleTextures()`
- `createStarField()`
- `createStarTextures()`
- `updateStars()`
- `renderStars()`

---
## Header: `Systems\UISystem.h`
### Class: `UISystem`
- `CleanUp()`
- `UpdateUIElements()`
- `updateGameUIElements()`
- `getGameUIElementsMutable()`
- `getGameUIElements()`
- `updatePauseUIElements()`
- `getPauseUIElementsMutable()`
- `getPauseUIElements()`
- `findElementWeAreIn()`
- `routeSliderInput()`
- `routeButtonClick()`
- `isSlider()`
- `isButton()`

---
## Header: `Utilities\Math\Vector2D.h`
### Class: `Vector2D`
- `operator+()`
- `operator-()`
- `operator*()`
- `operator/()`
- `operator+=()`
- `operator-=()`
- `operator*=()`
- `operator/=()`
- `magnitude()`
- `square_magnitude()`
- `dot()`
- `normalizeInPlace()`
- `normalize()`

---
## Header: `Utilities\Rendering\Colors.h`
### Class: `ColorLibrary`
- *No methods found*

---
## Header: `Utilities\UserInput\TransferInputs.h`
### Class: `TransferInputs`
- `operator=()`
- `resetAllMousePressedVars()`
- `resetJustPressed()`
- `resetAllKeyPressedVars()`

---
