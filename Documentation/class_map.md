# Codebase Structure Report

## Header: `Scenes/Scene.h`
### Class: `Scene`
- `populateMe()`
- `CleanUpSceneElements()`
- `getSceneElements()`

---
## Header: `Scenes/StartMenuScene/StartMenuScene.h`
### Class: `StartMenuScene`
- `populateMe()`

---
## Header: `Scenes/TestVisualScene/TestVisualScene.h`
### Class: `TestVisualScene`
- `populateMe()`

---
## Header: `Scenes/GameScene/GameScene.h`
### Class: `GameScene`
- `populateMe()`

---
## Header: `Scenes/PauseScene/PauseScene.h`
### Class: `PauseScene`
- `populateMe()`

---
## Header: `Core/InputState.h`
### Class: `InputState`
- `resetTransientFlags()`
- `resetFlagsForSceneChange()`
- `clearAllBodies()`

---
## Header: `Core/GameState.h`
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
- `getPlayMusic()`
- `invertPlayMusic()`
- `incrementMaxIDInstantiated()`
- `getMaxIDInstantiated()`

---
## Header: `Core/Game.h`
### Class: `Game`
- `StartGame()`
- `EndGame()`
- `Run()`
- `ProcessInput()`
- `IntegratePhysicsFrame()`
- `UpdateInstantiations()`
- `RenderFrame()`
- `PlayAudio()`
- `updateFPS()`
- `limitFrameRate()`

---
## Header: `Core/UIState.h`
### Class: `UIState`
- `getMutableInputState()`
- `getInputState()`
- `getFPS()`
- `setFPS()`
- `getAllUIVisibility()`
- `invertUIElementsVisibility()`
- `getRenderDebug()`
- `setRenderDebug()`
- `getTimeScaleFactor()`
- `getCurrentSceneID()`
- `setCurrentScene()`

---
## Header: `Utilities/UserInput/TransferInputs.h`
### Class: `TransferInputs`
- `operator=()`
- `resetAllInputsForSceneChange()`
- `resetAllMousePressedVars()`
- `resetJustPressed()`
- `resetAllKeyPressedVars()`

---
## Header: `Utilities/Math/Vector2D.h`
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
## Header: `Utilities/Rendering/Colors.h`
### Class: `ColorLibrary`
- *No methods found*

---
## Header: `Systems/InputSystem.h`
### Class: `InputSystem`
- `ProcessSystemInputFrame()`
- `CleanUp()`
- `routeSDL_EventInputInGame()`
- `routeSDL_EventInputInMenu()`
- `translateAndPassTransferInputsOff()`
- `translateAndPassMenuInputsOff()`

---
## Header: `Systems/RenderSystem.h`
### Class: `RenderSystem`
- `RenderFullFrame()`
- `CleanUp()`
- `getRenderer()`
- `getUIFontRegular()`
- `getUIFontTitle()`
- `renderGameFrame()`
- `renderNonGameFrame()`
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
## Header: `Systems/UISystem.h`
### Class: `UISystem`
- `CleanUp()`
- `UpdateUIElements()`
- `getScene()`
- `updateUISystemCurrentSceneID()`
- `updateGameUIElements()`
- `updateMenuUIElements()`
- `findElementWeAreIn()`
- `routeSliderInput()`
- `routeButtonClick()`
- `populateScenes()`
- `isSlider()`
- `isButton()`

---
## Header: `Systems/PhysicsSystem.h`
### Class: `SpawnLimiter`
- `canSpawn()`
- `reset()`

### Class: `PhysicsSystem`
- `UpdateSystemFrame()`
- `CleanUp()`
- `UpdateGravBodyInstantiations()`
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
## Header: `Systems/AudioSystem.h`
### Class: `AudioSystem`
- `ProcessSystemAudioFrame()`
- `LoadTrack()`
- `AddAllMusicToPlaylist()`
- `CleanUp()`

---
## Header: `Entities/Physics/GravitationalBody.h`
### Class: `GravitationalBody`
- *No methods found*

---
## Header: `Entities/VisualElements/TwinklingStars.h`
### Class: `TwinklingStar`
- *No methods found*

---
## Header: `Entities/UIElements/UIElement.h`
### Class: `UIElement`
- `renderMe()`
- `slideMe()`
- `clickMe()`
- `setPosition()`
- `getX()`
- `getY()`
- `setVisibility()`
- `checkAndReturnIfHit()`
- `getUIElementID()`

---
## Header: `Entities/UIElements/Buttons/Button.h`
### Class: `Button`
- `renderMe()`
- `getDisplayText()`
- `clickMe()`
- `getButtonState()`

---
## Header: `Entities/UIElements/Buttons/PlayGameButton/PlayGameButton.h`
### Class: `PlayGameButton`
- `clickMe()`

---
## Header: `Entities/UIElements/Buttons/ResumeButton/ResumeButton.h`
### Class: `ResumeButton`
- `clickMe()`

---
## Header: `Entities/UIElements/Checkboxes/Checkbox.h`
### Class: `Checkbox`
- *No methods found*

---
## Header: `Entities/UIElements/Overlay/FPSCounter.h`
### Class: `FPSCounter`
- `renderMe()`

---
## Header: `Entities/UIElements/Sliders/MassSlider.h`
### Class: `MassSlider`
- `getDisplayText()`
- `slideMe()`

---
## Header: `Entities/UIElements/Sliders/Slider.h`
### Class: `Slider`
- `renderMe()`
- `getDisplayText()`
- `slideMe()`
- `getSliderValue()`
- `getKnobPosition()`

---
## Header: `Entities/UIElements/Sliders/SimulationSpeedSlider.h`
### Class: `SimulationSpeedSlider`
- `getDisplayText()`

---
## Header: `Entities/UIElements/Sliders/RadiusSlider.h`
### Class: `RadiusSlider`
- `getDisplayText()`

---
