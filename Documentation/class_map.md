# Codebase Structure Report

## Header: `Scenes/Scene.hpp`
### Class: `Scene`
- `populateMe()`
- `CleanUpSceneElements()`
- `getSceneElements()`

---
## Header: `Scenes/StartMenuScene/StartMenuScene.hpp`
### Class: `StartMenuScene`
- `populateMe()`

---
## Header: `Scenes/TestVisualScene/TestVisualScene.hpp`
### Class: `TestVisualScene`
- `populateMe()`

---
## Header: `Scenes/GameScene/GameScene.hpp`
### Class: `GameScene`
- `populateMe()`

---
## Header: `Scenes/PauseScene/PauseScene.hpp`
### Class: `PauseScene`
- `populateMe()`

---
## Header: `Core/UIState.hpp`
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
- `QueueSoundEffect()`
- `HasPendingSoundEffects()`
- `PopNextSoundEffect()`
- `getRequestedMusicMode()`
- `setRequestedMusicMode()`
- `getPlayMusic()`
- `setPlayMusic()`
- `getPlaySoundEffects()`
- `setPlaySoundEffects()`

---
## Header: `Core/GameState.hpp`
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
- `incrementMaxIDInstantiated()`
- `getMaxIDInstantiated()`
- `getCameraStateMutable()`
- `getCameraState()`
- `getPlayerMutable()`
- `getPlayer()`

---
## Header: `Core/CameraState.hpp`
### Class: `CameraState`
- *No methods found*

---
## Header: `Core/Game.hpp`
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
## Header: `Core/InputState.hpp`
### Class: `InputState`
- `resetTransientFlags()`
- `resetFlagsForSceneChange()`
- `clearAllBodies()`

---
## Header: `Utilities/UserInput/TransferInputs.hpp`
### Class: `TransferInputs`
- `operator=()`
- `resetAllInputsForSceneChange()`
- `resetAllMousePressedVars()`
- `resetJustPressed()`
- `resetAllKeyPressedVars()`

---
## Header: `Utilities/Math/Vector2D.hpp`
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
## Header: `Utilities/Rendering/Colors.hpp`
### Class: `ColorLibrary`
- *No methods found*

---
## Header: `Utilities/Rendering/FontAtlasUtility.hpp`
### Class: `GlyphMetrics`
- *No methods found*

### Class: `FontAtlasUtility`
- `BuildAtlas()`
- `GetGlyph()`
- `CalculateTextWidth()`
- `GetFontHeight()`

---
## Header: `Utilities/Rendering/CameraData.hpp`
### Class: `CameraConstants`
- *No methods found*

---
## Header: `Utilities/Rendering/GPUTypes.hpp`
### Class: `UnifiedBodyVertex`
- *No methods found*

### Class: `TwinklingStarVertex`
- *No methods found*

### Class: `UIElementVertex`
- *No methods found*

### Class: `VelocityVectorVertex`
- *No methods found*

### Class: `StarshipVertex`
- *No methods found*

---
## Header: `Utilities/Physics/UniformParticleGrid.hpp`
### Class: `UniformParticleGrid`
- `build()`
- `queryCandidates()`

---
## Header: `Systems/PhysicsSystem.hpp`
### Class: `CollisionInfo`
- *No methods found*

### Class: `PhysicsSystem`
- `UpdateSystemFrame()`
- `CleanUp()`
- `UpdateGravBodyInstantiations()`
- `handleCollisions()`
- `handleMacroMacroCollisions()`
- `handleMacroParticleCollisions()`
- `handleParticleParticleCollisions()`
- `handleDynamicCollision()`
- `handleElasticCollisions()`
- `handleAccretion()`
- `promoteOversizedParticles()`
- `substituteWithParticles()`
- `substituteWithParticlesFromImpact()`
- `updateAllForces()`
- `updateGravityForSystem()`
- `calculateGravity()`
- `integrateForwardsVelocityVerletPhase1()`
- `applyVelocityVerletPhase1()`
- `integrateForwardsVelocityVerletPhase2()`
- `applyVelocityVerletPhase2()`
- `createMacroBody()`
- `createParticle()`
- `createParticleCluster()`
- `calculateTotalEnergy()`
- `updatePlayerPhysics()`
- `cleanupParticles()`
- `cleanupMacroBodies()`

---
## Header: `Systems/UISystem.hpp`
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
- `updateAllUILayouts()`
- `isSlider()`
- `isButton()`

---
## Header: `Systems/InputSystem.hpp`
### Class: `InputSystem`
- `ProcessSystemInputFrame()`
- `CleanUp()`
- `routeSDL_EventInputInGame()`
- `routeSDL_EventInputInMenu()`
- `translateAndPassTransferInputsOff()`
- `translateAndPassMenuInputsOff()`

---
## Header: `Systems/AudioSystem.hpp`
### Class: `AudioSystem`
- `ProcessSystemAudioFrame()`
- `loadAndPlayTrack()`
- `addAllSoundEffectsToLibrary()`
- `loadMusicLibrary()`
- `onTrackFinished()`
- `playSoundEffect()`
- `transitionMusicMode()`
- `stopCurrentMusic()`
- `prepareShufflePlaylist()`
- `playNextShuffleTrack()`
- `cleanupFinishedSFX()`
- `CleanUp()`
- `processMusic()`

---
## Header: `Systems/RenderSystem.hpp`
### Class: `RenderSystem`
- `RenderFullFrame()`
- `CleanUp()`
- `getUIFontRegular()`
- `getUIFontTitle()`
- `renderGameFrame()`
- `renderNonGameFrame()`
- `renderTestFrame()`
- `appendPreviewBodies()`
- `renderBodies()`
- `uploadUnifiedBodies()`
- `LoadShader()`
- `createUnifiedBodyGPUBufferAndPipeline()`
- `createUIGPUBufferAndPipeline()`
- `createVelocityVectorGPUBufferAndPipeline()`
- `createTwinklingStarGPUBufferAndPipeline()`
- `createStarshipGPUBufferAndPipeline()`
- `createFontAtlasTextureAndSampler()`
- `uploadUIVertices()`
- `renderUIElements()`
- `createTwinklingStarField()`
- `uploadTwinklingStarField()`
- `uploadStarship()`
- `renderTwinklingStarField()`
- `renderStarship()`
- `buildCameraConstants()`
- `buildVelocityVectorGeometry()`
- `uploadVelocityVectorVertices()`
- `renderVelocityVectors()`
- `getColorForProperty()`

---
## Header: `Entities/Physics/GravitationalBody.hpp`
### Class: `GravitationalBody`
- `toUnifiedVertex()`

---
## Header: `Entities/Physics/GravitationalBodyPair.hpp`
### Class: `GravitationalBodyPair`
- *No methods found*

---
## Header: `Entities/VisualElements/TwinklingStars.hpp`
### Class: `TwinklingStar`
- *No methods found*

---
## Header: `Entities/Sound/SoundEffect.hpp`
### Class: `SoundEffect`
- *No methods found*

---
## Header: `Entities/UIElements/UIElement.hpp`
### Class: `UIElement`
- `slideMe()`
- `clickMe()`
- `setPosition()`
- `getX()`
- `getY()`
- `setVisibility()`
- `isVisible()`
- `checkAndReturnIfHit()`
- `getUIElementID()`
- `buildGeometry()`
- `updateMe()`
- `updateLayout()`

---
## Header: `Entities/UIElements/Buttons/Button.hpp`
### Class: `Button`
- `buildGeometry()`
- `getDisplayText()`
- `clickMe()`
- `getButtonState()`

---
## Header: `Entities/UIElements/Buttons/VisorButton/VisorButton.hpp`
### Class: `VisorButton`
- *No methods found*

---
## Header: `Entities/UIElements/Buttons/PlayGameButton/PlayGameButton.hpp`
### Class: `PlayGameButton`
- `clickMe()`
- `updateLayout()`

---
## Header: `Entities/UIElements/Buttons/ResumeButton/ResumeButton.hpp`
### Class: `ResumeButton`
- `clickMe()`
- `updateLayout()`

---
## Header: `Entities/UIElements/DropDownMenu/DropDownMenu.hpp`
### Class: `DropDownMenu`
- `buildGeometry()`

---
## Header: `Entities/UIElements/Checkboxes/Checkbox.hpp`
### Class: `Checkbox`
- *No methods found*

---
## Header: `Entities/UIElements/Overlay/FPSCounter.hpp`
### Class: `FPSCounter`
- `buildGeometry()`
- `updateMe()`
- `updateLayout()`
- `getDisplayText()`

---
## Header: `Entities/UIElements/Sliders/Slider.hpp`
### Class: `Slider`
- `buildGeometry()`
- `getDisplayText()`
- `slideMe()`
- `getSliderValue()`
- `getKnobPosition()`
- `playTickSoundIfMoved()`

---
## Header: `Entities/UIElements/Sliders/RadiusSlider.hpp`
### Class: `RadiusSlider`
- `getDisplayText()`
- `updateLayout()`

---
## Header: `Entities/UIElements/Sliders/MassSlider.hpp`
### Class: `MassSlider`
- `getDisplayText()`
- `slideMe()`
- `updateLayout()`
- `playTickSoundIfMoved()`

---
## Header: `Entities/UIElements/Sliders/SimulationSpeedSlider.hpp`
### Class: `SimulationSpeedSlider`
- `getDisplayText()`
- `updateLayout()`

---
## Header: `Player/Starship.hpp`
### Class: `Starship`
- `integratePosition()`
- `applyVelocity()`
- `applyRotation()`
- `buildGeometry()`
- `getPointingVector()`

---
## Header: `Player/Player.hpp`
### Class: `Player`
- *No methods found*

---
