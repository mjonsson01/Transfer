# Codebase Structure Report

## Header: `Core/UIState.hpp`
### Class: `UIState`
- `getMutableDEPRECATED_InputState()`
- `getDEPRECATED_InputState()`
- `getFPS()`
- `setFPS()`
- `getRenderDebug()`
- `setRenderDebug()`
- `getTimeScaleFactor()`
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
## Header: `Core/DEPRECATED_InputState.hpp`
### Class: `DEPRECATED_InputState`
- `resetTransientFlags()`
- `resetFlagsForSceneChange()`
- `clearAllBodies()`

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
- `liveParticleCount()`
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
- `survivableFragmentCount()`

---
## Header: `Systems/InputSystem.hpp`
### Class: `InputSystem`
- `processSystemInputFrame()`
- `cleanUp()`
- `trackDragAnchor()`
- `updateSceneUI()`
- `updateCamera()`
- `translateGameInputs()`
- `translateMenuInputs()`
- `copySharedPointerState()`
- `updateVisor()`

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
- `renderGameFrame()`
- `renderNonGameFrame()`
- `appendPreviewBodies()`
- `renderBodies()`
- `uploadUnifiedBodies()`
- `LoadShader()`
- `createUnifiedBodyGPUBufferAndPipeline()`
- `createUIGPUBufferAndPipeline()`
- `createVelocityVectorGPUBufferAndPipeline()`
- `createTwinklingStarGPUBufferAndPipeline()`
- `createStarshipGPUBufferAndPipeline()`
- `createFontAtlasSampler()`
- `rebuildFontAtlas()`
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
