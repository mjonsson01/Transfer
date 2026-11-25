# Codebase Structure Report

## Header: `Core/UIState.hpp`
### Class: `UIState`
- `getFPS()`
- `setFPS()`
- `getRenderDebug()`
- `setRenderDebug()`
- `getTimeScaleFactor()`
- `setSimulationSpeed()`
- `spawnSettings()`
- `spawnSettings()`
- `QueueSoundEffect()`
- `HasPendingSoundEffects()`
- `PopNextSoundEffect()`
- `getRequestedMusicMode()`
- `setRequestedMusicMode()`
- `getPlayMusic()`
- `setPlayMusic()`
- `getPlaySoundEffects()`
- `setPlaySoundEffects()`
- `getShowTrajectoryTimeMarks()`
- `setShowTrajectoryTimeMarks()`
- `getShowProfiler()`
- `setShowProfiler()`
- `profilerLines()`
- `setProfilerLines()`
- `getIsPlanning()`
- `setIsPlanning()`

---
## Header: `Core/GameState.hpp`
### Class: `GameState`
- `IsPlaying()`
- `SetPlaying()`
- `getIsShuttingDownAudioSystem()`
- `setIsShuttingDownAudioSystem()`
- `getAlpha()`
- `setAlpha()`
- `getCameraStateMutable()`
- `getCameraState()`
- `getPlayerMutable()`
- `getPlayer()`
- `getUniverseStateMutable()`
- `getUniverseState()`
- `getTrajectoriesMutable()`
- `getTrajectories()`

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
- `publishProfilerReport()`

---
## Header: `Camera/CameraState.hpp`
### Class: `CameraState`
- *No methods found*

---
## Header: `Systems/Input/PlayerInput.hpp`
### Class: `ShipControls`
- *No methods found*

### Class: `SpawnProperties`
- *No methods found*

### Class: `PlayerCommands`
- *No methods found*

### Class: `PlayerInput`
- *No methods found*

---
## Header: `Systems/Input/InputSystem.hpp`
### Class: `InputSystem`
- `processSystemInputFrame()`
- `cleanUp()`
- `trackDragAnchor()`
- `updateSceneUI()`
- `updateCamera()`
- `translateGameInputs()`
- `translateMenuInputs()`
- `aimSpawn()`
- `updateVisor()`

---
## Header: `Systems/Audio/AudioSystem.hpp`
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
- `cleanUp()`
- `processMusic()`

---
## Header: `Systems/Audio/SoundEffect.hpp`
### Class: `SoundEffect`
- *No methods found*

---
## Header: `Systems/Rendering/CameraData.hpp`
### Class: `CameraConstants`
- *No methods found*

---
## Header: `Systems/Rendering/TwinklingStars.hpp`
### Class: `TwinklingStar`
- *No methods found*

---
## Header: `Systems/Rendering/GPUTypes.hpp`
### Class: `UnifiedBodyVertex`
- *No methods found*

### Class: `TwinklingStarVertex`
- *No methods found*

### Class: `VelocityVectorVertex`
- *No methods found*

### Class: `TrajectoryLineVertex`
- *No methods found*

### Class: `DebugLineVertex`
- *No methods found*

### Class: `StarshipVertex`
- *No methods found*

---
## Header: `Systems/Rendering/RenderSystem.hpp`
### Class: `RenderSystem`
- `renderFullFrame()`
- `cleanUp()`
- `getUIFontRegular()`
- `renderGameFrame()`
- `renderNonGameFrame()`
- `appendPreviewBodies()`
- `renderBodies()`
- `uploadUnifiedBodies()`
- `createUnifiedBodyBuffers()`
- `LoadShader()`
- `createUnifiedBodyGPUBufferAndPipeline()`
- `createUIGPUBufferAndPipeline()`
- `createVelocityVectorGPUBufferAndPipeline()`
- `createTwinklingStarGPUBufferAndPipeline()`
- `createStarshipGPUBufferAndPipeline()`
- `createDebugLineGPUBufferAndPipeline()`
- `buildDebugLines()`
- `addDebugLine()`
- `addDebugCircle()`
- `uploadDebugLines()`
- `renderDebugLines()`
- `createTrajectoryLineBuffersAndPipeline()`
- `buildTrajectoryLines()`
- `appendTrajectorySegments()`
- `appendTrajectoryTimeMarks()`
- `appendTrajectoryTip()`
- `uploadTrajectoryLines()`
- `renderTrajectoryLines()`
- `createFontAtlasSampler()`
- `loadSpriteTexture()`
- `createSpriteSampler()`
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
## Header: `Systems/Physics/GravitationalBody.hpp`
### Class: `GravitationalBody`
- `toUnifiedVertex()`

---
## Header: `Systems/Physics/UniverseState.hpp`
### Class: `UniverseState`
- *No methods found*

---
## Header: `Systems/Physics/UniformParticleGrid.hpp`
### Class: `UniformParticleGrid`
- `build()`
- `queryCandidates()`
- `firstEntryAtOrAfter()`

---
## Header: `Systems/Physics/PhysicsSystem.hpp`
### Class: `CollisionInfo`
- *No methods found*

### Class: `PhysicsSystem`
- `advanceUniverseState()`
- `cleanUp()`
- `updateGravBodyInstantiations()`
- `adjustThrottleWhilePlanning()`
- `updateTrajectories()`
- `integrateUniversePhysics()`
- `removeNewestSpawn()`
- `handleCollisionEvents()`
- `handleMacroMacroCollisions()`
- `handleMacroParticleCollisions()`
- `handleParticleParticleCollisions()`
- `handleShipCollisions()`
- `resolveShipContact()`
- `handleDynamicCollision()`
- `handleElasticCollisions()`
- `handleAccretion()`
- `promoteOversizedParticles()`
- `substituteWithParticles()`
- `substituteWithParticlesFromImpact()`
- `liveParticleCount()`
- `potentialParticleCount()`
- `hasRoomForParticles()`
- `recomputeAllForces()`
- `addGravityBetweenGravBodies()`
- `addGravityBetweenGravBodyPair()`
- `addGravityBetweenShipAndBodies()`
- `createMacroBody()`
- `createParticleCluster()`
- `applyShipInputPhysics()`
- `cleanupGravBodies()`
- `cleanupParticles()`
- `cleanupMacroBodies()`
- `startTrajectories()`
- `clearTrajectories()`
- `universeChangedSinceTraceStart()`
- `advanceIntegratedTrajectories()`
- `recordTrajectoryPoints()`
- `updateSpawnPreviewTrajectory()`
- `clearSpawnPreviewTrajectory()`
- `makeSpawnPreviewBody()`

---
## Header: `Systems/Physics/GravitationalBodyPair.hpp`
### Class: `GravitationalBodyPair`
- *No methods found*

---
## Header: `Systems/Physics/Trajectories.hpp`
### Class: `TracedPoint`
- *No methods found*

### Class: `Trajectories`
- *No methods found*

---
## Header: `Systems/Physics/SpawnSettings.hpp`
### Class: `SpawnSettings`
- `setFlag()`
- `flag()`

---
## Header: `Player/Starship.hpp`
### Class: `Starship`
- `applyRotation()`
- `applyThrust()`
- `adjustThrottle()`
- `throttle()`
- `applyVelocityVerletFirstHalfStep()`
- `setGravityForce()`
- `applyVelocityVerletSecondHalfStep()`
- `center()`
- `halfSize()`
- `mass()`
- `velocity()`
- `applyImpulse()`
- `moveBy()`
- `boundingRadius()`
- `setLastImpactEnergy()`
- `lastImpactEnergy()`
- `buildGeometry()`
- `getPointingVector()`
- `collisionPolygon()`
- `movementThisTick()`
- `rotatedOffset()`

---
## Header: `Player/Player.hpp`
### Class: `Player`
- *No methods found*

---
