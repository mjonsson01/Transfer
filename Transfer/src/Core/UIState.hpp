// File: Transfer/src/Core/UIState.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>

// Custom Imports
#include "DynamoEngine/Scenes/Scene.hpp"
#include "Systems/Physics/SpawnSettings.hpp"
#include "Systems/Audio/MusicModeEnum.hpp"
#include "Constants/GameSystemConstants.hpp"

// Standard Library Imports
#include <queue>
#include <string>
#include <vector>

// Which profiler the overlay shows (only one at a time: both together wouldn't fit on screen)
enum class ProfilerOverlay
{
    Off,
    Frame,   // F4
    Physics, // F5
};

class UIState
{
  public:
    UIState();
    ~UIState();
    float getFPS() { return framesPerSecond; }
    void setFPS(float fps) { framesPerSecond = fps; }
    bool getRenderDebug() { return renderDebug; }
    void setRenderDebug(bool rd) { renderDebug = rd; }
    // The simulation speed slider (1 = normal speed). Planning mode stops time without changing it.
    float getTimeScaleFactor() const { return static_cast<float>(m_simulation_speed); }
    void setSimulationSpeed(double speed) { m_simulation_speed = speed; }
    // How fast the simulation really runs (game seconds per real second, from the physics thread). Lower than the
    // slider only while the machine can't keep up (see MAX_PHYSICS_BACKLOG_TICKS); otherwise equal to it.
    double getEffectiveSimulationSpeed() const { return m_effective_simulation_speed; }
    void setEffectiveSimulationSpeed(double speed) { m_effective_simulation_speed = speed; }
    // What the next spawned planet or cluster will be like (edited by the spawn panel's checkboxes)
    SpawnSettings& spawnSettings() { return m_spawn_settings; }
    const SpawnSettings& spawnSettings() const { return m_spawn_settings; }

    void QueueSoundEffect(const std::string& soundName) { pendingSoundEffects.push(soundName); };
    bool HasPendingSoundEffects() const { return !pendingSoundEffects.empty(); };

    std::string PopNextSoundEffect()
    {
        if (pendingSoundEffects.empty())
            return "";

        std::string sound = pendingSoundEffects.front();
        pendingSoundEffects.pop();

        return sound;
    }
    MusicMode getRequestedMusicMode() const { return requestedMusicMode; }
    void setRequestedMusicMode(MusicMode musicMode) { requestedMusicMode = musicMode; }
    bool getPlayMusic() const { return playMusic; }
    void setPlayMusic(bool spm) { playMusic = spm; }
    bool getPlaySoundEffects() const { return playSoundEffects; }
    void setPlaySoundEffects(bool pse) { playSoundEffects = pse; }

    // Time marks on the trajectory lines (every TRAJECTORY_TIME_MARK_SECONDS); will become a setting
    bool getShowTrajectoryTimeMarks() const { return m_show_trajectory_time_marks; }
    void setShowTrajectoryTimeMarks(bool show) { m_show_trajectory_time_marks = show; }
    // The profiler overlay: off, F4 (the whole frame, system by system) or F5 (the physics tick, function by
    // function). Game writes each one's lines after every profiler window.
    ProfilerOverlay getProfilerOverlay() const { return m_profiler_overlay; }
    void setProfilerOverlay(ProfilerOverlay overlay) { m_profiler_overlay = overlay; }
    const std::vector<std::string>& frameProfilerLines() const { return m_frame_profiler_lines; }
    void setFrameProfilerLines(const std::vector<std::string>& lines) { m_frame_profiler_lines = lines; }
    const std::vector<std::string>& physicsProfilerLines() const { return m_physics_profiler_lines; }
    void setPhysicsProfilerLines(const std::vector<std::string>& lines) { m_physics_profiler_lines = lines; }
    // The scene Esc paused (Resume and Esc in the pause menu go back to it)
    DynamoEngine::SceneId getResumeScene() const { return m_resume_scene; }
    void setResumeScene(DynamoEngine::SceneId scene) { m_resume_scene = scene; }
    bool getIsPlanning() const { return m_is_planning; }
    void setIsPlanning(bool is_planning) { m_is_planning = is_planning; }

  private:
    double m_simulation_speed = 1.0;
    double m_effective_simulation_speed = 1.0;
    SpawnSettings m_spawn_settings;
    float framesPerSecond = TARGET_FPS;

    bool renderDebug = VIEW_DEBUG; // Toggles rendering of debug elements like
    // collision boxes, spawn areas, etc.

    // SOUND STUFF
    bool playMusic = false;
    bool playSoundEffects = false;
    MusicMode requestedMusicMode = MusicMode::NONE;
    std::queue<std::string> pendingSoundEffects;

    bool m_show_trajectory_time_marks = true;
    ProfilerOverlay m_profiler_overlay = ProfilerOverlay::Off;
    std::vector<std::string> m_frame_profiler_lines;
    std::vector<std::string> m_physics_profiler_lines;
    DynamoEngine::SceneId m_resume_scene = 0; // set every time a simulation scene is paused
    bool m_is_planning = false; // P: time stops, and I can trace trajectories
};
