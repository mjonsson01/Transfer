// File: Transfer/src/Core/UIState.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>

// Custom Imports
#include "Systems/Physics/SpawnSettings.hpp"
#include "Systems/Audio/MusicModeEnum.hpp"
#include "Constants/GameSystemConstants.hpp"

// Standard Library Imports
#include <queue>
#include <string>
#include <vector>

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
    // F4 profiler overlay: on/off, and the lines it shows (Game writes them after every profiler window)
    bool getShowProfiler() const { return m_show_profiler; }
    void setShowProfiler(bool show) { m_show_profiler = show; }
    const std::vector<std::string>& profilerLines() const { return m_profiler_lines; }
    void setProfilerLines(const std::vector<std::string>& lines) { m_profiler_lines = lines; }
    bool getIsPlanning() const { return m_is_planning; }
    void setIsPlanning(bool is_planning) { m_is_planning = is_planning; }

  private:
    double m_simulation_speed = 1.0;
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
    bool m_show_profiler = false;
    std::vector<std::string> m_profiler_lines;
    bool m_is_planning = false; // P: time stops, and I can trace trajectories
};
