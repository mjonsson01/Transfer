// File: Transfer/src/DynamoEngine/Debug/FrameProfiler.hpp

#pragma once

// SDL Imports
#include <SDL3/SDL_stdinc.h>

// Standard Library Imports
#include <cstddef>
#include <string>
#include <vector>

namespace DynamoEngine
{
// A rudimentary frame profiler: how long each named section of the frame takes, as an average and a "1% low" (the
// average of the slowest 1% of frames: the stutters an average hides). Times are collected for a window of frames,
// summarised into a report, then collected again from scratch.
//
// Use: addSection once per section (it gives back the section's number), then every frame record() each section's
// time and call endFrame(). After each full window, takeNewReport() returns true once and report() has the numbers.
class FrameProfiler
{
  public:
    struct SectionReport
    {
        std::string name;
        double average_ms = 0.0;
        double one_percent_low_ms = 0.0; // the average of the slowest 1% of frames (at least the single slowest)
        double worst_ms = 0.0;
        // How many frames of the window recorded this section at all. 0 = its code never ran (switched off, or not
        // reached), so its 0 ms is not a measurement; an overlay can leave such a section out.
        size_t frames_recorded = 0;
    };

    explicit FrameProfiler(size_t window_frames);

    // Returns the section's number, for record()
    int addSection(const std::string& name);
    // Adds to the section's time for the current frame: a section timed in several pieces records each piece
    void record(int section, double milliseconds);
    // Closes the frame: a section that recorded nothing this frame counts as 0 ms (and not as a recorded frame)
    void endFrame();

    // True once after each completed window (reading it clears it)
    bool takeNewReport();
    // The last completed window's numbers, one entry per section, in the order they were added
    const std::vector<SectionReport>& report() const { return m_report; }
    size_t windowFrames() const { return m_window_frames; }

    // Milliseconds between two SDL_GetPerformanceCounter() readings
    static double millisecondsBetween(Uint64 start, Uint64 end);

  private:
    struct Section
    {
        std::string name;
        double this_frame_ms = 0.0;
        bool recorded_this_frame = false;
        std::vector<double> frame_ms; // one entry per frame of the current window
        size_t frames_recorded = 0;   // frames of the current window that called record() for this section
    };

    void buildReport();

    size_t m_window_frames;
    size_t m_frames_in_window = 0;
    std::vector<Section> m_sections;
    std::vector<SectionReport> m_report;
    bool m_has_new_report = false;
};

// Times one block of code into a profiler section: the clock starts where the ProfileScope is declared, and the time
// is recorded when the enclosing { } block ends (C++ destroys the object there, and its destructor records it).
//     {
//         ProfileScope timer(profiler, gravity_section);
//         ...the code being timed...
//     } // <- recorded here
// With `enabled` false it records nothing (e.g. to time only some calls of a function used in several places).
class ProfileScope
{
  public:
    ProfileScope(FrameProfiler& profiler, int section, bool enabled = true);
    ~ProfileScope();

    // One timer per block: copying one would record the same time twice
    ProfileScope(const ProfileScope&) = delete;
    ProfileScope& operator=(const ProfileScope&) = delete;

  private:
    FrameProfiler& m_profiler;
    int m_section;
    bool m_enabled;
    Uint64 m_start;
};
} // namespace DynamoEngine
