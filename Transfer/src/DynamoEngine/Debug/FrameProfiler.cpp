// File: Transfer/src/DynamoEngine/Debug/FrameProfiler.cpp
#include "DynamoEngine/Debug/FrameProfiler.hpp"

// SDL Imports
#include <SDL3/SDL_timer.h>

// Standard Library Imports
#include <algorithm>
#include <cassert>
#include <functional>

namespace DynamoEngine
{
FrameProfiler::FrameProfiler(size_t window_frames) : m_window_frames(window_frames)
{
    assert(window_frames > 0 && "a window needs at least one frame");
}

int FrameProfiler::addSection(const std::string& name)
{
    Section section;
    section.name = name;
    section.frame_ms.reserve(m_window_frames);
    m_sections.push_back(section);
    return static_cast<int>(m_sections.size()) - 1;
}

void FrameProfiler::record(int section, double milliseconds)
{
    assert(section >= 0 && section < static_cast<int>(m_sections.size()) && "record() with an unknown section");
    m_sections[static_cast<size_t>(section)].this_frame_ms += milliseconds;
}

void FrameProfiler::endFrame()
{
    for (Section& section : m_sections)
    {
        section.frame_ms.push_back(section.this_frame_ms);
        section.this_frame_ms = 0.0;
    }
    m_frames_in_window += 1;

    if (m_frames_in_window == m_window_frames)
    {
        buildReport();
        for (Section& section : m_sections)
        {
            section.frame_ms.clear(); // clear() keeps the memory: no allocation in the next window
        }
        m_frames_in_window = 0;
        m_has_new_report = true;
    }
}

bool FrameProfiler::takeNewReport()
{
    bool had_new_report = m_has_new_report;
    m_has_new_report = false;
    return had_new_report;
}

double FrameProfiler::millisecondsBetween(Uint64 start, Uint64 end)
{
    return static_cast<double>(end - start) * 1000.0 / static_cast<double>(SDL_GetPerformanceFrequency());
}

void FrameProfiler::buildReport()
{
    m_report.clear();
    for (const Section& section : m_sections)
    {
        SectionReport section_report;
        section_report.name = section.name;
        if (!section.frame_ms.empty())
        {
            // Slowest first, so the slowest 1% are simply the first entries
            std::vector<double> slowest_first = section.frame_ms;
            std::sort(slowest_first.begin(), slowest_first.end(), std::greater<double>());

            double total = 0.0;
            for (double ms : slowest_first)
            {
                total += ms;
            }
            size_t slowest_count = std::max<size_t>(1, slowest_first.size() / 100);
            double slowest_total = 0.0;
            for (size_t i = 0; i < slowest_count; ++i)
            {
                slowest_total += slowest_first[i];
            }

            section_report.average_ms = total / static_cast<double>(slowest_first.size());
            section_report.one_percent_low_ms = slowest_total / static_cast<double>(slowest_count);
            section_report.worst_ms = slowest_first[0];
        }
        m_report.push_back(section_report);
    }
}

ProfileScope::ProfileScope(FrameProfiler& profiler, int section, bool enabled)
    : m_profiler(profiler), m_section(section), m_enabled(enabled), m_start(SDL_GetPerformanceCounter())
{
}

ProfileScope::~ProfileScope()
{
    if (m_enabled)
    {
        m_profiler.record(m_section, FrameProfiler::millisecondsBetween(m_start, SDL_GetPerformanceCounter()));
    }
}
} // namespace DynamoEngine
