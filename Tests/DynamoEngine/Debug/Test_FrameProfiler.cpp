// File: Tests/DynamoEngine/Debug/Test_FrameProfiler.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/Debug/FrameProfiler.hpp"

using namespace DynamoEngine;

// --------- TESTS --------- //

TEST(FrameProfiler, NoReportBeforeTheWindowIsFull)
{
    FrameProfiler profiler(10);
    int section = profiler.addSection("Physics");
    for (int frame = 0; frame < 9; ++frame)
    {
        profiler.record(section, 1.0);
        profiler.endFrame();
    }
    EXPECT_FALSE(profiler.takeNewReport());
}

TEST(FrameProfiler, ReportsOnceWhenTheWindowIsFull)
{
    FrameProfiler profiler(10);
    int section = profiler.addSection("Physics");
    for (int frame = 0; frame < 10; ++frame)
    {
        profiler.record(section, 1.0);
        profiler.endFrame();
    }
    EXPECT_TRUE(profiler.takeNewReport());
    EXPECT_FALSE(profiler.takeNewReport()); // reading it cleared it
}

TEST(FrameProfiler, AverageIsTheMeanOfTheFrames)
{
    FrameProfiler profiler(4);
    int section = profiler.addSection("Render");
    for (double ms : {1.0, 2.0, 3.0, 6.0})
    {
        profiler.record(section, ms);
        profiler.endFrame();
    }
    ASSERT_TRUE(profiler.takeNewReport());
    EXPECT_DOUBLE_EQ(profiler.report()[0].average_ms, 3.0);
    EXPECT_DOUBLE_EQ(profiler.report()[0].worst_ms, 6.0);
}

TEST(FrameProfiler, OnePercentLowIsTheMeanOfTheSlowestOnePercent)
{
    // 200 frames: the slowest 1% is the slowest 2 frames
    FrameProfiler profiler(200);
    int section = profiler.addSection("Physics");
    for (int frame = 0; frame < 200; ++frame)
    {
        double ms = 1.0;
        if (frame == 50)
            ms = 10.0;
        if (frame == 150)
            ms = 20.0;
        profiler.record(section, ms);
        profiler.endFrame();
    }
    ASSERT_TRUE(profiler.takeNewReport());
    EXPECT_DOUBLE_EQ(profiler.report()[0].one_percent_low_ms, 15.0);
}

TEST(FrameProfiler, OnePercentLowOfAShortWindowIsTheSlowestFrame)
{
    // Fewer than 100 frames: 1% rounds down to 0 frames, so it uses the single slowest one
    FrameProfiler profiler(5);
    int section = profiler.addSection("Audio");
    for (double ms : {1.0, 9.0, 2.0, 3.0, 1.0})
    {
        profiler.record(section, ms);
        profiler.endFrame();
    }
    ASSERT_TRUE(profiler.takeNewReport());
    EXPECT_DOUBLE_EQ(profiler.report()[0].one_percent_low_ms, 9.0);
}

TEST(FrameProfiler, RecordsInOneFrameAddUp)
{
    // A section timed in pieces (e.g. several physics ticks in one frame) counts as their sum for that frame
    FrameProfiler profiler(1);
    int section = profiler.addSection("Physics");
    profiler.record(section, 1.5);
    profiler.record(section, 2.5);
    profiler.endFrame();
    ASSERT_TRUE(profiler.takeNewReport());
    EXPECT_DOUBLE_EQ(profiler.report()[0].average_ms, 4.0);
}

TEST(FrameProfiler, SectionsAreIndependentAndKeepTheirOrder)
{
    FrameProfiler profiler(2);
    int input = profiler.addSection("Input");
    int render = profiler.addSection("Render");
    for (int frame = 0; frame < 2; ++frame)
    {
        profiler.record(input, 1.0);
        profiler.record(render, 5.0);
        profiler.endFrame();
    }
    ASSERT_TRUE(profiler.takeNewReport());
    ASSERT_EQ(profiler.report().size(), 2u);
    EXPECT_EQ(profiler.report()[0].name, "Input");
    EXPECT_DOUBLE_EQ(profiler.report()[0].average_ms, 1.0);
    EXPECT_EQ(profiler.report()[1].name, "Render");
    EXPECT_DOUBLE_EQ(profiler.report()[1].average_ms, 5.0);
}

TEST(FrameProfiler, ASectionWithoutARecordCountsAsZero)
{
    // Skipped work (e.g. no physics tick this frame) is a real 0 ms frame, not a missing one
    FrameProfiler profiler(2);
    int section = profiler.addSection("Physics");
    profiler.record(section, 4.0);
    profiler.endFrame();
    profiler.endFrame();
    ASSERT_TRUE(profiler.takeNewReport());
    EXPECT_DOUBLE_EQ(profiler.report()[0].average_ms, 2.0);
}

TEST(FrameProfiler, EachWindowStartsFromScratch)
{
    FrameProfiler profiler(2);
    int section = profiler.addSection("Physics");
    for (double ms : {10.0, 10.0, 1.0, 1.0})
    {
        profiler.record(section, ms);
        profiler.endFrame();
    }
    ASSERT_TRUE(profiler.takeNewReport());
    EXPECT_DOUBLE_EQ(profiler.report()[0].average_ms, 1.0); // the second window only
}
