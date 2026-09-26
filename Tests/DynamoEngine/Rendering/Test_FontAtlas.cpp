// File: Tests/DynamoEngine/Rendering/Test_FontAtlas.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/Rendering/FontAtlas.hpp"

using namespace DynamoEngine;

TEST(FontAtlas, EmptyAtlasHasNoSize)
{
    FontAtlas atlas;
    EXPECT_FLOAT_EQ(atlas.fontHeight(), 0.0f);
    EXPECT_FLOAT_EQ(atlas.measureTextWidth("hello"), 0.0f);
}

TEST(FontAtlas, CharactersOutsidePrintableAsciiGetEmptyMetrics)
{
    FontAtlas atlas;
    GlyphMetrics newline = atlas.glyph('\n'); // below the baked range
    EXPECT_FLOAT_EQ(newline.width, 0.0f);
    EXPECT_FLOAT_EQ(newline.advance_x, 0.0f);

    GlyphMetrics non_ascii =
        atlas.glyph(static_cast<char>(200)); // negative as a signed char: must not index out of bounds
    EXPECT_FLOAT_EQ(non_ascii.width, 0.0f);
}

TEST(FontAtlas, BuildAtlasWithNoFontReturnsNull)
{
    FontAtlas atlas;
    EXPECT_EQ(atlas.buildAtlas(nullptr), nullptr);
}