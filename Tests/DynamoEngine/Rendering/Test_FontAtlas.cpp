// File: Tests/DynamoEngine/Rendering/Test_FontAtlas.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// SDL Imports
#include <SDL3/SDL_surface.h>
#include <SDL3_ttf/SDL_ttf.h>

// Custom Imports
#include "DynamoEngine/Rendering/FontAtlas.hpp"

// Standard Library Imports
#include <cmath>

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
    EXPECT_EQ(atlas.buildAtlas(nullptr, 18.0f, 1.0f), nullptr);
}

// --------- BAKING A REAL FONT --------- //

// Opens the game's UI font (the path comes from CMakeLists.txt)
class FontAtlasBakeTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        ASSERT_TRUE(TTF_Init());
        font = TTF_OpenFont(TRANSFER_TEST_FONT_PATH, 18.0f);
        ASSERT_NE(font, nullptr);
    }
    void TearDown() override
    {
        if (font != nullptr)
        {
            TTF_CloseFont(font);
        }
        TTF_Quit();
    }

    // Bakes and immediately frees the surface: these tests only look at the metrics
    static bool bake(FontAtlas& atlas, TTF_Font* font, float pixel_scale)
    {
        SDL_Surface* surface = atlas.buildAtlas(font, 18.0f, pixel_scale);
        const bool baked = (surface != nullptr);
        SDL_DestroySurface(surface);
        return baked;
    }

    TTF_Font* font = nullptr;
};

TEST_F(FontAtlasBakeTest, MetricsStayInUIPointsAtAnyPixelScale)
{
    FontAtlas normal;
    FontAtlas retina;
    ASSERT_TRUE(bake(normal, font, 1.0f));
    ASSERT_TRUE(bake(retina, font, 2.0f)); // twice the pixels...

    // ...but the same size in UI points, so layout doesn't change (within a pixel of rounding)
    EXPECT_NEAR(retina.fontHeight(), normal.fontHeight(), 1.0f);
    EXPECT_NEAR(retina.measureTextWidth("Mass: 1.000000"), normal.measureTextWidth("Mass: 1.000000"), 1.0f);
}

TEST_F(FontAtlasBakeTest, GlyphsAreWholeScreenPixels)
{
    FontAtlas atlas;
    ASSERT_TRUE(bake(atlas, font, 1.25f));
    EXPECT_FLOAT_EQ(atlas.pixelScale(), 1.25f);

    // Converted back to screen pixels, every size is a whole number: glyphs stay on the pixel grid
    const GlyphMetrics letter = atlas.glyph('A');
    const float width_in_pixels = letter.width * 1.25f;
    const float advance_in_pixels = letter.advance_x * 1.25f;
    EXPECT_NEAR(width_in_pixels, std::round(width_in_pixels), 1e-4f);
    EXPECT_NEAR(advance_in_pixels, std::round(advance_in_pixels), 1e-4f);
}

TEST_F(FontAtlasBakeTest, EveryGlyphFitsAtLargeScales)
{
    FontAtlas atlas;
    // 72-pixel text, e.g. a 4K Retina screen with the window at full size
    ASSERT_TRUE(bake(atlas, font, 4.0f));
    EXPECT_GT(atlas.glyph('~').width, 0.0f); // '~' is baked last: if it's there, everything fit
}
