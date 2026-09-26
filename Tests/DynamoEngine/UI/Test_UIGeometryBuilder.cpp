// File: Tests/DynamoEngine/UI/Test_UIGeometryBuilder.cpp

// Test Framework Imports
#include <gtest/gtest.h>

// Custom Imports
#include "DynamoEngine/UI/UIGeometryBuilder.hpp"

// Standard Library Imports
#include <cstdint>
#include <vector>

using namespace DynamoEngine;

class UIGeometryBuilderTest : public ::testing::Test
{
  protected:
    std::vector<UIVertex> vertices;
    FontAtlas empty_atlas; // no font loaded: every glyph has zero size and zero advance
    UIGeometryBuilder builder{vertices, empty_atlas};
};

TEST_F(UIGeometryBuilderTest, RectIsTwoTrianglesCoveringItsCorners)
{
    builder.addRect(SDL_FRect{10.0f, 20.0f, 30.0f, 40.0f}, SDL_Color{255, 0, 0, 255});

    ASSERT_EQ(vertices.size(), 6u);
    EXPECT_FLOAT_EQ(vertices[0].x, 10.0f); // top-left
    EXPECT_FLOAT_EQ(vertices[0].y, 20.0f);
    EXPECT_FLOAT_EQ(vertices[1].x, 40.0f); // top-right
    EXPECT_FLOAT_EQ(vertices[1].y, 20.0f);
    EXPECT_FLOAT_EQ(vertices[2].x, 10.0f); // bottom-left
    EXPECT_FLOAT_EQ(vertices[2].y, 60.0f);
    EXPECT_FLOAT_EQ(vertices[4].x, 40.0f); // bottom-right
    EXPECT_FLOAT_EQ(vertices[4].y, 60.0f);
}

TEST_F(UIGeometryBuilderTest, RectIsSolidWithNormalizedColor)
{
    builder.addRect(SDL_FRect{0.0f, 0.0f, 1.0f, 1.0f}, SDL_Color{255, 0, 51, 255});

    for (const UIVertex& vertex : vertices)
    {
        EXPECT_EQ(vertex.mode, static_cast<uint32_t>(UIVertexMode::Solid));
        EXPECT_FLOAT_EQ(vertex.r, 1.0f);
        EXPECT_FLOAT_EQ(vertex.g, 0.0f);
        EXPECT_FLOAT_EQ(vertex.b, 0.2f); // 51 / 255
        EXPECT_FLOAT_EQ(vertex.a, 1.0f);
    }
}

TEST_F(UIGeometryBuilderTest, TextIsOneTexturedQuadPerCharacter)
{
    builder.addText("abc", Vector2F(0.0f, 0.0f));

    ASSERT_EQ(vertices.size(), 3u * 6u);
    for (const UIVertex& vertex : vertices)
    {
        EXPECT_EQ(vertex.mode, static_cast<uint32_t>(UIVertexMode::Textured));
    }
}

TEST_F(UIGeometryBuilderTest, BuilderAppendsWithoutClearing)
{
    vertices.resize(2); // vertices from an element drawn earlier this frame
    builder.addRect(SDL_FRect{0.0f, 0.0f, 1.0f, 1.0f}, SDL_Color{255, 255, 255, 255});
    EXPECT_EQ(vertices.size(), 2u + 6u);
}

TEST_F(UIGeometryBuilderTest, CenteredTextStartsAtAreaCenterWhenTextHasNoSize)
{
    // Empty atlas: text width and font height are 0, so "centered" lands exactly on the center point
    builder.addTextCentered("x", SDL_FRect{100.0f, 50.0f, 200.0f, 80.0f});

    ASSERT_EQ(vertices.size(), 6u);
    EXPECT_FLOAT_EQ(vertices[0].x, 200.0f); // 100 + 200/2
    EXPECT_FLOAT_EQ(vertices[0].y, 90.0f);  // 50 + 80/2
}