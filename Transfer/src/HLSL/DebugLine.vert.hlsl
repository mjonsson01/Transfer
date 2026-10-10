cbuffer CameraConstants : register(b0, space1)
{
    float screenWidth;
    float screenHeight;
    float zoom;
    float offsetX;
    float offsetY;
    uint viewMode;
    float rendering_alpha;
    float _padding1;
};

struct VertexOutput
{
    float4 clipPos : SV_POSITION;
    float4 color : COLOR0;
};

// One end of a debug line (RenderSystem's DebugLineVertex): world position now and one physics tick ago, and colour
VertexOutput main(
    float2 pos     : POSITION0,
    float2 prevPos : TEXCOORD0,
    float4 color   : TEXCOORD1)
{
    VertexOutput output;

    // Same interpolation and world -> screen steps as the ship sprite, so the outline sits exactly on it
    float2 interpPos = lerp(prevPos, pos, rendering_alpha);
    float2 screenPos = (interpPos + float2(offsetX, offsetY)) * zoom;
    float2 normalizedPos = (screenPos / float2(screenWidth, screenHeight)) * 2.0 - 1.0;
    normalizedPos.y *= -1.0;

    output.clipPos = float4(normalizedPos, 0.0, 1.0);
    output.color = color;
    return output;
}