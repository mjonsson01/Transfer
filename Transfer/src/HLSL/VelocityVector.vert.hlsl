cbuffer CameraConstants : register(b0, space1)
{
    float screenWidth;
    float screenHeight;
    float zoom;
    float offsetX;
    float offsetY;
    uint viewMode;
    float _padding0;
    float _padding1;
};

struct VertexOutput
{
    float4 clipPos : SV_POSITION;
    float4 color : COLOR0;
};

VertexOutput main(
    float2 pos   : POSITION0,
    float4 color : TEXCOORD0)
{
    VertexOutput output;

    float2 screenPos = (pos + float2(offsetX, offsetY)) * zoom;
    float2 normalizedPos = (screenPos / float2(screenWidth, screenHeight)) * 2.0 - 1.0;
    normalizedPos.y *= -1.0;

    output.clipPos = float4(normalizedPos, 0.0, 1.0);
    output.color = color;
    return output;
}