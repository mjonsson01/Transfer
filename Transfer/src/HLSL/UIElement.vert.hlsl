cbuffer ScreenConstants : register(b0, space1)
{
    float2 screenSize;
};

struct VertexOutput
{
    float4 clipPos : SV_POSITION;
    float2 uv : TEXCOORD0;
    float4 color : TEXCOORD1;
    nointerpolation uint mode : TEXCOORD2;
};

VertexOutput main(
    float2 pos   : POSITION0,
    float2 uv    : TEXCOORD0,
    float4 color : TEXCOORD1,
    uint zIndex  : TEXCOORD2,
    uint mode    : TEXCOORD3)
{
    VertexOutput output;
    // float2 ndc = (pos / float2(1280.0, 720.0)) * 2.0 - 1.0;
    // ndc.y *= -1.0;
    float2 ndc = (pos / screenSize) * 2.0 - 1.0;
    ndc.y *= -1.0;
    output.clipPos = float4(ndc, 0.0, 1.0);
    output.uv = uv;
    output.color = color;
    output.mode = mode;
    return output;
}