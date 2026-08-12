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
    float2 uv       : TEXCOORD1;
};


VertexOutput main(
    float2 pos        : POSITION0,
    float2 prevPos    : TEXCOORD0,
    float square_size : TEXCOORD1,
    float2 uv         : TEXCOORD2,
    uint vertexID   : SV_vertexID
    )
{
    VertexOutput output;
    float2 interpPos = lerp(prevPos, pos, rendering_alpha);
    float2 screenPos = (interpPos + float2(offsetX, offsetY)) * zoom;
    float2 normalizedPos = (screenPos / float2(screenWidth, screenHeight)) * 2 - 1.0;
    normalizedPos.y *= -1.0;
    output.clipPos = float4(normalizedPos, 0.0, 1.0);
    output.color = float4(1.0, 1.0, 1.0, 1.0);

    return output;
}