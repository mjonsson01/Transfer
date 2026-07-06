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
}

struct VertexOutput
{
    float4 clipPos : SV_POSITION;
    float2 localPos : TEXCOORD0;

    float alpha : TEXCOORD1;
    float twinkleSpeed : TEXCOORD2;
    uint seed : TEXCOORD3;
};

static float2 offsets[6] =
{
    float2(-1,-1),
    float2( 1,-1),
    float2(-1, 1),

    float2(-1, 1),
    float2( 1,-1),
    float2( 1, 1)
};

VertexOutput main(
    float2 pos          : POSITION0,
    float radius        : TEXCOORD0,
    float alpha         : TEXCOORD1,
    float twinkleSpeed  : TEXCOORD2,
    uint seed           : TEXCOORD3,
    uint vertexId       : SV_VertexID)
{
    VertexOutput output;

    float2 local = offsets[vertexId % 6];

    float2 worldPos = pos + local * radius;
    float2 screenPos = (worldPos + float2(offsetX, offsetY)) * zoom;

    float2 normalizedPos = (screenPos / float2(screenWidth, screenHeight)) * 2.0 - 1.0;
    normalizedPos.y *= -1.0;
    output.clipPos = float4(normalizedPos, 0.0, 1.0);

    output.localPos = local;

    output.alpha = alpha;
    output.twinkleSpeed = twinkleSpeed;
    output.seed = seed;

    return output;
}