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

    float2 finalPos = pos + local * radius;

    float2 normalizedPos =
        (finalPos / float2(1280.0, 720.0)) * 2.0 - 1.0;

    normalizedPos.y *= -1.0;

    output.clipPos = float4(normalizedPos, 0.0, 1.0);

    output.localPos = local;

    output.alpha = alpha;
    output.twinkleSpeed = twinkleSpeed;
    output.seed = seed;

    return output;
}