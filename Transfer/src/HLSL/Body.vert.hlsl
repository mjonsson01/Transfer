struct VertexOutput
{
    float4 clipPos : SV_POSITION;
    float4 color   : COLOR0;

    float2 localPos : TEXCOORD8; // Not a vertex attribute from the cpu side.
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
    float2 pos         : POSITION0,
    float2 prevPos     : TEXCOORD0,
    float  radius      : TEXCOORD1,
    float  logMass     : TEXCOORD2,
    float  temperature : TEXCOORD3,
    float  charge      : TEXCOORD4,
    uint   flags       : TEXCOORD5,
    uint   seed        : TEXCOORD6,
    uint vertexId      : SV_VertexID)
{
    VertexOutput output;

    float2 local = offsets[vertexId % 6];

    float2 finalPos = pos + local * radius;

    float2 normalizedPos =
        (finalPos / float2(1280.0,720.0))*2.0 - 1.0;

    normalizedPos.y *= -1.0;

    output.clipPos = float4(normalizedPos,0,1);

    output.color = float4(1,0,0,1);

    output.localPos = local;

    return output;
}