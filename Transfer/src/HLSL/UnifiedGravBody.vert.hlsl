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
}
struct VertexOutput
{
    float4 clipPos : SV_POSITION;
    float4 color : COLOR0;
    float2 localPos : TEXCOORD8; // Not a vertex attribute from the cpu side.


    float logMass : TEXCOORD9;
    float temperature : TEXCOORD10;
    float charge : TEXCOORD11;

    uint flags : TEXCOORD12;
    uint seed : TEXCOORD13;
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
    uint   vertexId    : SV_VertexID)
{
    VertexOutput output;

    float2 local = offsets[vertexId % 6];

    // float2 worldPos = pos + local * radius;
    float2 interpPos = lerp(prevPos, pos, rendering_alpha);
    float2 worldPos = interpPos + local * radius;

    float2 screenPos = (worldPos + float2(offsetX, offsetY))*zoom;

    float2 normalizedPos = (screenPos / (float2(screenWidth, screenHeight)))* 2 - 1.0;
    
    normalizedPos.y *= -1.0;

    output.clipPos = float4(normalizedPos, 0.0, 1.0);

    // Required for circle clipping in fragment shader
    output.localPos = local;

    // Pass all body attributes through to fragment shader
    output.logMass = logMass;
    output.temperature = temperature;
    output.charge = charge;
    output.flags = flags;
    output.seed = seed;

    output.color = float4(1.0, 0.0, 0.0, 1.0);

    return output;
}