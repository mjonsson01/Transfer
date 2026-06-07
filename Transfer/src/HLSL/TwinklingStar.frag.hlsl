struct VertexOutput
{
    float4 clipPos : SV_POSITION;
    float2 localPos : TEXCOORD0;

    float alpha : TEXCOORD1;
    float twinkleSpeed : TEXCOORD2;
    uint seed : TEXCOORD3;
};

float4 main(VertexOutput input) : SV_Target
{
    float dist = length(input.localPos);

    if (dist > 1.0)
    {
        discard;
    }

    return float4(
        1.0,
        1.0,
        1.0,
        input.alpha);
}