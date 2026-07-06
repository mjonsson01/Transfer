cbuffer TwinkleConstants : register(b0, space3)
{
    float elapsedTime;
};

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

    float phase = frac(float(input.seed % 10000) * 0.61803398875) * 6.28318530718;
    float twinkle = 0.5 + 0.5 * sin(elapsedTime * input.twinkleSpeed + phase);

    // Subtle per-star color variation, from neutral white to a cool blue-white.
    float colorRand = frac(float(input.seed % 9973) * 0.41421356);
    float3 warmWhite = float3(1.0, 1.0, 1.0);
    float3 coolBlueWhite = float3(0.75, 0.85, 1.0);
    float3 starColor = lerp(warmWhite, coolBlueWhite, colorRand);
    return float4(
        starColor,
        input.alpha * twinkle);
}