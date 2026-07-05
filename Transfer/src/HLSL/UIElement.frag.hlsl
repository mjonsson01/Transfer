struct VertexOutput
{
    float4 clipPos : SV_POSITION;
    float2 uv : TEXCOORD0;
    float4 color : TEXCOORD1;
    nointerpolation uint mode : TEXCOORD2;
};

Texture2D<float4> fontAtlas : register(t0, space2);
SamplerState fontSampler : register(s0, space2);

float4 main(VertexOutput input) : SV_Target
{
    if (input.mode == 1)
    {
        float alpha = fontAtlas.Sample(fontSampler, input.uv).a;
        return float4(input.color.rgb, input.color.a * alpha);
    }
    return input.color;
}
