struct VertexOutput
{
    float4 clipPos : SV_POSITION;
    float4 color : COLOR0;
    float2 uv : TEXCOORD1;
};

// The ship sprite (alpha premultiplied, mipmapped) and how to read it; bound by RenderSystem::renderStarship
Texture2D<float4> shipTexture : register(t0, space2);
SamplerState shipSampler : register(s0, space2);

float4 main(VertexOutput input) : SV_Target
{
    return shipTexture.Sample(shipSampler, input.uv);
}