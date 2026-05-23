struct VertexOutput
{
    float4 clipPos : SV_POSITION;
    float4 color   : COLOR0;
    float2 localPos : TEXCOORD0;
};

float4 main(VertexOutput input) : SV_Target
{
    float dist = length(input.localPos);

    if (dist > 1.0)
    {
        discard;
    }

    return input.color;
}