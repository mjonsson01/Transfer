struct VertexOutput
{
    float4 clipPos : SV_POSITION;
    float4 color : COLOR0;
    float2 uv : TEXCOORD1;
};

float4 main(VertexOutput input) : SV_Target
{
    return input.color;
}