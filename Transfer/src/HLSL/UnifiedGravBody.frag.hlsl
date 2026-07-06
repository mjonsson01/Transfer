struct VertexOutput
{
    float4 clipPos : SV_POSITION;
    float4 color   : COLOR0;
    float2 localPos : TEXCOORD0;

    float logMass   : TEXCOORD1;
    float temperature : TEXCOORD10;
    float charge : TEXCOORD11;

    uint flags : TEXCOORD12;
    uint seed : TEXCOORD13;
};

float4 main(VertexOutput input) : SV_Target
{
    float dist = length(input.localPos);

    if (dist > 1.0)
    {
        discard;
    }

    const uint viewMode = 1; // TBI: Mass View

    float4 color;

    if (viewMode == 0)
    {
        // "Real" view
        color = input.color;
    }
    else if (viewMode == 1)
    {
        //--------------------------------
        // Mass View
        //--------------------------------

        float massMag = saturate(abs(input.logMass) / 16.0);

        float3 c;

        if (input.logMass < 0.0)
        {
            // Negative mass:
            // white -> red
            c = lerp(
                float3(1.0, 1.0, 1.0),
                float3(1.0, 0.0, 0.0),
                massMag);
        }
        else
        {
            // Positive mass:
            // white -> blue
            c = lerp(
                float3(1.0, 1.0, 1.0),
                float3(0.0, 0.0, 1.0),
                massMag);
        }

        color = float4(c, 1.0);
    }
    else
    {
        color = float4(1, 0, 1, 1); // debug magenta
    }

    //------------------------------------
    // Opacity rules
    //------------------------------------

    bool isMacroGhost = (input.flags & (1u << 3)) != 0;
    bool isCollidable = (input.flags & (1u << 5)) != 0;
    bool isPreview = (input.flags & (1u<<15)) != 0;
    if (isMacroGhost || !isCollidable || isPreview)
    {
        color.a *= 0.7;
    }

    return color;
}