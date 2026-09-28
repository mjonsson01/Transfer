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

    uint viewMode : TEXCOORD14;
};

// Must match VisorView in the C++ code
static const uint VIEW_REALISTIC = 0;
static const uint VIEW_MASS = 1;
static const uint VIEW_CHARGE = 2;
static const uint VIEW_TEMPERATURE = 3;

float4 main(VertexOutput input) : SV_Target
{
    float dist = length(input.localPos);

    if (dist > 1.0)
    {
        discard;
    }

    uint viewMode = input.viewMode;

    float4 color;

    if (viewMode == VIEW_REALISTIC)
    {
        // "Real" view
        color = float4(0, 1.0, 0.0, 1.0);
    }
    else if (viewMode == VIEW_MASS)
    {
        //--------------------------------
        // Mass View
        //--------------------------------

        float massMag = saturate(abs(input.logMass) / 8.0);
        // Scale brightness of a pure hue instead of blending toward white, so
        // color stays fully saturated across the whole mass range instead of
        // washing out near zero mass.
        float brightness = lerp(0.35, 1.0, massMag);
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
    else if (viewMode == VIEW_CHARGE)
    {
        color = float4(1, 0, 0, 1);
    }
    else if (viewMode == VIEW_TEMPERATURE)
    {
        color = float4(0, 0, 1, 1);
    }
    else
    {
        color = float4(1, 0, 1, 1); // debug magenta: an unknown view mode is a bug, so make it obvious
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