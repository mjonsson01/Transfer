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
};

// The fade-in, computed here so a segment's vertices never change once uploaded (RenderSystem only uploads NEW
// segments each frame). Must match RenderSystem's TrajectoryConstants, field for field.
cbuffer TrajectoryConstants : register(b1, space1)
{
    float clock;          // Trajectories::clock: real seconds since the trace started
    float fadeInSeconds;  // TRAJECTORY_FADE_IN_SECONDS: how long a new piece takes to settle
    float fadeStartAlpha; // TRAJECTORY_FADE_START_ALPHA: how visible a brand-new piece is
    float _padding2;
};

struct VertexOutput
{
    float4 clipPos : SV_POSITION;
    float4 color : COLOR0;
};

// One end of a trajectory line segment (RenderSystem's TrajectoryLineVertex)
VertexOutput main(
    float2 pos     : POSITION0,
    float2 prevPos : TEXCOORD0,
    float4 color   : TEXCOORD1, // .a = the settled alpha
    float2 fade    : TEXCOORD2) // x = when this point was drawn (on the clock), y = 1 if it fades in, 0 if not
{
    VertexOutput output;

    // Same interpolation and world -> screen steps as the bodies, so the first point sits on the body as drawn
    float2 interpPos = lerp(prevPos, pos, rendering_alpha);
    float2 screenPos = (interpPos + float2(offsetX, offsetY)) * zoom;
    float2 normalizedPos = (screenPos / float2(screenWidth, screenHeight)) * 2.0 - 1.0;
    normalizedPos.y *= -1.0;

    float alpha = color.a;
    if (fade.y > 0.5)
    {
        float progress = saturate((clock - fade.x) / fadeInSeconds);
        alpha = lerp(fadeStartAlpha, color.a, progress);
    }

    output.clipPos = float4(normalizedPos, 0.0, 1.0);
    output.color = float4(color.rgb, alpha);
    return output;
}
