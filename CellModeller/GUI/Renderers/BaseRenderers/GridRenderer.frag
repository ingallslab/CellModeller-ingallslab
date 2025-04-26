#version 330 core
in vec3 nearPoint;
in vec3 farPoint;

uniform int minTick;
uniform int majTick;
uniform vec3 linecolor;

uniform mat4 view;
uniform mat4 proj;

vec4 grid(vec3 pos, float scale) {
    vec2 coord = pos.xy * scale;

    // T_T what is this black magic...
    vec2 derivative = fwidth(coord);
    vec2 grid = abs(fract(coord - 0.5) - 0.5) / derivative;
    float line = min(grid.x, grid.y);
    vec4 color = vec4(linecolor, 1.0 - min(line, 1.0));

    float minX = 5 * min(derivative.x, 1.0);
    bool xAxis = pos.x > -minX && pos.x < minX;

    float minY = 5 * min(derivative.y, 1.0);
    bool yAxis = pos.y > -minY && pos.y < minY;

    if (xAxis && yAxis)
        color.rgb = vec3(0.5, 0.5, 0.0);
    else if (yAxis)
        color.rgb = vec3(0.0, 0.5, 0.0);
    else if (xAxis)
        color.rgb = vec3(0.5, 0.0, 0.0);
    return color;
}

float computeDepth(vec3 pos) {
    vec4 clipSpacePos = proj * view * vec4(pos, 1.0);
    float depth = clipSpacePos.z / clipSpacePos.w;
    return (depth + 1) / 2;
}

float computeLinearDepth(vec3 pos) {
    vec4 clipSpacePos = proj * view * vec4(pos, 1.0);
    float dist = distance(clipSpacePos.xyz, nearPoint);
    float maxDist = distance(nearPoint, farPoint);
    return dist / maxDist;
}

void main() {
    float t = -nearPoint.z / (farPoint.z - nearPoint.z);
    vec3 pos = nearPoint + t * (farPoint - nearPoint);

    float linearDepth = computeLinearDepth(pos);
    float fading = max(0, 0.5 - linearDepth);

    gl_FragDepth = computeDepth(pos);
    gl_FragColor = grid(pos, 1 / float(minTick)) * float(t > 0);
    gl_FragColor += grid(pos, 1 / float(majTick)) * float(t > 0);
    gl_FragColor.a *= fading;
}
