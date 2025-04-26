#version 330 core
#define M_PI 3.1415926535897932384626433832795

layout(location = 0) in vec3 pos;
layout(location = 1) in vec3 dir;
layout(location = 2) in float radius;
layout(location = 3) in float length;
layout(location = 4) in vec3 color;
layout(location = 5) in int activeFlag;

out CellData {
    vec4 pos;
    mat3 rot;
    float outlineRadius;
    float radius;
    float length;
    float dist;
    vec3 color;
    int activeFlag;
} cell;

uniform bool outline;

uniform mat4 view;
uniform mat4 proj;

void main() {
    // i.e. rodrigues' rotation formula
    if (length(cross(dir, vec3(1, 0, 0))) == 0) {
        cell.rot = mat3(1.0);
    } else {
        vec3 rotAxis = cross(vec3(-1, 0, 0), dir);
        rotAxis = normalize(rotAxis);
        float rotAngle = acos(dir.x / length(dir));
        mat3 K = mat3(
                0, -rotAxis.z, rotAxis.y,
                rotAxis.z, 0, -rotAxis.x,
                -rotAxis.y, rotAxis.x, 0);
        cell.rot = mat3(1.0) + sin(rotAngle) * K + (1 - cos(rotAngle)) * K * K;
    }

    cell.pos = vec4(pos, 1.0);
    cell.dist = length(proj * view * cell.pos);
    cell.outlineRadius = radius;
    cell.activeFlag = activeFlag;
    cell.length = length;

    if (outline) {
        cell.radius = radius;
        cell.color = vec3(1.0, 1.0, 1.0);
    } else {
        cell.radius = radius * 0.8;
        cell.color = color;
    }
}
