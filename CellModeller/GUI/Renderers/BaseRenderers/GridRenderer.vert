#version 330 core
out vec3 nearPoint;
out vec3 farPoint;

uniform mat4 viewInv;
uniform mat4 projInv;

vec3 gridPlane[6] = vec3[](
        vec3(1, 1, 0), vec3(-1, -1, 0), vec3(-1, 1, 0),
        vec3(-1, -1, 0), vec3(1, 1, 0), vec3(1, -1, 0)
    );

vec3 UnprojectPoint(vec3 point) {
    vec4 unprojectedPoint = viewInv * projInv * vec4(point, 1.0);
    return unprojectedPoint.xyz / unprojectedPoint.w;
}

void main() {
    vec3 p = gridPlane[gl_VertexID];
    nearPoint = vec3(p.xy, 0.0);
    farPoint = vec3(p.xy, 1.0);
    nearPoint = UnprojectPoint(nearPoint);
    farPoint = UnprojectPoint(farPoint);
    gl_Position = vec4(p, 1.0);
}
