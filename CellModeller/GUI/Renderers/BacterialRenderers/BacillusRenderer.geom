#version 330 core
#define M_PI 3.1415926535897932384626433832795

layout(points) in;
layout(triangle_strip, max_vertices = 224) out;
// vertex count = 4N(N+1)
//  3   4   5   6   7
// 48  80 120 168 224

in CellData {
    vec4 pos;
    mat3 rot;
    float outlineRadius;
    float radius;
    float length;
    float dist;
    vec3 color;
    int activeFlag;
} cell[];

out vec3 fColor;

uniform mat4 view;
uniform mat4 proj;

vec3 projSphere(float phi, float theta) {
    return vec3(
        cos(phi),
        sin(phi) * cos(theta),
        sin(phi) * sin(theta)
    );
}

void hemisphere(int N, vec4 pos, float offset, float radius, bool top) {
    int start, end;
    float orientation;
    if (top) {
        start = 1;
        end = int(floor(N / 2));
        orientation = 1;
    } else {
        if (N % 2 == 0)
            start = int(floor(N / 2) + 1);
        else
            start = int(floor(N / 2) + 2);
        end = N;
        orientation = -1;
    }
    int sectorStart = 0, sectorEnd = 2 * N, stepDir = 1;

    float phi_1 = M_PI * float(start - 1) / float(N);
    for (int i = start; i <= end; i++) {
        float phi_2 = M_PI * float(i) / float(N);

        for (int j = sectorStart; abs(j - sectorEnd) != 0; j += stepDir) {
            float theta = 2 * M_PI * float(j) / float(2 * N);
            gl_Position = proj * view * (pos
                        + vec4(cell[0].rot * (
                                vec3(orientation * offset, 0.0, 0.0)
                                    + radius * projSphere(phi_1, theta)),
                            0.0));
            EmitVertex();
            gl_Position = proj * view * (pos
                        + vec4(cell[0].rot * (
                                vec3(orientation * offset, 0.0, 0.0)
                                    + radius * projSphere(phi_2, theta)),
                            0.0));
            EmitVertex();
        }

        phi_1 = phi_2;
    }
}

void cylinderpipe(int N, vec4 pos, float offset, float radius) {
    float phi_1 = M_PI * floor(N / 2) / float(N);
    float phi_2;
    if (N % 2 == 0)
        phi_2 = M_PI * floor(N / 2) / float(N);
    else
        phi_2 = M_PI * (floor(N / 2) + 1) / float(N);
    int sectorStart = 0, sectorEnd = 2 * N, stepDir = 1;

    for (int j = sectorStart; abs(j - sectorEnd) != 0; j += stepDir) {
        float theta = 2 * M_PI * float(j) / float(2 * N);
        gl_Position = proj * view * (pos
                    + vec4(cell[0].rot * (
                            vec3(offset, 0.0, 0.0)
                                + radius * projSphere(phi_1, theta)),
                        0.0));
        EmitVertex();
        gl_Position = proj * view * (pos
                    + vec4(cell[0].rot * (
                            vec3(-offset, 0.0, 0.0)
                                + radius * projSphere(phi_2, theta)),
                        0.0));
        EmitVertex();
    }
}

void main() {
    if (cell[0].activeFlag == 0)
        return;
    fColor = cell[0].color;

    // there is a reason for this scaling i promise
    // (no there isn't)
    int lod = 3;
    if (cell[0].dist < 48 * cell[0].outlineRadius) {
        lod = 7;
    } else if (cell[0].dist < 80 * cell[0].outlineRadius) {
        lod = 6;
    } else if (cell[0].dist < 120 * cell[0].outlineRadius) {
        lod = 5;
    } else if (cell[0].dist < 168 * cell[0].outlineRadius) {
        lod = 4;
    }

    hemisphere(lod, cell[0].pos, cell[0].length / 2, cell[0].radius, true);
    cylinderpipe(lod, cell[0].pos, cell[0].length / 2, cell[0].radius);
    hemisphere(lod, cell[0].pos, cell[0].length / 2, cell[0].radius, false);
    EndPrimitive();
}
