#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define WIDTH 400
#define HEIGHT 400
#define ARRAY_LEN(arr) (sizeof(arr) / sizeof((arr)[0]))

extern void cocoa_start(int width, int height, void (*callback)(void));
extern void draw_pixel(int x, int y, uint8_t r, uint8_t g, uint8_t b, uint8_t a);
extern void present_frame(void);

typedef struct
{
    int x, y;
} Vector2D;

typedef struct
{
    float x, y, z;
} Vector3D;

typedef struct
{
    Vector3D p[3];
} Triangle;

typedef struct
{
    Triangle *tris;
} Mesh;

#define V3(x, y, z) ((Vector3D){x, y, z})
#define TRI(a, b, c) ((Triangle){.p = {a, b, c}})

typedef struct
{
    float m[4][4];
} Matrix4x4;

void make_rotation_z(Matrix4x4 *m, float theta);
void make_rotation_x(Matrix4x4 *m, float theta);

void MultiplyMatrixVector(const Vector3D *i, Vector3D *o, const Matrix4x4 *m)
{
    o->x = i->x * m->m[0][0] + i->y * m->m[1][0] + i->z * m->m[2][0] + m->m[3][0];
    o->y = i->x * m->m[0][1] + i->y * m->m[1][1] + i->z * m->m[2][1] + m->m[3][1];
    o->z = i->x * m->m[0][2] + i->y * m->m[1][2] + i->z * m->m[2][2] + m->m[3][2];
    float w = i->x * m->m[0][3] + i->y * m->m[1][3] + i->z * m->m[2][3] + m->m[3][3];

    if (w != 0.0f)
    {
        o->x /= w;
        o->y /= w;
        o->z /= w;
    }
}

void draw_line(int x0, int y0, int x1, int y1, uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
    int dx = abs(x1 - x0);
    int dy = abs(y1 - y0);
    int sx = (x0 < x1) ? 1 : -1;
    int sy = (y0 < y1) ? 1 : -1;
    int err = dx - dy;

    while (1)
    {
        draw_pixel(x0, y0, r, g, b, a);
        if (x0 == x1 && y0 == y1)
            break;
        int e2 = err * 2;
        if (e2 > -dy)
        {
            err -= dy;
            x0 += sx;
        }
        if (e2 < dx)
        {
            err += dx;
            y0 += sy;
        }
    }
}

void draw_trangle(int x0, int y0, int x1, int y1, int x2, int y2, uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
    draw_line(x0, y0, x1, y1, r, g, b, a);
    draw_line(x1, y1, x2, y2, r, g, b, a);
    draw_line(x2, y2, x0, y0, r, g, b, a);
}

int frame = 0;

void update_frame()
{
    // for (int y = 0; y < HEIGHT; y++)
    // {
    //     for (int x = 0; x < WIDTH; x++)
    //     {
    //         draw_pixel(x, y, (x + frame) % 255, (y + frame) % 255, 128, 255);
    //     }
    // }
    // present_frame();
    // frame++;

    for (int y = 0; y < HEIGHT; y++)
    {
        for (int x = 0; x < WIDTH; x++)
        {
            draw_pixel(x, y, 0, 0, 0, 155); // Clear the screen to black
        }
    }

    static Triangle cube_tris[] = {
        // Front face (z = 0)
        TRI(V3(0, 0, 0), V3(1, 0, 0), V3(1, 1, 0)),
        TRI(V3(0, 0, 0), V3(1, 1, 0), V3(0, 1, 0)),

        // Back face (z = 1)
        TRI(V3(0, 0, 1), V3(1, 1, 1), V3(1, 0, 1)),
        TRI(V3(0, 0, 1), V3(0, 1, 1), V3(1, 1, 1)),

        // Left face (x = 0)
        TRI(V3(0, 0, 0), V3(0, 1, 1), V3(0, 0, 1)),
        TRI(V3(0, 0, 0), V3(0, 1, 0), V3(0, 1, 1)),

        // Right face (x = 1)
        TRI(V3(1, 0, 0), V3(1, 0, 1), V3(1, 1, 1)),
        TRI(V3(1, 0, 0), V3(1, 1, 1), V3(1, 1, 0)),

        // Top face (y = 1)
        TRI(V3(0, 1, 0), V3(1, 1, 1), V3(1, 1, 0)),
        TRI(V3(0, 1, 0), V3(0, 1, 1), V3(1, 1, 1)),

        // Bottom face (y = 0)
        TRI(V3(0, 0, 0), V3(1, 0, 1), V3(1, 0, 0)),
        TRI(V3(0, 0, 0), V3(0, 0, 1), V3(1, 0, 1)),
    };

    static Mesh cube = {
        .tris = cube_tris};

    float fNear = 0.1f;
    float fFar = 1000.0f;
    float fFov = 145.0f;
    float fAspectRation = (float)WIDTH / (float)HEIGHT;
    float fFovRad = 1.0f / tanf(fFov * 0.5f * (3.14159265358979323846f / 180.0f));

    Matrix4x4 projection = {0};
    projection.m[0][0] = fAspectRation * fFovRad;
    projection.m[1][1] = fFovRad;
    projection.m[2][2] = fFar / (fFar - fNear);
    projection.m[3][2] = (-fFar * fNear) / (fFar - fNear);
    projection.m[2][3] = 1.0f;
    projection.m[3][3] = 0.0f;

    float theta = frame * 0.02f; // Adjust speed as desired

    Matrix4x4 rotZ, rotX;
    make_rotation_z(&rotZ, theta);
    make_rotation_x(&rotX, theta * 0.5f);

    for (int i = 0; i < ARRAY_LEN(cube_tris); i++)
    {
        Triangle *tri = &cube.tris[i];
        Triangle rotated, projected;

        // Rotate in Z
        MultiplyMatrixVector(&tri->p[0], &rotated.p[0], &rotZ);
        MultiplyMatrixVector(&tri->p[1], &rotated.p[1], &rotZ);
        MultiplyMatrixVector(&tri->p[2], &rotated.p[2], &rotZ);

        // Rotate in X
        MultiplyMatrixVector(&rotated.p[0], &rotated.p[0], &rotX);
        MultiplyMatrixVector(&rotated.p[1], &rotated.p[1], &rotX);
        MultiplyMatrixVector(&rotated.p[2], &rotated.p[2], &rotX);

        // Offset into the screen (so the cube is in front of the camera)
        rotated.p[0].z += 2.0f;
        rotated.p[1].z += 2.0f;
        rotated.p[2].z += 2.0f;

        // Project
        MultiplyMatrixVector(&rotated.p[0], &projected.p[0], &projection);
        MultiplyMatrixVector(&rotated.p[1], &projected.p[1], &projection);
        MultiplyMatrixVector(&rotated.p[2], &projected.p[2], &projection);

        // Map to screen
        projected.p[0].x = (projected.p[0].x + 1.0f) * 0.5f * WIDTH;
        projected.p[0].y = (projected.p[0].y + 1.0f) * 0.5f * HEIGHT;
        projected.p[1].x = (projected.p[1].x + 1.0f) * 0.5f * WIDTH;
        projected.p[1].y = (projected.p[1].y + 1.0f) * 0.5f * HEIGHT;
        projected.p[2].x = (projected.p[2].x + 1.0f) * 0.5f * WIDTH;
        projected.p[2].y = (projected.p[2].y + 1.0f) * 0.5f * HEIGHT;

        draw_trangle(
            projected.p[0].x, projected.p[0].y,
            projected.p[1].x, projected.p[1].y,
            projected.p[2].x, projected.p[2].y,
            255, 255, 255, 255);
    }

    present_frame();
    frame++;
}

void make_rotation_z(Matrix4x4 *m, float theta)
{
    m->m[0][0] = cosf(theta);
    m->m[0][1] = sinf(theta);
    m->m[0][2] = 0.0f;
    m->m[0][3] = 0.0f;
    m->m[1][0] = -sinf(theta);
    m->m[1][1] = cosf(theta);
    m->m[1][2] = 0.0f;
    m->m[1][3] = 0.0f;
    m->m[2][0] = 0.0f;
    m->m[2][1] = 0.0f;
    m->m[2][2] = 1.0f;
    m->m[2][3] = 0.0f;
    m->m[3][0] = 0.0f;
    m->m[3][1] = 0.0f;
    m->m[3][2] = 0.0f;
    m->m[3][3] = 1.0f;
}

void make_rotation_x(Matrix4x4 *m, float theta)
{
    m->m[0][0] = 1.0f;
    m->m[0][1] = 0.0f;
    m->m[0][2] = 0.0f;
    m->m[0][3] = 0.0f;
    m->m[1][0] = 0.0f;
    m->m[1][1] = cosf(theta);
    m->m[1][2] = sinf(theta);
    m->m[1][3] = 0.0f;
    m->m[2][0] = 0.0f;
    m->m[2][1] = -sinf(theta);
    m->m[2][2] = cosf(theta);
    m->m[2][3] = 0.0f;
    m->m[3][0] = 0.0f;
    m->m[3][1] = 0.0f;
    m->m[3][2] = 0.0f;
    m->m[3][3] = 1.0f;
}

int main()
{
    cocoa_start(WIDTH, HEIGHT, update_frame);
    return 0;
}
