#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

extern void cocoa_start(int width, int height, void (*callback)(void));
extern void draw_pixel(int x, int y, uint8_t r, uint8_t g, uint8_t b, uint8_t a);
extern void present_frame(void);

#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 800

typedef struct
{
    float x, y, z;
} Vector3;

typedef struct
{
    Vector3 points[3];
} Triangle;

typedef struct
{
    Triangle *triangles;
    size_t triangleCount;
} Mesh;

typedef struct
{
    float m[4][4];
} Matrix4x4;

Mesh cubeMesh;
Matrix4x4 projectionMatrix;
Vector3 Camera;
float theta = 0.0f;

void MultiplyMatrixVector(const Vector3 *in, Vector3 *out, const Matrix4x4 *m)
{
    float x = in->x, y = in->y, z = in->z;
    out->x = x * m->m[0][0] + y * m->m[1][0] + z * m->m[2][0] + m->m[3][0];
    out->y = x * m->m[0][1] + y * m->m[1][1] + z * m->m[2][1] + m->m[3][1];
    out->z = x * m->m[0][2] + y * m->m[1][2] + z * m->m[2][2] + m->m[3][2];
    float w = x * m->m[0][3] + y * m->m[1][3] + z * m->m[2][3] + m->m[3][3];

    if (w != 0.0f)
    {
        out->x /= w;
        out->y /= w;
        out->z /= w;
    }
}

void DrawLine(int x0, int y0, int x1, int y1)
{
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy, e2;

    while (1)
    {
        if (x0 >= 0 && x0 < SCREEN_WIDTH && y0 >= 0 && y0 < SCREEN_HEIGHT)
            draw_pixel(x0, y0, 255, 255, 255, 255);

        if (x0 == x1 && y0 == y1)
            break;

        e2 = 2 * err;
        if (e2 >= dy)
        {
            err += dy;
            x0 += sx;
        }
        if (e2 <= dx)
        {
            err += dx;
            y0 += sy;
        }
    }
}

void DrawTriangle(int x1, int y1, int x2, int y2, int x3, int y3)
{
    DrawLine(x1, y1, x2, y2);
    DrawLine(x2, y2, x3, y3);
    DrawLine(x3, y3, x1, y1);
}

void FillTriangle(int x1, int y1, int x2, int y2, int x3, int y3, uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
    // Sort the points by y-coordinate ascending (y1 <= y2 <= y3)
    if (y1 > y2)
    {
        int t = y1;
        y1 = y2;
        y2 = t;
        t = x1;
        x1 = x2;
        x2 = t;
    }
    if (y1 > y3)
    {
        int t = y1;
        y1 = y3;
        y3 = t;
        t = x1;
        x1 = x3;
        x3 = t;
    }
    if (y2 > y3)
    {
        int t = y2;
        y2 = y3;
        y3 = t;
        t = x2;
        x2 = x3;
        x3 = t;
    }

    int total_height = y3 - y1;
    for (int i = 0; i < total_height; i++)
    {
        int second_half = i > y2 - y1 || y2 == y1;
        int segment_height = second_half ? y3 - y2 : y2 - y1;
        float alpha = total_height == 0 ? 0.0f : (float)i / total_height;
        float beta = segment_height == 0 ? 0.0f : (float)(i - (second_half ? y2 - y1 : 0)) / segment_height;

        int ax = x1 + (int)((x3 - x1) * alpha);
        int bx = second_half
                     ? x2 + (int)((x3 - x2) * beta)
                     : x1 + (int)((x2 - x1) * beta);
        int ay = y1 + i;
        int by = ay;
        if (ax > bx)
        {
            int t = ax;
            ax = bx;
            bx = t;
        }
        for (int j = ax; j <= bx; j++)
        {
            if (j >= 0 && j < SCREEN_WIDTH && ay >= 0 && ay < SCREEN_HEIGHT)
                draw_pixel(j, ay, r, g, b, a);
        }
    }
}

void ClearScreen()
{
    for (int y = 0; y < SCREEN_HEIGHT; y++)
        for (int x = 0; x < SCREEN_WIDTH; x++)
            draw_pixel(x, y, 0, 0, 0, 0); // black
}

void InitEngine()
{
    cubeMesh.triangleCount = 12;
    cubeMesh.triangles = malloc(sizeof(Triangle) * cubeMesh.triangleCount);

    Vector3 vertices[] = {
        {-0.5f, -0.5f, -0.5f}, {-0.5f, 0.5f, -0.5f}, {0.5f, 0.5f, -0.5f}, {0.5f, -0.5f, -0.5f}, {0.5f, -0.5f, 0.5f}, {0.5f, 0.5f, 0.5f}, {-0.5f, 0.5f, 0.5f}, {-0.5f, -0.5f, 0.5f}};

    Triangle tris[] = {
        {vertices[0], vertices[1], vertices[2]}, {vertices[0], vertices[2], vertices[3]}, // SOUTH
        {vertices[3], vertices[2], vertices[5]},
        {vertices[3], vertices[5], vertices[4]}, // EAST
        {vertices[4], vertices[5], vertices[6]},
        {vertices[4], vertices[6], vertices[7]}, // NORTH
        {vertices[7], vertices[6], vertices[1]},
        {vertices[7], vertices[1], vertices[0]}, // WEST
        {vertices[1], vertices[6], vertices[5]},
        {vertices[1], vertices[5], vertices[2]}, // TOP
        {vertices[4], vertices[7], vertices[0]},
        {vertices[4], vertices[0], vertices[3]} // BOTTOM
    };

    memcpy(cubeMesh.triangles, tris, sizeof(tris));

    float nearPlane = 0.1f;
    float farPlane = 1000.0f;
    float fov = 90.0f;
    float aspectRatio = (float)SCREEN_HEIGHT / (float)SCREEN_WIDTH;
    float fovRad = 1.0f / tanf(fov * 0.5f * (3.14159f / 180.0f));

    memset(&projectionMatrix, 0, sizeof(projectionMatrix));
    projectionMatrix.m[0][0] = aspectRatio * fovRad;
    projectionMatrix.m[1][1] = fovRad;
    projectionMatrix.m[2][2] = farPlane / (farPlane - nearPlane);
    projectionMatrix.m[3][2] = (-farPlane * nearPlane) / (farPlane - nearPlane);
    projectionMatrix.m[2][3] = 1.0f;
    projectionMatrix.m[3][3] = 0.0f;
}

void Tick()
{
    static int frameCount = 0;
    printf("Frame: %d\n", frameCount++);
    ClearScreen();

    theta += 0.05f;

    Matrix4x4 rotationZ = {0};
    Matrix4x4 rotationX = {0};

    rotationZ.m[0][0] = cosf(theta);
    rotationZ.m[0][1] = sinf(theta);
    rotationZ.m[1][0] = -sinf(theta);
    rotationZ.m[1][1] = cosf(theta);
    rotationZ.m[2][2] = 1.0f;
    rotationZ.m[3][3] = 1.0f;

    rotationX.m[0][0] = 1.0f;
    rotationX.m[1][1] = cosf(theta * 0.5f);
    rotationX.m[1][2] = sinf(theta * 0.5f);
    rotationX.m[2][1] = -sinf(theta * 0.5f);
    rotationX.m[2][2] = cosf(theta * 0.5f);
    rotationX.m[3][3] = 1.0f;

    for (size_t i = 0; i < cubeMesh.triangleCount; i++)
    {
        Triangle projected, translated, rotatedZ, rotatedZX;

        for (int j = 0; j < 3; j++)
            MultiplyMatrixVector(&cubeMesh.triangles[i].points[j], &rotatedZ.points[j], &rotationZ);

        for (int j = 0; j < 3; j++)
            MultiplyMatrixVector(&rotatedZ.points[j], &rotatedZX.points[j], &rotationX);

        translated = rotatedZX;
        for (int j = 0; j < 3; j++)
            translated.points[j].z += 3.0f;

        Vector3 normal, line1, line2;
        line1.x = translated.points[1].x - translated.points[0].x;
        line1.y = translated.points[1].y - translated.points[0].y;
        line1.z = translated.points[1].z - translated.points[0].z;

        line2.x = translated.points[2].x - translated.points[0].x;
        line2.y = translated.points[2].y - translated.points[0].y;
        line2.z = translated.points[2].z - translated.points[0].z;

        normal.x = line1.y * line2.z - line1.z * line2.y;
        normal.y = line1.z * line2.x - line1.x * line2.z;
        normal.z = line1.x * line2.y - line1.y * line2.x;

        float length = sqrtf(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
        normal.x /= length;
        normal.y /= length;
        normal.z /= length;

        if (normal.x * (translated.points[0].x - Camera.x) +
                normal.y * (translated.points[0].y - Camera.y) +
                normal.z * (translated.points[0].z - Camera.z) >
            0.0f)
            continue; // Culling based on camera position

        Vector3 lightDirection = {0.0f, 0.0f, -1.0f};
        float lightLength = sqrtf(lightDirection.x * lightDirection.x +
                                  lightDirection.y * lightDirection.y +
                                  lightDirection.z * lightDirection.z);
        lightDirection.x /= lightLength;
        lightDirection.y /= lightLength;
        lightDirection.z /= lightLength;

        float dotProduct = normal.x * lightDirection.x + normal.y * lightDirection.y + normal.z * lightDirection.z;

        for (int j = 0; j < 3; j++)
            MultiplyMatrixVector(&translated.points[j], &projected.points[j], &projectionMatrix);

        for (int j = 0; j < 3; j++)
        {
            projected.points[j].x += 1.0f;
            projected.points[j].y += 1.0f;
            projected.points[j].x *= 0.5f * SCREEN_WIDTH;
            projected.points[j].y *= 0.5f * SCREEN_HEIGHT;
        }

        // Clamp dotProduct to [0, 1]
        float shade = dotProduct;
        if (shade < 0)
            shade = 0;
        if (shade > 1)
            shade = 1;

        // Calculate shaded color (white base)
        uint8_t r = (uint8_t)(255 * shade);
        uint8_t g = (uint8_t)(255 * shade);
        uint8_t b = (uint8_t)(255 * shade);

        FillTriangle(
            (int)projected.points[0].x, (int)projected.points[0].y,
            (int)projected.points[1].x, (int)projected.points[1].y,
            (int)projected.points[2].x, (int)projected.points[2].y, r, g, b, 255);
    }

    present_frame();
}

int main()
{
    InitEngine();
    cocoa_start(SCREEN_WIDTH, SCREEN_HEIGHT, Tick);
    return 0;
}