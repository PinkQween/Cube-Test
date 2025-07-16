// main.c - 3D Software Renderer Demo
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>

// --- External Cocoa Bridge Functions ---
extern void cocoa_start(int width, int height, int maxFPS, void (*callback)(void));
extern void draw_pixel(int x, int y, uint8_t r, uint8_t g, uint8_t b, uint8_t a);
extern void present_frame(void);

#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 800
#define LINE_MAX_LEN 128

// Ambient light intensity (minimum light)
static const float AMBIENT_LIGHT_INTENSITY = 0.8f;

// --- Math Types ---
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

// --- Matrix Type ---
typedef struct
{
    float m[4][4];
} Matrix4x4;

// --- Vector3 Math Utilities ---
static inline Vector3 vec3_add(Vector3 a, Vector3 b)
{
    return (Vector3){a.x + b.x, a.y + b.y, a.z + b.z};
}
static inline Vector3 vec3_sub(Vector3 a, Vector3 b)
{
    return (Vector3){a.x - b.x, a.y - b.y, a.z - b.z};
}
static inline Vector3 vec3_mul(Vector3 a, Vector3 b)
{
    return (Vector3){a.x * b.x, a.y * b.y, a.z * b.z};
}
static inline Vector3 vec3_div(Vector3 a, Vector3 b)
{
    return (Vector3){a.x / b.x, a.y / b.y, a.z / b.z};
}
static inline Vector3 vec3_scale(Vector3 a, float s)
{
    return (Vector3){a.x * s, a.y * s, a.z * s};
}
static inline float vec3_dot(Vector3 a, Vector3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
static inline Vector3 vec3_cross(Vector3 a, Vector3 b)
{
    return (Vector3){
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x};
}
static inline float vec3_length(Vector3 a)
{
    return sqrtf(a.x * a.x + a.y * a.y + a.z * a.z);
}
static inline Vector3 vec3_normalize(Vector3 a)
{
    float len = vec3_length(a);
    if (len == 0.0f)
        return (Vector3){0, 0, 0};
    return vec3_scale(a, 1.0f / len);
}
static inline Vector3 vec3_copy(Vector3 a)
{
    return a;
}
// --- Matrix4x4 Math Utilities ---
// Set matrix to identity
static inline void mat4x4_identity(Matrix4x4 *m)
{
    memset(m, 0, sizeof(Matrix4x4));
    m->m[0][0] = 1.0f;
    m->m[1][1] = 1.0f;
    m->m[2][2] = 1.0f;
    m->m[3][3] = 1.0f;
}

static inline void mat4x4_multiplyMatrix(Matrix4x4 *result, const Matrix4x4 *a, const Matrix4x4 *b)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            result->m[i][j] = a->m[i][0] * b->m[0][j] +
                              a->m[i][1] * b->m[1][j] +
                              a->m[i][2] * b->m[2][j] +
                              a->m[i][3] * b->m[3][j];
        }
    }
}

// Create rotation matrix around X axis
static inline void mat4x4_make_rotation_x(Matrix4x4 *m, float angle)
{
    mat4x4_identity(m);
    m->m[1][1] = cosf(angle);
    m->m[1][2] = sinf(angle);
    m->m[2][1] = -sinf(angle);
    m->m[2][2] = cosf(angle);
}
// Create rotation matrix around Z axis
static inline void mat4x4_make_rotation_z(Matrix4x4 *m, float angle)
{
    mat4x4_identity(m);
    m->m[0][0] = cosf(angle);
    m->m[0][1] = sinf(angle);
    m->m[1][0] = -sinf(angle);
    m->m[1][1] = cosf(angle);
}
// Create translation matrix
static inline void mat4x4_make_translation(Matrix4x4 *m, float x, float y, float z)
{
    mat4x4_identity(m);
    m->m[3][0] = x;
    m->m[3][1] = y;
    m->m[3][2] = z;
}
// Create projection matrix
static inline void mat4x4_make_projection(Matrix4x4 *m, float fov, float aspect, float near, float far)
{
    float fovRad = 1.0f / tanf(fov * 0.5f * (3.14159265f / 180.0f));
    memset(m, 0, sizeof(Matrix4x4));
    m->m[0][0] = aspect * fovRad;
    m->m[1][1] = fovRad;
    m->m[2][2] = far / (far - near);
    m->m[3][2] = (-far * near) / (far - near);
    m->m[2][3] = 1.0f;
    m->m[3][3] = 0.0f;
}

static inline void mat4x4_make_rotation_y(Matrix4x4 *m, float angle)
{
    mat4x4_identity(m);
    m->m[0][0] = cosf(angle);
    m->m[0][2] = sinf(angle);
    m->m[2][0] = -sinf(angle);
    m->m[1][1] = 1.0f;
    m->m[2][2] = cosf(angle);
}

// --- OBJ File Loader ---
bool LoadFromObjectFile(Mesh *m, const char *filename)
{
    FILE *f = fopen(filename, "r");
    if (!f)
        return false;
    Vector3 *verts = NULL;
    size_t vertCount = 0;
    m->triangles = NULL;
    m->triangleCount = 0;
    char line[LINE_MAX_LEN];
    while (fgets(line, sizeof(line), f))
    {
        if (line[0] == 'v' && line[1] == ' ')
        {
            Vector3 v;
            if (sscanf(line, "v %f %f %f", &v.x, &v.y, &v.z) == 3)
            {
                verts = realloc(verts, sizeof(Vector3) * (vertCount + 1));
                verts[vertCount++] = v;
            }
        }
        if (line[0] == 'f' && line[1] == ' ')
        {
            int fIdx[3];
            if (sscanf(line, "f %d %d %d", &fIdx[0], &fIdx[1], &fIdx[2]) == 3)
            {
                Triangle t = {
                    .points[0] = verts[fIdx[0] - 1],
                    .points[1] = verts[fIdx[1] - 1],
                    .points[2] = verts[fIdx[2] - 1]};
                m->triangles = realloc(m->triangles, sizeof(Triangle) * (m->triangleCount + 1));
                m->triangles[m->triangleCount++] = t;
            }
        }
    }
    free(verts);
    fclose(f);
    return true;
}

// --- Global State ---
Mesh cubeMesh;
Matrix4x4 projectionMatrix;
Vector3 Camera = {0};
float theta = 0.0f;

// --- Matrix-Vector Multiplication ---
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

// --- Line and Triangle Drawing ---
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

// --- Triangle Rasterization ---
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
        int bx = second_half ? x2 + (int)((x3 - x2) * beta) : x1 + (int)((x2 - x1) * beta);
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

// --- Screen Clear ---
void ClearScreen()
{
    for (int y = 0; y < SCREEN_HEIGHT; y++)
        for (int x = 0; x < SCREEN_WIDTH; x++)
            draw_pixel(x, y, 0, 0, 0, 0); // black
}

// --- Engine Initialization ---
void InitEngine()
{
    if (!LoadFromObjectFile(&cubeMesh, "./teapot.obj"))
    {
        fprintf(stderr, "Failed to load ./teapot.obj\n");
        exit(1);
    }
    float nearPlane = 0.1f, farPlane = 1000.0f, fov = 90.0f;
    float aspectRatio = (float)SCREEN_HEIGHT / (float)SCREEN_WIDTH;
    mat4x4_make_projection(&projectionMatrix, fov, aspectRatio, nearPlane, farPlane);
}

// --- Triangle Raster Info for Sorting ---
typedef struct
{
    Triangle tri;
    float avgDepth;
    float shade;
} TriangleToRaster;

// --- Painter's Algorithm Sort ---
int cmp(const void *a, const void *b)
{
    float z1 = ((const TriangleToRaster *)a)->avgDepth;
    float z2 = ((const TriangleToRaster *)b)->avgDepth;
    if (z1 < z2)
        return 1;
    if (z1 > z2)
        return -1;
    return 0;
}

// --- Main Render Loop ---
void Tick()
{
    static int frameCount = 0;
    printf("Frame: %d\n", frameCount++);
    ClearScreen();
    theta += 0.005f;

    // --- Rotation Matrices ---

    Matrix4x4 rotationZ, rotationX, rotationY, translation, world, temp;
    mat4x4_make_rotation_z(&rotationZ, theta);
    mat4x4_make_rotation_x(&rotationX, theta * 0.5f);
    mat4x4_make_translation(&translation, 0.0f, 0.0f, 8.0f);
    mat4x4_multiplyMatrix(&temp, &rotationZ, &rotationX);
    mat4x4_multiplyMatrix(&world, &temp, &translation);

    // --- Allocate Raster List ---

    TriangleToRaster *trisToRaster = malloc(sizeof(TriangleToRaster) * cubeMesh.triangleCount);
    size_t trisToRasterCount = 0;

    for (size_t i = 0; i < cubeMesh.triangleCount; i++)
    {
        Triangle projected, translated, rotatedZ, rotatedZX;
        for (int j = 0; j < 3; j++)
            MultiplyMatrixVector(&cubeMesh.triangles[i].points[j], &rotatedZ.points[j], &world);
        translated = rotatedZ; // Now 'rotatedZ' is actually world-transformed

        // --- Backface Culling ---

        Vector3 line1 = vec3_sub(translated.points[1], translated.points[0]);
        Vector3 line2 = vec3_sub(translated.points[2], translated.points[0]);
        Vector3 normal = vec3_normalize(vec3_cross(line1, line2));
        if (vec3_dot(normal, vec3_sub(translated.points[0], Camera)) > 0.0f)
            continue; // Cull

        // --- Lighting ---

        Vector3 lightDirection = vec3_normalize((Vector3){0.0f, 0.0f, -1.0f});
        float dotProduct = vec3_dot(normal, lightDirection);
        // --- Project to 2D ---
        for (int j = 0; j < 3; j++)
            MultiplyMatrixVector(&translated.points[j], &projected.points[j], &projectionMatrix);
        for (int j = 0; j < 3; j++)
        {
            projected.points[j].x += 1.0f;
            projected.points[j].y += 1.0f;
            projected.points[j].x *= 0.5f * SCREEN_WIDTH;
            projected.points[j].y *= 0.5f * SCREEN_HEIGHT;
        }

        // --- Shading ---

        float shade = vec3_dot(normal, lightDirection);
        if (shade < AMBIENT_LIGHT_INTENSITY)
            shade = AMBIENT_LIGHT_INTENSITY;
        if (shade > 1.0f)
            shade = 1.0f;

        // --- Depth for Sorting ---

        float avgDepth = (projected.points[0].z + projected.points[1].z + projected.points[2].z) / 3.0f;
        trisToRaster[trisToRasterCount].tri = projected;
        trisToRaster[trisToRasterCount].avgDepth = avgDepth;
        trisToRaster[trisToRasterCount].shade = shade;
        trisToRasterCount++;
    }

    // --- Sort and Rasterize ---

    qsort(trisToRaster, trisToRasterCount, sizeof(TriangleToRaster), cmp);
    for (size_t i = 0; i < trisToRasterCount; i++)
    {
        Triangle *t = &trisToRaster[i].tri;
        float shade = trisToRaster[i].shade;
        uint8_t r = (uint8_t)(255 * shade);
        uint8_t g = (uint8_t)(255 * shade);
        uint8_t b = (uint8_t)(255 * shade);
        FillTriangle(
            (int)t->points[0].x, (int)t->points[0].y,
            (int)t->points[1].x, (int)t->points[1].y,
            (int)t->points[2].x, (int)t->points[2].y, r, g, b, 255);
    }
    free(trisToRaster);
    present_frame();
}

// --- Program Entry ---
int main()
{
    InitEngine();
    cocoa_start(SCREEN_WIDTH, SCREEN_HEIGHT, 0, Tick);
    return 0;
}