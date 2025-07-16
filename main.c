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
extern int cocoa_is_key_down(int keycode); // key state from cocoa_bridge

#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 800
#define LINE_MAX_LEN 128

// macOS virtual key codes for WASD and arrows
#define KEY_W 13
#define KEY_A 0
#define KEY_S 1
#define KEY_D 2
#define KEY_UP 126
#define KEY_DOWN 125
#define KEY_LEFT 123
#define KEY_RIGHT 124

// Ambient light intensity (minimum light)
static const float AMBIENT_LIGHT_INTENSITY = 0.8f;

// --- Math Types ---
typedef struct
{
    float x, y, z;
} Vector3;

typedef struct
{
    float x, y, z, w;
} Vector4;

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
static inline Vector4 vec4_from_vec3(Vector3 a, float w)
{
    return (Vector4){a.x, a.y, a.z, w};
}
static inline Vector3 vec3_from_vec4(Vector4 a)
{
    return (Vector3){a.x, a.y, a.z};
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
Vector3 Camera = {0.0f, 0.0f, 0.0f};
Vector3 LookDir = {0.0f, 0.0f, 1.0f};
float Yaw = 0.0f;
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
// --- Vector3 Math Utilities ---
// --- Matrix-Vector Multiplication (Vector4) ---
void MultiplyMatrixVector4(const Vector4 *in, Vector4 *out, const Matrix4x4 *m)
{
    float x = in->x, y = in->y, z = in->z, w = in->w;
    out->x = x * m->m[0][0] + y * m->m[1][0] + z * m->m[2][0] + w * m->m[3][0];
    out->y = x * m->m[0][1] + y * m->m[1][1] + z * m->m[2][1] + w * m->m[3][1];
    out->z = x * m->m[0][2] + y * m->m[1][2] + z * m->m[2][2] + w * m->m[3][2];
    out->w = x * m->m[0][3] + y * m->m[1][3] + z * m->m[2][3] + w * m->m[3][3];
}

// --- Camera/View Matrix Utilities ---
Matrix4x4 Matrix_PointAt(Vector3 pos, Vector3 target, Vector3 up)
{
    Vector3 newForward = vec3_normalize(vec3_sub(target, pos));
    Vector3 a = vec3_scale(newForward, vec3_dot(up, newForward));
    Vector3 newUp = vec3_normalize(vec3_sub(up, a));
    Vector3 newRight = vec3_cross(newUp, newForward);
    Matrix4x4 matrix;
    mat4x4_identity(&matrix);
    matrix.m[0][0] = newRight.x;
    matrix.m[0][1] = newRight.y;
    matrix.m[0][2] = newRight.z;
    matrix.m[1][0] = newUp.x;
    matrix.m[1][1] = newUp.y;
    matrix.m[1][2] = newUp.z;
    matrix.m[2][0] = newForward.x;
    matrix.m[2][1] = newForward.y;
    matrix.m[2][2] = newForward.z;
    matrix.m[3][0] = pos.x;
    matrix.m[3][1] = pos.y;
    matrix.m[3][2] = pos.z;
    matrix.m[3][3] = 1.0f;
    return matrix;
}

Matrix4x4 Matrix_QuickInverse(const Matrix4x4 *m)
{
    Matrix4x4 inv;
    mat4x4_identity(&inv);
    inv.m[0][0] = m->m[0][0];
    inv.m[0][1] = m->m[1][0];
    inv.m[0][2] = m->m[2][0];
    inv.m[1][0] = m->m[0][1];
    inv.m[1][1] = m->m[1][1];
    inv.m[1][2] = m->m[2][1];
    inv.m[2][0] = m->m[0][2];
    inv.m[2][1] = m->m[1][2];
    inv.m[2][2] = m->m[2][2];
    inv.m[3][0] = -(m->m[3][0] * inv.m[0][0] + m->m[3][1] * inv.m[1][0] + m->m[3][2] * inv.m[2][0]);
    inv.m[3][1] = -(m->m[3][0] * inv.m[0][1] + m->m[3][1] * inv.m[1][1] + m->m[3][2] * inv.m[2][1]);
    inv.m[3][2] = -(m->m[3][0] * inv.m[0][2] + m->m[3][1] * inv.m[1][2] + m->m[3][2] * inv.m[2][2]);
    return inv;
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

// --- Triangle Clipping Utilities ---
// Returns number of output triangles (0, 1, or 2)
static inline float dist_to_plane(const Vector3 *p, Vector3 plane_p, Vector3 plane_n)
{
    return plane_n.x * p->x + plane_n.y * p->y + plane_n.z * p->z - vec3_dot(plane_n, plane_p);
}
int Triangle_ClipAgainstPlane(Vector3 plane_p, Vector3 plane_n, Triangle *in_tri, Triangle *out_tri1, Triangle *out_tri2)
{
    plane_n = vec3_normalize(plane_n);
    Vector3 *inside_points[3];
    int nInsidePointCount = 0;
    Vector3 *outside_points[3];
    int nOutsidePointCount = 0;
    float d0 = dist_to_plane(&in_tri->points[0], plane_p, plane_n);
    float d1 = dist_to_plane(&in_tri->points[1], plane_p, plane_n);
    float d2 = dist_to_plane(&in_tri->points[2], plane_p, plane_n);
    float t;
    Vector3 dir;
    if (d0 >= 0)
    {
        inside_points[nInsidePointCount++] = &in_tri->points[0];
    }
    else
    {
        outside_points[nOutsidePointCount++] = &in_tri->points[0];
    }
    if (d1 >= 0)
    {
        inside_points[nInsidePointCount++] = &in_tri->points[1];
    }
    else
    {
        outside_points[nOutsidePointCount++] = &in_tri->points[1];
    }
    if (d2 >= 0)
    {
        inside_points[nInsidePointCount++] = &in_tri->points[2];
    }
    else
    {
        outside_points[nOutsidePointCount++] = &in_tri->points[2];
    }
    if (nInsidePointCount == 0)
        return 0;
    if (nInsidePointCount == 3)
    {
        *out_tri1 = *in_tri;
        return 1;
    }
    if (nInsidePointCount == 1 && nOutsidePointCount == 2)
    {
        out_tri1->points[0] = *inside_points[0];
        dir = vec3_sub(*outside_points[0], *inside_points[0]);
        t = -dist_to_plane(inside_points[0], plane_p, plane_n) / (dist_to_plane(outside_points[0], plane_p, plane_n) - dist_to_plane(inside_points[0], plane_p, plane_n));
        out_tri1->points[1] = vec3_add(*inside_points[0], vec3_scale(dir, t));
        dir = vec3_sub(*outside_points[1], *inside_points[0]);
        t = -dist_to_plane(inside_points[0], plane_p, plane_n) / (dist_to_plane(outside_points[1], plane_p, plane_n) - dist_to_plane(inside_points[0], plane_p, plane_n));
        out_tri1->points[2] = vec3_add(*inside_points[0], vec3_scale(dir, t));
        return 1;
    }
    if (nInsidePointCount == 2 && nOutsidePointCount == 1)
    {
        out_tri1->points[0] = *inside_points[0];
        out_tri1->points[1] = *inside_points[1];
        dir = vec3_sub(*outside_points[0], *inside_points[0]);
        t = -dist_to_plane(inside_points[0], plane_p, plane_n) / (dist_to_plane(outside_points[0], plane_p, plane_n) - dist_to_plane(inside_points[0], plane_p, plane_n));
        out_tri1->points[2] = vec3_add(*inside_points[0], vec3_scale(dir, t));
        out_tri2->points[0] = *inside_points[1];
        out_tri2->points[1] = out_tri1->points[2];
        dir = vec3_sub(*outside_points[0], *inside_points[1]);
        t = -dist_to_plane(inside_points[1], plane_p, plane_n) / (dist_to_plane(outside_points[0], plane_p, plane_n) - dist_to_plane(inside_points[1], plane_p, plane_n));
        out_tri2->points[2] = vec3_add(*inside_points[1], vec3_scale(dir, t));
        return 2;
    }
    return 0;
}

// --- Main Render Loop ---
void Tick()
{
    static int frameCount = 0;
    printf("Frame: %d\n", frameCount++);
    ClearScreen();
       theta += 0.005f;

    // --- Camera Movement/Input ---
    float moveSpeed = 0.1f;
    float rotSpeed = 0.03f;
    Vector3 forward = vec3_normalize(LookDir);
    Vector3 right = vec3_normalize(vec3_cross(forward, (Vector3){0, 1, 0}));
    if (cocoa_is_key_down(KEY_W) || cocoa_is_key_down(KEY_UP))
        Camera = vec3_add(Camera, vec3_scale(forward, moveSpeed));
    if (cocoa_is_key_down(KEY_S) || cocoa_is_key_down(KEY_DOWN))
        Camera = vec3_sub(Camera, vec3_scale(forward, moveSpeed));
    if (cocoa_is_key_down(KEY_A) || cocoa_is_key_down(KEY_LEFT))
        Camera = vec3_sub(Camera, vec3_scale(right, moveSpeed));
    if (cocoa_is_key_down(KEY_D) || cocoa_is_key_down(KEY_RIGHT))
        Camera = vec3_add(Camera, vec3_scale(right, moveSpeed));
    if (cocoa_is_key_down(12))
        Yaw -= rotSpeed; // Q
    if (cocoa_is_key_down(14))
        Yaw += rotSpeed; // E
    LookDir.x = sinf(Yaw);
    LookDir.z = cosf(Yaw);

    // --- Rotation Matrices ---
    Matrix4x4 rotationZ, rotationX, translation, world, temp;
    mat4x4_make_rotation_z(&rotationZ, theta);
    mat4x4_make_rotation_x(&rotationX, theta * 0.5f + 3.14159265f);
    mat4x4_make_translation(&translation, 0.0f, 0.0f, 8.0f);
    mat4x4_multiplyMatrix(&temp, &rotationZ, &rotationX);
    mat4x4_multiplyMatrix(&world, &temp, &translation);

    // --- Camera/View Matrix ---
    Vector3 up = {0, 1, 0};
    Vector3 target = vec3_add(Camera, LookDir);
    Matrix4x4 matCamera = Matrix_PointAt(Camera, target, up);
    Matrix4x4 matView = Matrix_QuickInverse(&matCamera);

    // --- Allocate Raster List ---
    TriangleToRaster *trisToRaster = malloc(sizeof(TriangleToRaster) * cubeMesh.triangleCount * 8);
    size_t trisToRasterCount = 0;

    for (size_t i = 0; i < cubeMesh.triangleCount; i++)
    {
        // --- 1. World Transform ---
        Vector4 triInput[3];
        for (int j = 0; j < 3; j++)
        {
            triInput[j] = (Vector4){cubeMesh.triangles[i].points[j].x, cubeMesh.triangles[i].points[j].y, cubeMesh.triangles[i].points[j].z, 1.0f};
            MultiplyMatrixVector4(&triInput[j], &triInput[j], &world);
        }
        // --- 2. Calculate Normal ---
        Vector3 line1 = vec3_sub((Vector3){triInput[1].x, triInput[1].y, triInput[1].z}, (Vector3){triInput[0].x, triInput[0].y, triInput[0].z});
        Vector3 line2 = vec3_sub((Vector3){triInput[2].x, triInput[2].y, triInput[2].z}, (Vector3){triInput[0].x, triInput[0].y, triInput[0].z});
        Vector3 normal = vec3_normalize(vec3_cross(line1, line2));
        // --- 3. Backface Culling ---
        Vector3 cameraRay = vec3_sub((Vector3){triInput[0].x, triInput[0].y, triInput[0].z}, Camera);
        if (vec3_dot(normal, cameraRay) < 0.0f)
        {
            // --- 4. Lighting ---
            Vector3 light_direction = vec3_normalize((Vector3){0.0f, 1.0f, -1.0f});
            float shade = fmaxf(0.1f, vec3_dot(light_direction, normal));
            // --- 5. View Transform ---
            Vector4 triViewed[3];
            for (int j = 0; j < 3; j++)
                MultiplyMatrixVector4(&triInput[j], &triViewed[j], &matView);
            // --- 6. Near Plane Clipping ---
            Triangle triViewedSimple = {
                .points[0] = (Vector3){triViewed[0].x, triViewed[0].y, triViewed[0].z},
                .points[1] = (Vector3){triViewed[1].x, triViewed[1].y, triViewed[1].z},
                .points[2] = (Vector3){triViewed[2].x, triViewed[2].y, triViewed[2].z}};
            Triangle clipped[2];
            int nClippedTriangles = Triangle_ClipAgainstPlane((Vector3){0, 0, 0.1f}, (Vector3){0, 0, 1}, &triViewedSimple, &clipped[0], &clipped[1]);
            for (int n = 0; n < nClippedTriangles; n++)
            {
                // --- 7. Projection ---
                Vector4 triProjected[3];
                for (int j = 0; j < 3; j++)
                {
                    Vector4 v = {clipped[n].points[j].x, clipped[n].points[j].y, clipped[n].points[j].z, 1.0f};
                    MultiplyMatrixVector4(&v, &triProjected[j], &projectionMatrix);
                }
                // --- 8. Perspective Divide ---
                for (int j = 0; j < 3; j++)
                {
                    if (fabsf(triProjected[j].w) > 1e-5f)
                    {
                        triProjected[j].x /= triProjected[j].w;
                        triProjected[j].y /= triProjected[j].w;
                        triProjected[j].z /= triProjected[j].w;
                    }
                }
                // --- 9. Screen Transform ---
                for (int j = 0; j < 3; j++)
                {
                    triProjected[j].x *= -1.0f;
                    triProjected[j].y *= -1.0f;
                    triProjected[j].x += 1.0f;
                    triProjected[j].y += 1.0f;
                    triProjected[j].x *= 0.5f * SCREEN_WIDTH;
                    triProjected[j].y *= 0.5f * SCREEN_HEIGHT;
                }
                // --- 10. Store for sorting ---
                Triangle triToRaster = {
                    .points[0] = (Vector3){triProjected[0].x, triProjected[0].y, triProjected[0].z},
                    .points[1] = (Vector3){triProjected[1].x, triProjected[1].y, triProjected[1].z},
                    .points[2] = (Vector3){triProjected[2].x, triProjected[2].y, triProjected[2].z}};
                trisToRaster[trisToRasterCount].tri = triToRaster;
                trisToRaster[trisToRasterCount].avgDepth = (triProjected[0].z + triProjected[1].z + triProjected[2].z) / 3.0f;
                trisToRaster[trisToRasterCount].shade = shade;
                trisToRasterCount++;
            }
        }
    }
    // --- Painter's Algorithm Sort ---
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

// (Optional) Input handling stubs for camera movement can be added here if platform allows