#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define WIDTH 400
#define HEIGHT 400

extern void cocoa_start(int width, int height, void (*callback)(void));
extern void draw_pixel(int x, int y, uint8_t r, uint8_t g, uint8_t b, uint8_t a);
extern void present_frame(void);

int frame = 0;

void draw_line(int x1, int y1, int x2, int y2, uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
    // Simple Bresenham's line algorithm
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    int sx = (x1 < x2) ? 1 : -1;
    int sy = (y1 < y2) ? 1 : -1;
    int err = dx - dy;

    while (1)
    {
        draw_pixel(x1, y1, r, g, b, a);
        if (x1 == x2 && y1 == y2) break;
        int e2 = err * 2;
        if (e2 > -dy) { err -= dy; x1 += sx; }
        if (e2 < dx) { err += dx; y1 += sy; }
    }
}

void update_frame()
{
    for (int y = 0; y < HEIGHT; y++)
    {
        for (int x = 0; x < WIDTH; x++)
        {
            draw_pixel(x, y, (x + frame) % 255, (y + frame) % 255, 128, 255);
        }
    }
    present_frame();
    frame++;
}

int main()
{
    cocoa_start(WIDTH, HEIGHT, update_frame);
    return 0;
}
