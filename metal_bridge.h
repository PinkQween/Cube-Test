#ifndef METAL_BRIDGE_H
#define METAL_BRIDGE_H
#ifdef __cplusplus
extern "C"
{
#endif
    void metal_start(int width, int height, int maxFPS, void (*callback)(void));
    void metal_draw_mesh(const float *vertices, const float *uvs, int vertex_count, const uint8_t *texture, int tex_w, int tex_h);

    // Keyboard input handling
    void metal_set_key_callback(void (*key_callback)(int key, bool pressed));

#ifdef __cplusplus
}
#endif
#endif // METAL_BRIDGE_H