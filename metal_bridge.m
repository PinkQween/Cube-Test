#import <Cocoa/Cocoa.h>
#import <Metal/Metal.h>
#import <QuartzCore/CAMetalLayer.h>
#include "metal_bridge.h"

static void (*g_frame_callback)(void) = NULL;
static void (*g_key_callback)(int key, bool pressed) = NULL;

@interface MetalView : NSView
@property (nonatomic, strong) CAMetalLayer *metalLayer;
@end
@implementation MetalView
+ (Class)layerClass { return [CAMetalLayer class]; }
- (instancetype)initWithFrame:(NSRect)frame {
    self = [super initWithFrame:frame];
    if (self) {
        self.wantsLayer = YES;
        self.layer = [CAMetalLayer layer];
        self.metalLayer = (CAMetalLayer *)self.layer;
    }
    return self;
}

- (BOOL)acceptsFirstResponder {
    return YES;
}

- (void)keyDown:(NSEvent *)event {
    if (g_key_callback) {
        g_key_callback((int)event.keyCode, true);
    }
}

- (void)keyUp:(NSEvent *)event {
    if (g_key_callback) {
        g_key_callback((int)event.keyCode, false);
    }
}
@end

static NSWindow *window = nil;
static id<MTLDevice> device = nil;
static id<MTLCommandQueue> commandQueue = nil;
static id<MTLRenderPipelineState> pipelineState = nil;
static id<MTLBuffer> vertexBuffer = nil;
static id<MTLTexture> texture = nil;
static id<MTLTexture> depthTexture = nil;
static int currentVertexCount = 0;

void metal_start(int width, int height, int maxFPS, void (*callback)(void)) {
    @autoreleasepool {
        g_frame_callback = callback;
        device = MTLCreateSystemDefaultDevice();
        if (!device) {
            NSLog(@"Metal is not supported on this device");
            return;
        }
        
        commandQueue = [device newCommandQueue];
        
        // Create shaders inline
        NSString *shaderSource = @"#include <metal_stdlib>\n"
                                 "using namespace metal;\n"
                                 "struct VertexIn {\n"
                                 "    float3 position [[attribute(0)]];\n"
                                 "    float2 uv [[attribute(1)]];\n"
                                 "};\n"
                                 "struct VertexOut {\n"
                                 "    float4 position [[position]];\n"
                                 "    float2 uv;\n"
                                 "};\n"
                                 "vertex VertexOut vertex_main(VertexIn in [[stage_in]]) {\n"
                                 "    VertexOut out;\n"
                                 "    out.position = float4(in.position, 1.0);\n"
                                 "    out.uv = in.uv;\n"
                                 "    return out;\n"
                                 "}\n"
                                 "fragment float4 fragment_main(VertexOut in [[stage_in]], texture2d<float> texture [[texture(0)]]) {\n"
                                 "    constexpr sampler textureSampler(mag_filter::linear, min_filter::linear);\n"
                                 "    return texture.sample(textureSampler, in.uv);\n"
                                 "}\n";
        
        NSError *error = nil;
        id<MTLLibrary> library = [device newLibraryWithSource:shaderSource options:nil error:&error];
        if (!library) {
            NSLog(@"Failed to create library: %@", error);
            return;
        }
        
        id<MTLFunction> vertexFunction = [library newFunctionWithName:@"vertex_main"];
        id<MTLFunction> fragmentFunction = [library newFunctionWithName:@"fragment_main"];
        
        if (!vertexFunction || !fragmentFunction) {
            NSLog(@"Failed to get shader functions");
            return;
        }
        
        // Create render pipeline
        MTLRenderPipelineDescriptor *pipelineDescriptor = [[MTLRenderPipelineDescriptor alloc] init];
        pipelineDescriptor.vertexFunction = vertexFunction;
        pipelineDescriptor.fragmentFunction = fragmentFunction;
        pipelineDescriptor.colorAttachments[0].pixelFormat = MTLPixelFormatBGRA8Unorm;
        
        // Enable depth testing
        pipelineDescriptor.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float;
        
        // Set up vertex descriptor
        MTLVertexDescriptor *vertexDescriptor = [[MTLVertexDescriptor alloc] init];
        vertexDescriptor.attributes[0].format = MTLVertexFormatFloat3;
        vertexDescriptor.attributes[0].offset = 0;
        vertexDescriptor.attributes[0].bufferIndex = 0;
        vertexDescriptor.attributes[1].format = MTLVertexFormatFloat2;
        vertexDescriptor.attributes[1].offset = 3 * sizeof(float);
        vertexDescriptor.attributes[1].bufferIndex = 0;
        vertexDescriptor.layouts[0].stride = 5 * sizeof(float);
        vertexDescriptor.layouts[0].stepRate = 1;
        vertexDescriptor.layouts[0].stepFunction = MTLVertexStepFunctionPerVertex;
        
        pipelineDescriptor.vertexDescriptor = vertexDescriptor;
        
        pipelineState = [device newRenderPipelineStateWithDescriptor:pipelineDescriptor error:&error];
        if (!pipelineState) {
            NSLog(@"Failed to create pipeline state: %@", error);
            return;
        }
        
        // Create vertex buffer
        vertexBuffer = [device newBufferWithLength:100000 * sizeof(float) * 5 options:MTLResourceStorageModeShared];
        
        // Create depth texture
        MTLTextureDescriptor *depthTextureDescriptor = [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:MTLPixelFormatDepth32Float
                                                                                                           width:width
                                                                                                          height:height
                                                                                                       mipmapped:NO];
        depthTextureDescriptor.usage = MTLTextureUsageRenderTarget;
        depthTexture = [device newTextureWithDescriptor:depthTextureDescriptor];
        
        NSApplication *app = [NSApplication sharedApplication];
        [app setActivationPolicy:NSApplicationActivationPolicyRegular];
        NSRect frame = NSMakeRect(0, 0, width, height);
        window = [[NSWindow alloc] initWithContentRect:frame
                                             styleMask:(NSWindowStyleMaskTitled | NSWindowStyleMaskClosable | NSWindowStyleMaskResizable)
                                               backing:NSBackingStoreBuffered
                                                 defer:NO];
        [window setTitle:@"Metal Renderer"]; 
        MetalView *view = [[MetalView alloc] initWithFrame:frame];
        [window setContentView:view];
        
        // Setup Metal layer
        view.metalLayer.device = device;
        view.metalLayer.pixelFormat = MTLPixelFormatBGRA8Unorm;
        view.metalLayer.framebufferOnly = YES;
        view.metalLayer.drawableSize = CGSizeMake(width, height);
        
        [window makeKeyAndOrderFront:nil];
        [view becomeFirstResponder];
        [app activateIgnoringOtherApps:YES];
        
        // Timer for frame callback
        NSTimer *timer = [NSTimer scheduledTimerWithTimeInterval:(maxFPS > 0 ? 1.0/maxFPS : 1.0/60.0)
                                                          repeats:YES
                                                            block:^(NSTimer * _Nonnull timer) {
            if (g_frame_callback) g_frame_callback();
        }];
        [app run];
    }
}

void metal_set_key_callback(void (*key_callback)(int key, bool pressed)) {
    g_key_callback = key_callback;
}

void metal_draw_mesh(const float *vertices, const float *uvs, int vertex_count, const uint8_t *texture_data, int tex_w, int tex_h) {
    if (!device || !commandQueue || !pipelineState) return;
    
    // Update vertex buffer - vertices now contain both position and UV data
    float *vertexData = [vertexBuffer contents];
    memcpy(vertexData, vertices, vertex_count * 5 * sizeof(float));
    currentVertexCount = vertex_count;
    
    // Create texture if needed
    if (!texture && texture_data && tex_w > 0 && tex_h > 0) {
        MTLTextureDescriptor *textureDescriptor = [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:MTLPixelFormatRGBA8Unorm
                                                                                                      width:tex_w
                                                                                                     height:tex_h
                                                                                                  mipmapped:NO];
        texture = [device newTextureWithDescriptor:textureDescriptor];
        
        // Convert RGB to RGBA
        uint8_t *rgba_data = malloc(tex_w * tex_h * 4);
        for (int i = 0; i < tex_w * tex_h; i++) {
            rgba_data[i * 4 + 0] = texture_data[i * 3 + 0]; // R
            rgba_data[i * 4 + 1] = texture_data[i * 3 + 1]; // G
            rgba_data[i * 4 + 2] = texture_data[i * 3 + 2]; // B
            rgba_data[i * 4 + 3] = 255;                     // A
        }
        
        [texture replaceRegion:MTLRegionMake2D(0, 0, tex_w, tex_h)
                   mipmapLevel:0
                     withBytes:rgba_data
                   bytesPerRow:tex_w * 4];
        
        free(rgba_data);
    }
    
    // Get drawable
    MetalView *view = (MetalView *)window.contentView;
    id<CAMetalDrawable> drawable = [view.metalLayer nextDrawable];
    if (!drawable) return;
    
    // Create render pass descriptor
    MTLRenderPassDescriptor *renderPassDescriptor = [MTLRenderPassDescriptor renderPassDescriptor];
    renderPassDescriptor.colorAttachments[0].texture = drawable.texture;
    renderPassDescriptor.colorAttachments[0].loadAction = MTLLoadActionClear;
    renderPassDescriptor.colorAttachments[0].clearColor = MTLClearColorMake(0.2, 0.2, 0.2, 1.0);
    renderPassDescriptor.colorAttachments[0].storeAction = MTLStoreActionStore;
    renderPassDescriptor.depthAttachment.texture = depthTexture;
    renderPassDescriptor.depthAttachment.loadAction = MTLLoadActionClear;
    renderPassDescriptor.depthAttachment.clearDepth = 1.0;
    renderPassDescriptor.depthAttachment.storeAction = MTLStoreActionStore;
    
    // Create command buffer and encoder
    id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
    id<MTLRenderCommandEncoder> renderEncoder = [commandBuffer renderCommandEncoderWithDescriptor:renderPassDescriptor];
    
    [renderEncoder setRenderPipelineState:pipelineState];
    [renderEncoder setVertexBuffer:vertexBuffer offset:0 atIndex:0];
    if (texture) {
        [renderEncoder setFragmentTexture:texture atIndex:0];
    }
    
    // Set viewport
    MTLViewport viewport = {0.0, 0.0, (double)drawable.texture.width, (double)drawable.texture.height, 0.0, 1.0};
    [renderEncoder setViewport:viewport];
    
    [renderEncoder drawPrimitives:MTLPrimitiveTypeTriangle vertexStart:0 vertexCount:currentVertexCount];
    [renderEncoder endEncoding];
    
    [commandBuffer presentDrawable:drawable];
    [commandBuffer commit];
} 