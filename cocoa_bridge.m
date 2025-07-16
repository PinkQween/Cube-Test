#import <Cocoa/Cocoa.h>
#include <stdint.h>
#include <stdlib.h>

static uint8_t *pixelBuffer = NULL;
static int bufferWidth = 0, bufferHeight = 0;
static int target_fps = 60;
static void (*frame_callback)(void) = NULL;

// --- Minimal Key State Tracking ---
#define COCOA_KEY_MAX 256
static uint8_t key_state[COCOA_KEY_MAX] = {0};

// Expose to C
int cocoa_is_key_down(int keycode) {
    if (keycode < 0 || keycode >= COCOA_KEY_MAX) return 0;
    return key_state[keycode];
}

@interface CView : NSView
@end

@implementation CView
- (BOOL)isFlipped { return YES; }

- (void)drawRect:(NSRect)dirtyRect {
    if (!pixelBuffer) return;

    CGColorSpaceRef cs = CGColorSpaceCreateDeviceRGB();
    CGContextRef ctx = CGBitmapContextCreate(pixelBuffer,
                                             bufferWidth,
                                             bufferHeight,
                                             8,
                                             bufferWidth * 4,
                                             cs,
                                             kCGBitmapByteOrder32Host | kCGImageAlphaPremultipliedLast);
    CGImageRef img = CGBitmapContextCreateImage(ctx);
    CGContextDrawImage([[NSGraphicsContext currentContext] CGContext], self.bounds, img);
    CGImageRelease(img);
    CGContextRelease(ctx);
    CGColorSpaceRelease(cs);
}

- (BOOL)acceptsFirstResponder { return YES; }
- (void)keyDown:(NSEvent *)event {
    unsigned short code = [event keyCode];
    if (code < COCOA_KEY_MAX) key_state[code] = 1;
}
- (void)keyUp:(NSEvent *)event {
    unsigned short code = [event keyCode];
    if (code < COCOA_KEY_MAX) key_state[code] = 0;
}
@end

@interface CAppDelegate : NSObject <NSApplicationDelegate>
@end

@implementation CAppDelegate {
    NSTimer *timer;
    NSWindow *window;
    CView *view;
}

- (void)applicationDidFinishLaunching:(NSNotification *)notification {
    window = [[NSWindow alloc]
              initWithContentRect:NSMakeRect(100, 100, bufferWidth, bufferHeight)
                        styleMask:(NSWindowStyleMaskTitled | NSWindowStyleMaskClosable)
                          backing:NSBackingStoreBuffered
                            defer:NO];
    view = [[CView alloc] initWithFrame:NSMakeRect(0, 0, bufferWidth, bufferHeight)];
    [window setContentView:view];
    [window makeKeyAndOrderFront:nil];
    [NSApp activateIgnoringOtherApps:YES];
    [window makeFirstResponder:view];
    [view becomeFirstResponder];

    double interval = (target_fps > 0) ? (1.0 / target_fps) : 0.0;

    if (interval > 0.0) {
        timer = [NSTimer scheduledTimerWithTimeInterval:interval
                                                 target:self
                                               selector:@selector(tick)
                                               userInfo:nil
                                                repeats:YES];
    } else {
        // Unlimited FPS simulation — high-speed timer
        timer = [NSTimer scheduledTimerWithTimeInterval:(1.0 / 1000.0)
                                                 target:self
                                               selector:@selector(tick)
                                               userInfo:nil
                                                repeats:YES];
    }
}

- (void)tick {
    if (frame_callback) frame_callback();
    [view setNeedsDisplay:YES];
}

@end

// C-callable entry point
void cocoa_start(int width, int height, int maxFPS, void (*callback)(void)) {
    bufferWidth = width;
    bufferHeight = height;
    frame_callback = callback;
    target_fps = maxFPS;
    pixelBuffer = calloc(width * height * 4, 1);

    @autoreleasepool {
        NSApplication *app = [NSApplication sharedApplication];
        CAppDelegate *delegate = [[CAppDelegate alloc] init];
        [app setDelegate:delegate];
        [app setActivationPolicy:NSApplicationActivationPolicyRegular];
        [app run];
    }
}

void draw_pixel(int x, int y, uint8_t a, uint8_t b, uint8_t g, uint8_t r) {
    if (!pixelBuffer || x < 0 || y < 0 || x >= bufferWidth || y >= bufferHeight) return;
    int offset = (y * bufferWidth + x) * 4;
    pixelBuffer[offset + 0] = r;
    pixelBuffer[offset + 1] = g;
    pixelBuffer[offset + 2] = b;
    pixelBuffer[offset + 3] = a;
}

void present_frame(void) {
    // no-op, handled by view refresh after tick
}
