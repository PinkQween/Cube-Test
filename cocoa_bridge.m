#import <Cocoa/Cocoa.h>
#include <stdint.h>
#include <stdlib.h>

static uint8_t *pixelBuffer = NULL;
static int bufferWidth = 0, bufferHeight = 0;
static void (*frame_callback)(void) = NULL;

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

    timer = [NSTimer scheduledTimerWithTimeInterval:1.0/60.0
                                              target:self
                                            selector:@selector(tick)
                                            userInfo:nil
                                             repeats:YES];
}
- (void)tick {
    if (frame_callback) frame_callback();
    [view setNeedsDisplay:YES];
}
@end

// C-callable entry point
void cocoa_start(int width, int height, void (*callback)(void)) {
    bufferWidth = width;
    bufferHeight = height;
    frame_callback = callback;
    pixelBuffer = calloc(width * height * 4, 1);

    @autoreleasepool {
        NSApplication *app = [NSApplication sharedApplication];
        CAppDelegate *delegate = [[CAppDelegate alloc] init];
        [app setDelegate:delegate];
        [app setActivationPolicy:NSApplicationActivationPolicyRegular];
        [app run];
    }
}

void draw_pixel(int x, int y, uint8_t r, uint8_t g, uint8_t b, uint8_t a) {
    if (!pixelBuffer || x < 0 || y < 0 || x >= bufferWidth || y >= bufferHeight) return;
    int offset = (y * bufferWidth + x) * 4;
    pixelBuffer[offset + 0] = r;
    pixelBuffer[offset + 1] = g;
    pixelBuffer[offset + 2] = b;
    pixelBuffer[offset + 3] = a;
}

void present_frame(void) {
    // handled by view refresh after tick
}