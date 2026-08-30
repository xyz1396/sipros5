#define _GNU_SOURCE

#include <dlfcn.h>
#include <stddef.h>

typedef void (*GlfwWindowHintFunction)(int, int);
typedef void (*GlProcedure)(void);
typedef GlProcedure (*GetProcedureFunction)(const char*);

enum {
    GlfwContextCreationApi = 0x0002200B,
    GlfwNativeContextApi = 0x00036001,
    GlfwEglContextApi = 0x00036002,
};

/*
 * Moba/X may advertise GLX while returning no framebuffer configurations.
 * Switch the context requested by siproswf to EGL, which Mesa can render in
 * software and present through the forwarded X11 display.
 */
void glfwWindowHint(int hint, int value) {
    static GlfwWindowHintFunction next_window_hint;
    if (next_window_hint == NULL) {
        next_window_hint =
            (GlfwWindowHintFunction)dlsym(RTLD_NEXT, "glfwWindowHint");
    }
    if (hint == GlfwContextCreationApi && value == GlfwNativeContextApi) {
        value = GlfwEglContextApi;
    }
    if (next_window_hint != NULL) {
        next_window_hint(hint, value);
    }
}

/*
 * Conda's ImGui OpenGL loader can select its GLX lookup path even when GLFW
 * owns an EGL context. Route that lookup through the already-loaded EGL
 * implementation so all OpenGL entry points use the current EGL context.
 */
static GlProcedure egl_get_procedure(const unsigned char* name) {
    static void* egl_library;
    static GetProcedureFunction get_procedure;
    if (egl_library == NULL) {
        egl_library =
            dlopen("libEGL.so.1", RTLD_LAZY | RTLD_LOCAL | RTLD_NOLOAD);
    }
    if (egl_library == NULL) {
        egl_library = dlopen("libEGL.so.1", RTLD_LAZY | RTLD_LOCAL);
    }
    if (get_procedure == NULL && egl_library != NULL) {
        get_procedure =
            (GetProcedureFunction)dlsym(egl_library, "eglGetProcAddress");
    }
    return get_procedure != NULL
               ? get_procedure((const char*)name)
               : (GlProcedure)0;
}

GlProcedure glXGetProcAddressARB(const unsigned char* name) {
    return egl_get_procedure(name);
}

GlProcedure glXGetProcAddress(const unsigned char* name) {
    return egl_get_procedure(name);
}
