/* gl_probe.c -- can this process really render with raylib? Asked once, in a child process.
 *
 * WHY A CHILD. Raylib 5.5's InitWindow cannot fail safely. On a Mac where the window can be
 * created but no usable OpenGL context comes with it (seen from a Claude Code tool process: a GUI
 * session exists, CGL even reports an accelerated pixel format), InitWindow goes on to rlglInit,
 * whose first rlLoadTexture calls a GL entry point that was never loaded -- a jump to address 0:
 *
 *     frame #0: 0x0000000000000000
 *     frame #1: libraylib.550.dylib`rlLoadTexture + 65
 *     frame #2: libraylib.550.dylib`rlglInit + 54
 *     frame #3: libraylib.550.dylib`InitWindow + 388
 *     frame #4: Mathilda`graphics_render_rgba ... render.c
 *
 * Nothing in-process can see that coming: GLFW's init and window creation both "succeeded", and
 * the crash is inside the library before control returns. Every environment heuristic tried (an
 * Aqua session, a pixel format) said yes in that process. The only honest test is to DO it -- open
 * the window, render to a texture, read it back -- somewhere a crash costs nothing.
 *
 * HOW. The parent posix_spawn()s its OWN executable with MATHILDA_GL_PROBE=1 in the environment.
 * A constructor in this file (so it runs in Mathilda and in every test binary linking the kernel,
 * before main) sees the variable, performs exactly the export path's raylib sequence, and _exit()s
 * 0 only if all of it worked. spawn rather than fork: a forked child of a process that has touched
 * CoreFoundation may not use it again without exec, and Cocoa is all CoreFoundation. The answer is
 * cached for the session, so the cost (one short-lived process, a few hundred ms) is paid at most
 * once, and only by the first raster export. A probe that does not finish in 20 s is killed and
 * counts as "no". macOS only; elsewhere the caller's own check stands. */
#if defined(__APPLE__) && defined(USE_GRAPHICS) && !defined(_WIN32)
#define _DARWIN_C_SOURCE 1          /* posix_spawn, waitpid, kill: before any #include */
#endif

#include "render_common.h"

#if defined(__APPLE__) && defined(USE_GRAPHICS) && !defined(_WIN32)

#include <raylib.h>
#include <spawn.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <mach-o/dyld.h>
#include "mtime.h"

extern char** environ;

#define GL_PROBE_ENV "MATHILDA_GL_PROBE"

/* The child: the same calls graphics_render_rgba makes, in the same order, on a tiny canvas. */
static int gl_probe_run(void) {
    SetTraceLogLevel(LOG_NONE);
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(64, 64, "mathilda-gl-probe");
    if (!IsWindowReady()) return 2;
    SetWindowState(FLAG_WINDOW_HIDDEN);
    RenderTexture2D rt = LoadRenderTexture(64, 64);
    if (rt.id == 0) { CloseWindow(); return 3; }
    BeginTextureMode(rt);
    ClearBackground(RAYWHITE);
    DrawCircle(32, 32, 10.0f, RED);
    EndTextureMode();
    Image im = LoadImageFromTexture(rt.texture);
    int ok = im.data != NULL && im.width == 64 && im.height == 64;
    UnloadImage(im);
    UnloadRenderTexture(rt);
    CloseWindow();
    return ok ? 0 : 4;
}

__attribute__((constructor))
static void gl_probe_child_entry(void) {
    const char* v = getenv(GL_PROBE_ENV);
    if (v && strcmp(v, "1") == 0) _exit(gl_probe_run());
}

bool gl_probe_usable(void) {
    static int cached = -1;
    if (cached >= 0) return cached != 0;
    cached = 0;

    char path[4096];
    uint32_t sz = (uint32_t)sizeof path;
    if (_NSGetExecutablePath(path, &sz) != 0) return false;

    /* environ + MATHILDA_GL_PROBE=1 */
    size_t ne = 0;
    while (environ && environ[ne]) ne++;
    char** envp = malloc(sizeof(char*) * (ne + 2));
    if (!envp) return false;
    for (size_t i = 0; i < ne; i++) envp[i] = environ[i];
    envp[ne] = (char*)GL_PROBE_ENV "=1";
    envp[ne + 1] = NULL;

    posix_spawn_file_actions_t fa;
    posix_spawn_file_actions_init(&fa);
    posix_spawn_file_actions_addopen(&fa, 0, "/dev/null", O_RDONLY, 0);
    posix_spawn_file_actions_addopen(&fa, 1, "/dev/null", O_WRONLY, 0);
    posix_spawn_file_actions_addopen(&fa, 2, "/dev/null", O_WRONLY, 0);
    char* argv[2];
    argv[0] = path;
    argv[1] = NULL;
    pid_t pid = 0;
    int rc = posix_spawn(&pid, path, &fa, NULL, argv, envp);
    posix_spawn_file_actions_destroy(&fa);
    free(envp);
    if (rc != 0) return false;

    int st = 0;
    double deadline = mth_mono_seconds() + 20.0;
    for (;;) {
        pid_t r = waitpid(pid, &st, WNOHANG);
        if (r == pid) break;
        if (r < 0) return false;
        if (mth_mono_seconds() > deadline) {
            kill(pid, SIGKILL);
            waitpid(pid, &st, 0);
            return false;
        }
        mth_sleep_seconds(0.005);
    }
    cached = (WIFEXITED(st) && WEXITSTATUS(st) == 0) ? 1 : 0;
    return cached != 0;
}

#else

/* Not macOS (or no renderer): nothing to probe here; the caller's own check stands. */
bool gl_probe_usable(void) { return true; }

#endif
