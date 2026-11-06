/*
 * tc_unwind_safe.c -- may TimeConstrained's SIGPROF handler unwind the code
 * it interrupted?  See tc_unwind_safe.h.
 *
 * TimeConstrained preempts a long-running C builtin by siglongjmp'ing out of
 * an asynchronous SIGPROF handler.  The signal can land anywhere, and a jump
 * out of malloc / realloc / free leaves the allocator's lock owned by a frame
 * that no longer exists: the next allocation traps (macOS libmalloc,
 * _os_unfair_lock_recursive_abort) or deadlocks (glibc).  A stdio call is the
 * same hazard with the FILE lock.  The handler therefore asks, before
 * jumping, whether any frame of the interrupted stack is executing inside the
 * C library, and if so defers and tries again a millisecond later.  The walk
 * follows the frame-pointer chain, which macOS guarantees on both arm64 and
 * x86_64 (its system libraries keep frame pointers).  glibc's allocator lives
 * entirely in libc itself, so on Linux the interrupted PC alone already
 * catches the allocator; the chain walk is a bonus where frame pointers
 * exist.
 *
 * Everything the handler touches is fixed at init: the code ranges of the
 * library images (found by image name on macOS, by the image that holds
 * `malloc` / `fwrite` on Linux) and the main thread's stack bounds, which
 * keep the chain walk from reading outside the stack.  On any other platform,
 * or before init, the check answers "clear" and the handler behaves as it
 * always did.
 */
#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE   /* dl_iterate_phdr, pthread_getattr_np, REG_RIP */
#endif

#include "tc_unwind_safe.h"

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(__APPLE__) && (defined(__x86_64__) || defined(__aarch64__) || defined(__arm64__))
#define TC_UNWIND_MACOS 1
#include <pthread.h>
#include <sys/ucontext.h>
#include <mach-o/dyld.h>
#include <mach-o/loader.h>
#elif defined(__linux__) && (defined(__x86_64__) || defined(__aarch64__))
#define TC_UNWIND_LINUX 1
#include <pthread.h>
#include <ucontext.h>
#include <link.h>
#endif

#if defined(TC_UNWIND_MACOS) || defined(TC_UNWIND_LINUX)

#define TC_MAX_RANGES 8
#define TC_MAX_FRAMES 128

typedef struct { uintptr_t lo, hi; } CodeRange;

static CodeRange g_ranges[TC_MAX_RANGES];
static int       g_nranges = 0;
static uintptr_t g_stack_lo = 0, g_stack_hi = 0;   /* main thread's stack */
static int       g_ready = 0;

static void add_range(uintptr_t lo, uintptr_t hi) {
    if (g_nranges < TC_MAX_RANGES && hi > lo) {
        g_ranges[g_nranges].lo = lo;
        g_ranges[g_nranges].hi = hi;
        g_nranges++;
    }
}

static int in_ranges(uintptr_t pc) {
    for (int i = 0; i < g_nranges; i++)
        if (pc >= g_ranges[i].lo && pc < g_ranges[i].hi) return 1;
    return 0;
}

#ifdef TC_UNWIND_MACOS
/* The libsystem images an unwind must not leave: the allocator, the
 * platform layer it calls into (memmove, os_unfair_lock), libc (stdio), and
 * pthread (mutexes) -- and, in a sanitizer build, the sanitizer runtime,
 * whose malloc replaces the system one and holds its own locks (a jump out
 * of ASan's allocator left its stack-depot lock taken, and the next
 * allocation spun forever). */
static const char* const k_images[] = {
    "/libsystem_malloc.dylib", "/libsystem_platform.dylib",
    "/libsystem_c.dylib", "/libsystem_pthread.dylib",
    "/libasan.", "/libubsan.", "/libclang_rt.asan_", "/libclang_rt.ubsan_"
};

static void find_ranges(void) {
    uint32_t n = _dyld_image_count();
    for (uint32_t i = 0; i < n; i++) {
        const char* name = _dyld_get_image_name(i);
        if (!name) continue;
        int want = 0;
        for (size_t k = 0; k < sizeof k_images / sizeof k_images[0]; k++)
            if (strstr(name, k_images[k])) { want = 1; break; }
        if (!want) continue;
        const struct mach_header_64* h =
            (const struct mach_header_64*)_dyld_get_image_header(i);
        if (!h || h->magic != MH_MAGIC_64) continue;
        intptr_t slide = _dyld_get_image_vmaddr_slide(i);
        const struct load_command* lc = (const struct load_command*)(h + 1);
        for (uint32_t c = 0; c < h->ncmds; c++) {
            if (lc->cmd == LC_SEGMENT_64) {
                const struct segment_command_64* sg = (const struct segment_command_64*)lc;
                if (strncmp(sg->segname, "__TEXT", 16) == 0)
                    add_range((uintptr_t)(sg->vmaddr + slide),
                              (uintptr_t)(sg->vmaddr + slide + sg->vmsize));
            }
            lc = (const struct load_command*)((const char*)lc + lc->cmdsize);
        }
    }
}

static void find_stack(void) {
    pthread_t self = pthread_self();
    uintptr_t top = (uintptr_t)pthread_get_stackaddr_np(self);
    size_t size = pthread_get_stacksize_np(self);
    g_stack_hi = top;
    g_stack_lo = top - size;
}

/* pc, fp and sp of the interrupted context. */
static void read_context(const void* uc_, uintptr_t* pc, uintptr_t* fp,
                         uintptr_t* sp, uintptr_t* lr) {
    const ucontext_t* uc = (const ucontext_t*)uc_;
#if defined(__x86_64__)
    *pc = (uintptr_t)uc->uc_mcontext->__ss.__rip;
    *fp = (uintptr_t)uc->uc_mcontext->__ss.__rbp;
    *sp = (uintptr_t)uc->uc_mcontext->__ss.__rsp;
    *lr = 0;
#else
    *pc = (uintptr_t)uc->uc_mcontext->__ss.__pc;
    *fp = (uintptr_t)uc->uc_mcontext->__ss.__fp;
    *sp = (uintptr_t)uc->uc_mcontext->__ss.__sp;
    *lr = (uintptr_t)uc->uc_mcontext->__ss.__lr;
#endif
}
#endif /* TC_UNWIND_MACOS */

#ifdef TC_UNWIND_LINUX
/* The object whose executable segment holds `malloc` (glibc: libc.so, which
 * also holds stdio and, since 2.34, pthread). */
static int phdr_cb(struct dl_phdr_info* info, size_t size, void* data) {
    (void)size;
    const uintptr_t* targets = (const uintptr_t*)data;
    for (int k = 0; k < info->dlpi_phnum; k++) {
        const ElfW(Phdr)* ph = &info->dlpi_phdr[k];
        if (ph->p_type != PT_LOAD || !(ph->p_flags & PF_X)) continue;
        uintptr_t lo = (uintptr_t)(info->dlpi_addr + ph->p_vaddr);
        uintptr_t hi = lo + (uintptr_t)ph->p_memsz;
        for (int t = 0; targets[t]; t++)
            if (targets[t] >= lo && targets[t] < hi) { add_range(lo, hi); break; }
    }
    return 0;
}

static void find_ranges(void) {
    uintptr_t targets[] = { (uintptr_t)&malloc, (uintptr_t)&fwrite,
                            (uintptr_t)&pthread_mutex_lock, 0 };
    dl_iterate_phdr(phdr_cb, targets);
}

static void find_stack(void) {
    pthread_attr_t a;
    void* addr = NULL;
    size_t size = 0;
    if (pthread_getattr_np(pthread_self(), &a) != 0) return;
    if (pthread_attr_getstack(&a, &addr, &size) == 0) {
        g_stack_lo = (uintptr_t)addr;
        g_stack_hi = (uintptr_t)addr + size;
    }
    pthread_attr_destroy(&a);
}

static void read_context(const void* uc_, uintptr_t* pc, uintptr_t* fp,
                         uintptr_t* sp, uintptr_t* lr) {
    const ucontext_t* uc = (const ucontext_t*)uc_;
#if defined(__x86_64__)
    *pc = (uintptr_t)uc->uc_mcontext.gregs[REG_RIP];
    *fp = (uintptr_t)uc->uc_mcontext.gregs[REG_RBP];
    *sp = (uintptr_t)uc->uc_mcontext.gregs[REG_RSP];
    *lr = 0;
#else
    *pc = (uintptr_t)uc->uc_mcontext.pc;
    *fp = (uintptr_t)uc->uc_mcontext.regs[29];
    *sp = (uintptr_t)uc->uc_mcontext.sp;
    *lr = (uintptr_t)uc->uc_mcontext.regs[30];
#endif
}
#endif /* TC_UNWIND_LINUX */

void tc_unwind_safe_init(void) {
    if (g_ready) return;
    find_ranges();
    find_stack();
    g_ready = g_nranges > 0;
}

int tc_interrupted_in_libc(const void* ucontext) {
    if (!g_ready || !ucontext) return 0;
    uintptr_t pc, fp, sp, lr;
    read_context(ucontext, &pc, &fp, &sp, &lr);
    /* Not on the main thread's stack: a worker thread took the signal.  The
     * jump targets the main thread's frame, so it must not be taken here. */
    if (g_stack_hi && (sp < g_stack_lo || sp >= g_stack_hi)) return 1;
    if (in_ranges(pc) || in_ranges(lr)) return 1;
    /* x86-64: a leaf the signal caught before it pushed a frame has its
     * caller's return address at [sp]. */
#if defined(__x86_64__)
    if (sp && (sp & 7u) == 0 && in_ranges(*(const uintptr_t*)sp)) return 1;
#endif
    /* Frame-pointer chain: [fp] = caller's fp, [fp + 8] = return address.
     * Each fp must lie on the stack, above the last, and be aligned. */
    if (!g_stack_hi) return 0;          /* bounds unknown: do not walk */
    uintptr_t prev = sp;
    for (int i = 0; i < TC_MAX_FRAMES; i++) {
        if (fp == 0 || (fp & 7u) != 0 || fp < prev) break;
        if (fp < g_stack_lo || fp + 16 > g_stack_hi) break;
        uintptr_t ret = ((const uintptr_t*)fp)[1];
        if (in_ranges(ret)) return 1;
        prev = fp;
        fp = ((const uintptr_t*)fp)[0];
    }
    return 0;
}

#else /* no context access on this platform */

void tc_unwind_safe_init(void) {}
int tc_interrupted_in_libc(const void* ucontext) { (void)ucontext; return 0; }

#endif
