

#pragma once
#include <stdio.h>
#include "ndarray.h"   /* is_packed_list / ndarray_to_nested_list */
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <stdlib.h>
#include <fcntl.h>

extern struct Expr* parse_expression(const char*);
extern struct Expr* evaluate(struct Expr*);
extern char* expr_to_string(struct Expr*);
extern char* expr_to_string_fullform(struct Expr*);
extern void expr_free(struct Expr*);

/* Load the src/internal/init.m bootstrap (deriv.m rules, CRC integral tables) the
 * way the REPL / -file mode does, so DSolve / Integrate tests exercise the same
 * environment as production.  Without it the suite tests a non-production DSolve
 * (which, e.g., took a different, working path where init.m makes it hang — the
 * Lagrange/NthAlgebraic case).  chdir's to the repo root so Get[] resolves. */
static inline void test_load_init_m(void) {
    const char* ups[] = { ".", "..", "../..", "../../..", "../../../..", NULL };
    for (int i = 0; ups[i]; i++) {
        char path[256];
        snprintf(path, sizeof(path), "%s/src/internal/init.m", ups[i]);
        if (access(path, F_OK) == 0) { (void)!chdir(ups[i]); break; }
    }
    struct Expr* c = parse_expression("Get[\"src/internal/init.m\"]");
    if (c) { struct Expr* r = evaluate(c); expr_free(c); if (r) expr_free(r); }
}

/* Mute Mathilda's messages for the rest of the run, the way a Quiet around
 * everything would: for a test whose calls emit expected diagnostics. Stderr
 * itself stays open, so every FAIL line still shows -- pointing stderr at
 * /dev/null, as these tests used to, hid the failures along with the noise.
 * A Check[] inside still sees its messages, and they still count as fired.
 * Idempotent. */
void mth_msg_suppress_push(void);
static inline void test_mute_messages(void) {
    static int muted = 0;
    if (!muted) { mth_msg_suppress_push(); muted = 1; }
}

/* Point stderr (fd 2) at `path`, truncated, until test_stderr_restore; for a
 * test that reads back the messages a call printed. Returns the saved
 * descriptor, or -1 with stderr untouched. The swap is at the descriptor, so
 * the stderr FILE stays valid whatever happens -- and restoring by dup2, not
 * by reopening /dev/tty, keeps stderr working without a terminal: under CI or
 * a pipe that reopen fails, and every later FAIL line was silently lost. */
static inline int test_stderr_redirect(const char* path) {
    fflush(stderr);
    int saved = dup(STDERR_FILENO);
    if (saved < 0) return -1;
    int fd = open(path, O_RDWR | O_CREAT | O_TRUNC, 0644);
    if (fd < 0 || dup2(fd, STDERR_FILENO) < 0) {
        if (fd >= 0) close(fd);
        close(saved);
        return -1;
    }
    close(fd);
    return saved;
}

/* Undo test_stderr_redirect (a no-op for -1). */
static inline void test_stderr_restore(int saved) {
    if (saved < 0) return;
    fflush(stderr);
    dup2(saved, STDERR_FILENO);
    close(saved);
}

static inline void assert_eval_eq(const char* input, const char* expected, int is_fullform) {
    struct Expr* parsed = parse_expression(input);
    if (!parsed) { fprintf(stderr, "FAIL: does not parse: %s\n", input); abort(); }
    struct Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* str = is_fullform ? expr_to_string_fullform(evaluated) : expr_to_string(evaluated);
    if (strcmp(str, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected, str);
        abort();   /* not assert(): -DNDEBUG would make a failed check pass */
    }
    free(str);
    expr_free(evaluated);
}

/* Like assert_eval_eq but only requires that the printed result *starts
 * with* the expected prefix. Useful for arbitrary-precision outputs whose
 * low-order digits depend on exact bit-level rounding. */
static inline void assert_eval_startswith(const char* input, const char* prefix) {
    struct Expr* parsed = parse_expression(input);
    if (!parsed) { fprintf(stderr, "FAIL: does not parse: %s\n", input); abort(); }
    struct Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* str = expr_to_string(evaluated);
    size_t plen = strlen(prefix);
    int ok = (strncmp(str, prefix, plen) == 0);
    if (!ok) {
        fprintf(stderr, "FAIL: %s\n  Expected prefix: %s\n  Actual:          %s\n",
                input, prefix, str);
        abort();
    }
    free(str);
    expr_free(evaluated);
}

#define TEST(name) printf("Running test: %s\n", #name); name()

/* IMPORTANT: do NOT use the libc <assert.h> assert() here.  CMake
 * builds with CMAKE_BUILD_TYPE=Release pass -DNDEBUG, which silently
 * compiles every assert() to (void)0.  We need test failures to abort
 * the binary unconditionally so CTest sees the non-zero exit code. */
/* Normalise a result to an ORDINARY List so a test can walk data.function.args.
 *
 * A packed list is an EXPR_NDARRAY with no args[] array, and it is a List by
 * every language-level measure (`Head` is List, `ListQ` is True) -- so a C test
 * that asserts `type == EXPR_FUNCTION` is asserting a storage decision, not a
 * value. That decision moves: PACK_MIN_ELEMENTS went from 250 to 4 on
 * 2026-08-02 and five test files began failing on results that had simply
 * started packing, none of which had changed value. Consumes `r` and returns a
 * List the caller owns; a non-packed argument is returned unchanged. */
static inline Expr* test_delist(Expr* r) {
    if (!r || !is_packed_list(r)) return r;
    Expr* out = ndarray_to_nested_list(r);
    expr_free(r);
    return out;
}

#define ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "FAIL: assertion failed: %s\n    at %s:%d\n", \
                #cond, __FILE__, __LINE__); \
        exit(1); \
    } \
} while(0)

#define ASSERT_STR_EQ(a, b) do { \
    const char* _a_ = (a); \
    const char* _b_ = (b); \
    if (strcmp(_a_, _b_) != 0) { \
        fprintf(stderr, "FAIL: ASSERT_STR_EQ\n  expected: %s\n  actual:   %s\n" \
                        "    at %s:%d\n", \
                _b_, _a_, __FILE__, __LINE__); \
        exit(1); \
    } \
} while(0)

#define ASSERT_MSG(cond, fmt, ...) do { \
    if (!(cond)) { \
        fprintf(stderr, "Assertion failed: " fmt "\n", ##__VA_ARGS__); \
        fprintf(stderr, "    at %s:%d\n", __FILE__, __LINE__); \
        exit(1); \
    } \
} while(0)

#ifdef __GNUC__
__attribute__((constructor))
static void set_timeout() {
    /* Line-buffer stdout: under ctest stdout is a pipe (block-buffered), and a
     * suite killed by SIGALRM loses the buffer -- the log then cannot name the
     * test that was running (dsolve_tests on CI, 2026-11-10). */
    setvbuf(stdout, NULL, _IOLBF, 0);
    alarm(120);
}
#endif

