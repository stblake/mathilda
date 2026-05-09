

#pragma once
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <stdlib.h>

extern struct Expr* parse_expression(const char*);
extern struct Expr* evaluate(struct Expr*);
extern char* expr_to_string(struct Expr*);
extern char* expr_to_string_fullform(struct Expr*);
extern void expr_free(struct Expr*);

static inline void assert_eval_eq(const char* input, const char* expected, int is_fullform) {
    struct Expr* parsed = parse_expression(input);
    assert(parsed != NULL);
    struct Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* str = is_fullform ? expr_to_string_fullform(evaluated) : expr_to_string(evaluated);
    if (strcmp(str, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected, str);
    }
    assert(strcmp(str, expected) == 0);
    free(str);
    expr_free(evaluated);
}

/* Like assert_eval_eq but only requires that the printed result *starts
 * with* the expected prefix. Useful for arbitrary-precision outputs whose
 * low-order digits depend on exact bit-level rounding. */
static inline void assert_eval_startswith(const char* input, const char* prefix) {
    struct Expr* parsed = parse_expression(input);
    assert(parsed != NULL);
    struct Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* str = expr_to_string(evaluated);
    size_t plen = strlen(prefix);
    int ok = (strncmp(str, prefix, plen) == 0);
    if (!ok) {
        fprintf(stderr, "FAIL: %s\n  Expected prefix: %s\n  Actual:          %s\n",
                input, prefix, str);
    }
    assert(ok);
    free(str);
    expr_free(evaluated);
}

#define TEST(name) printf("Running test: %s\n", #name); name()

/* IMPORTANT: do NOT use the libc <assert.h> assert() here.  CMake
 * builds with CMAKE_BUILD_TYPE=Release pass -DNDEBUG, which silently
 * compiles every assert() to (void)0.  We need test failures to abort
 * the binary unconditionally so CTest sees the non-zero exit code. */
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
    alarm(60);
}
#endif

