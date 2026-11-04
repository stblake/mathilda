/* test_dense_heads.h -- the shared chk helpers of the dense linear-algebra
 * head suites (test_cholesky.c, test_hessenberg.c, test_eigensystem.c,
 * test_singularvaluelist.c, test_minors_adjugate.c, test_orthogonalize.c).
 *
 * Failures are counted rather than asserted: the CMake Release build passes
 * -DNDEBUG, which would no-op assert().  Each suite defines `failures`
 * through this header and returns 1 from main when it is non-zero. */
#ifndef MATHILDA_TEST_DENSE_HEADS_H
#define MATHILDA_TEST_DENSE_HEADS_H

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

#if defined(__GNUC__)
#define MATHILDA_TEST_UNUSED __attribute__((unused))
#else
#define MATHILDA_TEST_UNUSED
#endif

static int failures = 0;

/* Evaluate `input`; fail unless it prints exactly `expected`. */
static void chk(const char* input, const char* expected) {
    struct Expr* p = parse_expression(input);
    if (!p) { fprintf(stderr, "FAIL(parse): %s\n", input); failures++; return; }
    struct Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected, s);
        failures++;
    }
    free(s);
    expr_free(e);
}

/* The InputForm of `input`. */
static void chk_if(const char* input, const char* expected) {
    char buf[8192];
    snprintf(buf, sizeof buf, "InputForm[%s]", input);
    chk(buf, expected);
}

static void chk_true(const char* input) { chk(input, "True"); }

/* `call` (a CompileDiagnostics[...] expression) reports Compiled -> True. */
static void chk_lowers(const char* call) {
    char buf[4096];
    snprintf(buf, sizeof buf, "(\"Compiled\" /. %s)", call);
    chk(buf, "True");
}

/* `call` reports Compiled -> False (a documented decline). */
MATHILDA_TEST_UNUSED static void chk_no_lower(const char* call) {
    char buf[4096];
    snprintf(buf, sizeof buf, "(\"Compiled\" /. %s)", call);
    chk(buf, "False");
}

/* `expr` emits message `msg` (Head::tag): Check returns $Failed. */
static void chk_msg(const char* expr, const char* msg) {
    char buf[8192];
    snprintf(buf, sizeof buf, "Check[%s, $Failed, %s]", expr, msg);
    chk(buf, "$Failed");
}

#endif /* MATHILDA_TEST_DENSE_HEADS_H */
