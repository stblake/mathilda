/* bronstein_test_utils.h -- shared helpers for the Bronstein integration
 * suites (tests/test_bronstein_*.c), which exercise src/internal/bronstein/
 * (and the Trager routines it generalises) function by function.
 *
 * br_init() starts a production-like kernel (init.m), loads the Bronstein
 * package (which loads Trager first) and puts Bronstein`Private` and
 * Trager`Private` on $ContextPath, so a suite calls internal functions by
 * their short names (kc[...], brCatch[...]).
 *
 * Every assertion is a HARD failure (exit 1 with the offending input): the
 * libc assert() used by assert_eval_eq is a no-op under -DNDEBUG.
 *
 * BR_ANTIDERIV(f, method) integrates f with Integrate[..., Method -> method]
 * and requires PossibleZeroQ[D[F, x] - f] (and no Integrate left in F).
 * PossibleZeroQ is the verifier (no numeric fallback): a false negative is a
 * PossibleZeroQ bug to fix there, not a reason to weaken the test.
 */
#ifndef BRONSTEIN_TEST_UTILS_H
#define BRONSTEIN_TEST_UTILS_H

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"
#include "print.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int br_assertions = 0;
static int br_failures = 0;   /* nonzero only under BR_KEEP_GOING */

/* Evaluate `input` and return its printed form (caller frees). */
static char* br_eval(const char* input) {
    Expr* parsed = parse_expression(input);
    if (!parsed) {
        fprintf(stderr, "FAIL (parse): %s\n", input);
        exit(1);
    }
    Expr* v = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string(v);
    expr_free(v);
    return s;
}

/* Hard-assert that `input` prints as `expected`. */
static void br_eq(const char* input, const char* expected, const char* file, int line) {
    char* got = br_eval(input);
    br_assertions++;
    if (strcmp(got, expected) != 0) {
        fprintf(stderr, "FAIL %s:%d\n  input:    %s\n  expected: %s\n  got:      %s\n",
                file, line, input, expected, got);
        free(got);
        /* BR_KEEP_GOING=1 reports every mismatch of a run (still exits 1). */
        if (getenv("BR_KEEP_GOING")) { br_failures++; return; }
        exit(1);
    }
    free(got);
}
#define BR_EQ(in, exp) br_eq((in), (exp), __FILE__, __LINE__)
#define BR_TRUE(in) br_eq((in), "True", __FILE__, __LINE__)
#define BR_FALSE(in) br_eq((in), "False", __FILE__, __LINE__)

/* Evaluate for side effects (definitions); the printed value is ignored. */
static void br_run(const char* input) {
    char* s = br_eval(input);
    free(s);
}

/* Assert that evaluating `input` throws the tragerFail of class `cls` and
 * reason `reason` (both strings): brFailTag[brCatch[input]] == {cls, reason}. */
static void br_throws(const char* input, const char* cls, const char* reason,
                      const char* file, int line) {
    char buf[8192], exp[512];
    snprintf(buf, sizeof buf, "brFailTag[brCatch[%s]]", input);
    snprintf(exp, sizeof exp, "{\"%s\", \"%s\"}", cls, reason);
    br_eq(buf, exp, file, line);
}
#define BR_THROWS(in, cls, reason) br_throws((in), (cls), (reason), __FILE__, __LINE__)

/* Integrate f by `method` and assert the result is a verified antiderivative. */
static void br_antideriv(const char* f, const char* method, const char* file, int line) {
    char buf[8192];
    snprintf(buf, sizeof buf,
             "Module[{brF = Integrate[%s, x, Method -> \"%s\"]},"
             " FreeQ[brF, Integrate] && TrueQ[PossibleZeroQ[D[brF, x] - (%s)]]]",
             f, method, f);
    br_eq(buf, "True", file, line);
}
#define BR_ANTIDERIV(f, method) br_antideriv((f), (method), __FILE__, __LINE__)

/* Start the kernel and load the Bronstein package. */
static void br_init(void) {
    symtab_init();
    core_init();
    test_load_init_m();          /* chdir to the repo root; LoadModule resolves */
    test_mute_messages();        /* expected diagnostics stay quiet; FAILs still print */
    br_run("LoadModule[\"bronstein/Bronstein.m\"]");
    br_run("$ContextPath = Join[$ContextPath, {\"Bronstein`Private`\", \"Trager`Private`\"}]; True");
}

#define BR_DONE(name) do { \
    if (br_failures) { fprintf(stderr, "%s: %d of %d assertions FAILED\n", (name), \
                               br_failures, br_assertions); exit(1); } \
    printf("%s: all %d assertions passed\n", (name), br_assertions); } while (0)

#endif /* BRONSTEIN_TEST_UTILS_H */
