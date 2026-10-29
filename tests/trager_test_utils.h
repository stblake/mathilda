/* trager_test_utils.h -- shared helpers for the Trager integration suites
 * (tests/test_trager_*.c), which exercise src/internal/trager/ function by
 * function.
 *
 * tr_init() starts a production-like kernel (init.m), loads the package and
 * puts Trager`Private` on $ContextPath, so a suite can call the package's
 * internal functions by their short names (canonicalForm[...], hnf[...]).
 *
 * Every assertion is a HARD failure (exit 1 with the offending input): the
 * libc assert() used by assert_eval_eq is a no-op under -DNDEBUG.
 *
 * Correctness of an antiderivative is asserted with PossibleZeroQ:
 * TR_ANTIDERIV(f, method) integrates f with Integrate[..., Method -> method]
 * and requires PossibleZeroQ[D[F, x] - f] (and no Integrate left in F).
 * PossibleZeroQ has false negatives on nested radicals (e.g. the answer for
 * ((2x+1)/(x+3))^(1/6), zero to 40 digits), so a False is re-checked by the
 * package's numericZeroQ: |D[F, x] - f| < 10^-30 at >= 4 rational points with
 * 60-digit arithmetic -- the same backstop the pipeline itself uses.
 */
#ifndef TRAGER_TEST_UTILS_H
#define TRAGER_TEST_UTILS_H

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

static int tr_assertions = 0;

/* Evaluate `input` and return its printed form (caller frees). */
static char* tr_eval(const char* input) {
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
static void tr_eq(const char* input, const char* expected, const char* file, int line) {
    char* got = tr_eval(input);
    tr_assertions++;
    if (strcmp(got, expected) != 0) {
        fprintf(stderr, "FAIL %s:%d\n  input:    %s\n  expected: %s\n  got:      %s\n",
                file, line, input, expected, got);
        free(got);
        exit(1);
    }
    free(got);
}
#define TR_EQ(in, exp) tr_eq((in), (exp), __FILE__, __LINE__)
#define TR_TRUE(in) tr_eq((in), "True", __FILE__, __LINE__)
#define TR_FALSE(in) tr_eq((in), "False", __FILE__, __LINE__)

/* Integrate f by `method` and assert the result is a verified antiderivative:
 * Integrate is gone and PossibleZeroQ[D[F, x] - f] is True. */
static void tr_antideriv(const char* f, const char* method, const char* file, int line) {
    char buf[8192];
    snprintf(buf, sizeof buf,
             "Module[{trF = Integrate[%s, x, Method -> \"%s\"], trD},"
             " trD = D[trF, x] - (%s);"
             " FreeQ[trF, Integrate] && (TrueQ[PossibleZeroQ[trD]] ||"
             " Trager`Private`numericZeroQ[trD, x])]",
             f, method, f);
    tr_eq(buf, "True", file, line);
}
#define TR_ANTIDERIV(f, method) tr_antideriv((f), (method), __FILE__, __LINE__)

/* Start the kernel and load the Trager package. */
static void tr_init(void) {
    symtab_init();
    core_init();
    test_load_init_m();          /* chdir to the repo root; LoadModule resolves */
    test_mute_messages();        /* expected diagnostics stay quiet; FAILs still print */
    char* r = tr_eval("LoadModule[\"trager/Trager.m\"]");
    free(r);
    r = tr_eval("$ContextPath = Append[$ContextPath, \"Trager`Private`\"]; True");
    free(r);
}

#endif /* TRAGER_TEST_UTILS_H */
