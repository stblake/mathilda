/*
 * test_trigexp_zero.c -- exact trig/exp-kernel zero test (FP2) and the
 * trig-log canonicalization (FP1).
 *
 * FP2 (src/simp/simp_trigexp_zero.c) proves a rational function of a single
 * kernel t = E^(i x) identically zero by exact rational point-evaluation on a
 * Nullstellensatz grid — closing the Sec^n/Csc^n and symbolic-parameter Risch
 * diff-back identities that the general Simplify search cannot reduce / hangs
 * on (SIMPLIFY_GAPS.md Families 1 & 3). Reached from both Simplify (top-level
 * fast path) and PossibleZeroQ (a symbolic zero_test.c stage).
 *
 * FP1 (src/simp/simp_trig_pi.c) normalizes Log[Sec[u]^2] -> -Log[Cos[u]^2] etc.
 * so the log-fusion pass can cancel it (SIMPLIFY_GAPS.md Family 2 / D2).
 */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "simp_trigexp_zero.h"
#include "parse.h"
#include "eval.h"
#include <time.h>

#define TEST(func) do { printf("Running test: %s\n", #func); fflush(stdout); func(); } while(0)

static double seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

/* Run the C primitive directly on a parsed+evaluated expression. */
static TrigExpZeroResult tez(const char* src) {
    Expr* e = evaluate(parse_expression(src));
    TrigExpZeroResult r = trigexp_rational_is_zero(e);
    expr_free(e);
    return r;
}

/* ---- FP2: the C primitive verdicts ---- */

void test_tez_multiple_angle_identities(void) {
    /* Multiple-angle identities that are exactly zero. */
    ASSERT(tez("Sin[3 x] - 3 Sin[x] + 4 Sin[x]^3") == TRIGEXP_ZERO_TRUE);
    ASSERT(tez("Cos[3 x] - 4 Cos[x]^3 + 3 Cos[x]") == TRIGEXP_ZERO_TRUE);
    ASSERT(tez("Sin[2 x] - 2 Sin[x] Cos[x]") == TRIGEXP_ZERO_TRUE);
    ASSERT(tez("Cos[2 x] - 1 + 2 Sin[x]^2") == TRIGEXP_ZERO_TRUE);
    ASSERT(tez("Sin[x]^2 + Cos[x]^2 - 1") == TRIGEXP_ZERO_TRUE);
}

void test_tez_genuine_nonzero(void) {
    /* Not identities: the primitive must NOT claim zero. */
    ASSERT(tez("Sin[3 x] - 3 Sin[x]") == TRIGEXP_ZERO_FALSE);
    ASSERT(tez("Cos[2 x] - Cos[x]") == TRIGEXP_ZERO_FALSE);
    ASSERT(tez("Sec[x]^2 - Tan[x]") == TRIGEXP_ZERO_FALSE);
}

void test_tez_declines_gracefully(void) {
    /* Non-single-kernel / non-rational forms: decline (UNKNOWN), never crash. */
    ASSERT(tez("Tan[Log[x]]/x") == TRIGEXP_ZERO_UNKNOWN);   /* nested kernel */
    ASSERT(tez("Sin[x] + Sin[y]") == TRIGEXP_ZERO_UNKNOWN); /* two kernel vars */
    ASSERT(tez("x + 1") == TRIGEXP_ZERO_UNKNOWN);           /* no kernel */
}

void test_tez_symbolic_parameters(void) {
    /* Family 3: an identity with symbolic parameters a, b decided exactly.
     * Cos[2x] = 1 - 2 Sin[x]^2 scaled by an arbitrary rational parameter. */
    ASSERT(tez("a (Cos[2 x] - 1 + 2 Sin[x]^2)") == TRIGEXP_ZERO_TRUE);
    ASSERT(tez("a Sin[2 x] - 2 a Sin[x] Cos[x] + b - b") == TRIGEXP_ZERO_TRUE);
    /* A symbolic-parameter non-identity must not read as zero. */
    ASSERT(tez("a Sin[2 x] - Sin[x] Cos[x]") == TRIGEXP_ZERO_FALSE);
}

/* ---- FP2 through Simplify and PossibleZeroQ ---- */

void test_tez_secant_diffback_simplify(void) {
    /* The multiple-angle Sec^3 antiderivative diff-back (SIMPLIFY_GAPS.md
     * Family 1): identically 0 but the general Simplify search runs >40 s
     * without terminating. The fast path must reduce it to 0 quickly. */
    const char* src =
        "Simplify[D[(8 Sin[x] + 4 Sin[5 x] + 12 Sin[3 x] "
        "+ (-10 - 15 Cos[2 x] - 6 Cos[4 x] - Cos[6 x]) Log[2 - 2 Sin[x]] "
        "+ (10 + Cos[6 x] + 6 Cos[4 x] + 15 Cos[2 x]) Log[2 + 2 Sin[x]]) "
        "/ (40 + 4 Cos[6 x] + 24 Cos[4 x] + 60 Cos[2 x]), x] - Sec[x]^3]";
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    double t0 = seconds();
    Expr* r = evaluate(parsed);
    double dt = seconds() - t0;
    expr_free(parsed);
    char* s = expr_to_string(r);
    ASSERT_STR_EQ(s, "0");
    ASSERT_MSG(dt < 5.0, "Sec^3 diff-back must certify in well under 5 s");
    free(s);
    expr_free(r);
}

void test_tez_possiblezeroq(void) {
    /* PossibleZeroQ decides the same identities exactly (symbolic stage). */
    assert_eval_eq("PossibleZeroQ[Sin[3 x] - 3 Sin[x] + 4 Sin[x]^3]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sin[3 x] - 3 Sin[x]]", "False", 0);
    /* No regression on pure-rational PossibleZeroQ. */
    assert_eval_eq("PossibleZeroQ[(x + 1)^2 - x^2 - 2 x - 1]", "True", 0);
    assert_eval_eq("PossibleZeroQ[x^2 - 1]", "False", 0);
}

void test_tez_small_cases_fast(void) {
    /* Typical cases must resolve in milliseconds, not seconds. */
    double t0 = seconds();
    Expr* r = evaluate(parse_expression("Simplify[Sin[3 x] - 3 Sin[x] + 4 Sin[x]^3]"));
    double dt = seconds() - t0;
    char* s = expr_to_string(r);
    ASSERT_STR_EQ(s, "0");
    ASSERT_MSG(dt < 0.5, "small multiple-angle identity must be fast");
    free(s);
    expr_free(r);
}

/* ---- FP1: trig-log canonicalization ---- */

void test_fp1_log_reciprocal_squared(void) {
    assert_eval_eq("Simplify[Log[Sec[x]^2] + Log[Cos[x]^2]]", "0", 0);
    assert_eval_eq("Simplify[Log[Csc[x]^2] + Log[Sin[x]^2]]", "0", 0);
    assert_eval_eq("Simplify[1/2 Log[1 + Tan[x]^2] + 1/2 Log[Cos[x]^2]]", "0", 0);
    assert_eval_eq("Simplify[1/2 Log[1 + Cot[x]^2] + 1/2 Log[Sin[x]^2]]", "0", 0);
}

void test_fp1_no_regression(void) {
    /* Baselines the pythag rules already handle must be preserved. */
    assert_eval_eq("Simplify[1 + Tan[x]^2]", "Sec[x]^2", 0);
    assert_eval_eq("Simplify[Sin[x]^2 + Cos[x]^2]", "1", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_tez_multiple_angle_identities);
    TEST(test_tez_genuine_nonzero);
    TEST(test_tez_declines_gracefully);
    TEST(test_tez_symbolic_parameters);
    TEST(test_tez_secant_diffback_simplify);
    TEST(test_tez_possiblezeroq);
    TEST(test_tez_small_cases_fast);
    TEST(test_fp1_log_reciprocal_squared);
    TEST(test_fp1_no_regression);

    printf("All trigexp_zero tests passed!\n");
    return 0;
}
