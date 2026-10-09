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

/* Constant (root-of-unity) phase angle-expansion: a circular trig head with an
 * affine argument k x + c Pi (c rational) reduces through the exact zero test. */
void test_tez_constant_phase(void) {
    assert_eval_eq("Simplify[Tan[Pi/2 - x] - Cot[x]]", "0", 0);
    assert_eval_eq("Simplify[Tan[x] + Tan[x + Pi/3] + Tan[x + 2 Pi/3] - 3 Tan[3 x]]", "0", 0);
    assert_eval_eq("Simplify[Sin[x] Sin[Pi/3 - x] Sin[Pi/3 + x] - 1/4 Sin[3 x]]", "0", 0);
    /* PossibleZeroQ must agree (previously returned a WRONG False). */
    assert_eval_eq("PossibleZeroQ[Tan[Pi/2 - x] - Cot[x]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Tan[x] + Tan[x + Pi/3] + Tan[x + 2 Pi/3] - 3 Tan[3 x]]", "True", 0);
    /* Soundness: genuine non-identities with a constant phase stay non-zero. */
    assert_eval_eq("PossibleZeroQ[Tan[x + Pi/3] - Tan[x]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Sin[x + Pi/5] - Cos[x]]", "False", 0);
}

/* Hyperbolic analogue: an affine IMAGINARY phase k x + i c Pi on a hyperbolic
 * head (Cosh[i c Pi] = Cos[c Pi], Sinh[i c Pi] = i Sin[c Pi]). */
void test_tez_hyperbolic_phase(void) {
    assert_eval_eq("Simplify[Tanh[I Pi/2 - x] + Coth[x]]", "0", 0);
    assert_eval_eq("Simplify[Tanh[x + I Pi] - Tanh[x]]", "0", 0);
    assert_eval_eq("Simplify[Tanh[x] + Tanh[x + I Pi/3] + Tanh[x + 2 I Pi/3] - 3 Tanh[3 x]]", "0", 0);
    /* Reciprocal-difference half-angle: Coth - Csch = Tanh[x/2]. */
    assert_eval_eq("Simplify[Tanh[x/2] - (Cosh[x] - 1)/Sinh[x]]", "0", 0);
    assert_eval_eq("Simplify[Coth[x] - Csch[x] - Tanh[x/2]]", "0", 0);
    /* Soundness: not identities. */
    assert_eval_eq("PossibleZeroQ[Tanh[x + I Pi/3] - Tanh[x]]", "False", 0);
    assert_eval_eq("Simplify[Coth[x] + Csch[x] - Tanh[x/2]] === 0", "False", 0);
}

/* Representative hyperbolic structural identities (Osborn's rule: note the +
 * signs in the Cosh addition and odd-power formulas). Guards the bulk of the
 * hyperbolic stress corpus. */
void test_tez_hyperbolic_structural(void) {
    assert_eval_eq("Simplify[Cosh[x]^2 - Sinh[x]^2 - 1]", "0", 0);
    assert_eval_eq("Simplify[Sech[x]^2 + Tanh[x]^2 - 1]", "0", 0);
    assert_eval_eq("Simplify[Coth[x]^2 - Csch[x]^2 - 1]", "0", 0);
    assert_eval_eq("Simplify[Sinh[2 x] - 2 Sinh[x] Cosh[x]]", "0", 0);
    assert_eval_eq("Simplify[Cosh[2 x] - (Cosh[x]^2 + Sinh[x]^2)]", "0", 0);
    assert_eval_eq("Simplify[Cosh[x + y] - (Cosh[x] Cosh[y] + Sinh[x] Sinh[y])]", "0", 0);
    assert_eval_eq("Simplify[Tanh[x + y] - (Tanh[x] + Tanh[y])/(1 + Tanh[x] Tanh[y])]", "0", 0);
    assert_eval_eq("Simplify[Sinh[3 x] - (3 Sinh[x] + 4 Sinh[x]^3)]", "0", 0);
    assert_eval_eq("Simplify[Cosh[3 x] - (4 Cosh[x]^3 - 3 Cosh[x])]", "0", 0);
    assert_eval_eq("Simplify[Sinh[x]^6 + Cosh[x]^6 - (15 Cosh[2 x] + Cosh[6 x])/16]", "0", 0);
    assert_eval_eq("Simplify[(1 + Tanh[x/2]^2)/(1 - Tanh[x/2]^2) - Cosh[x]]", "0", 0);
    assert_eval_eq("Simplify[Cosh[ArcSinh[x]] - Sqrt[1 + x^2]]", "0", 0);
    assert_eval_eq("Simplify[Sinh[2 ArcSinh[x]] - 2 x Sqrt[1 + x^2]]", "0", 0);
    assert_eval_eq("Simplify[Cosh[I x] - Cos[x]]", "0", 0);
    assert_eval_eq("Simplify[Sinh[I x] - I Sin[x]]", "0", 0);
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
    TEST(test_tez_constant_phase);
    TEST(test_tez_hyperbolic_phase);
    TEST(test_tez_hyperbolic_structural);

    printf("All trigexp_zero tests passed!\n");
    return 0;
}
