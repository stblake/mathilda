/*
 * test_possiblezeroq_expcombine_stress.c
 *
 * Stress corpus for the exponential-combining normalisation (Stage 0.5 of
 * src/zero_test.c) added for POSSIBLE_ZEROQ_IMPROVEMENTS.md #1.
 *
 * The failure mode: a Gaussian x Erf/Erfi residual whose same-base
 * exponentials (E^(-g) from the solution, E^(+g) from differentiating the
 * error function) sit in SEPARATE summands.  The Times like-base collector
 * never combines them to E^0 = 1, so every Schwartz-Zippel sample numericalises
 * as a tiny*huge catastrophic cancellation and the precision ladder climbs to
 * 1000 bits per sample -- effectively a hang (> 20 s, flooding $IterationLimit).
 * The recogniser now ExpandAlls an input carrying a NON-LINEAR-exponent
 * exponential first, distributing the sums so the exponentials become adjacent
 * factors and collapse before the ladder ever runs.
 *
 * Two guards, both hard (exit(1), so they gate ctest even under -DNDEBUG):
 *   - correctness: each True identity decides True, each matched non-identity
 *     decides False;
 *   - anti-hang: each True identity is decided within a generous per-case time
 *     budget (a re-introduced hang trips this, and the 120 s test_utils.h alarm
 *     is the backstop that turns a true infinite hang into a ctest failure).
 *
 * The "True" cases are ODE integrating-factor residuals and error-function
 * derivative identities; the "False" cases add a clean, non-overflowing term to
 * a zero residual so the sampler sees a genuine non-zero.  Every case here was
 * confirmed against the running binary before being asserted.
 *
 * KNOWN LIMITATION (documented, NOT asserted): a genuine non-identity that
 * differs from a zero residual only INSIDE the E*Erf[imaginary] tiny*huge terms
 * (e.g. D[b,x]-b for the body b below) is mis-classified True.  That is a
 * pre-existing Schwartz-Zippel sampler weakness (the imaginary-argument Erf
 * overflows to Inf at the sampler's moderate range, degrading to UNKNOWN->True),
 * independent of this normalisation, which only ever transforms the input to a
 * value-equal form.  So are the flat-product false positives such as
 * PossibleZeroQ[D[Erf[x],x] - 3/Sqrt[Pi] E^(-x^2)].
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "parse.h"
#include "symtab.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Hard-asserting PZQ verdict check (exit(1) on mismatch; NDEBUG-proof). */
static void pzq(const char* input, const char* expected) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    Expr* evald = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string(evald);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n", input, expected, s);
        free(s); expr_free(evald); exit(1);
    }
    free(s);
    expr_free(evald);
}

/* Verdict AND wall-clock budget: the anti-hang guard (the bare form of the
 * flagship case took > 20 s). */
static void pzq_timed(const char* input, const char* expected, double budget_s) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    clock_t t0 = clock();
    Expr* evald = evaluate(parsed);
    double dt = (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
    expr_free(parsed);
    char* s = expr_to_string(evald);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n", input, expected, s);
        free(s); expr_free(evald); exit(1);
    }
    if (dt > budget_s) {
        fprintf(stderr, "FAIL (too slow): %s\n  verdict %s ok, but took %.2fs (budget %.2fs)\n",
                input, s, dt, budget_s);
        free(s); expr_free(evald); exit(1);
    }
    free(s);
    expr_free(evald);
}

/* The body b that solves y'' + x y' + y == 0. */
#define BODY(c1, c2) \
    "E^(-1/2 x^2) (" c1 " + I (" c2 " Sqrt[Pi] Erf[-(I x)/Sqrt[2]]) / Sqrt[2])"
#define RESID(c1, c2) \
    "With[{b = " BODY(c1, c2) "}, D[b, {x, 2}] + x D[b, x] + b]"

/* --- True: y'' + x y' + y == 0 integrating-factor residuals --- */
static void test_ode_residuals_true(void) {
    pzq_timed("PossibleZeroQ[" RESID("C[1]", "C[2]") "]", "True", 5.0);
    pzq_timed("PossibleZeroQ[" RESID("3", "5") "]", "True", 5.0);
    pzq_timed("PossibleZeroQ[" RESID("C[1]", "0") "]", "True", 5.0);
    pzq_timed("PossibleZeroQ[" RESID("0", "C[2]") "]", "True", 5.0);
    pzq_timed("PossibleZeroQ[" RESID("1", "1") "]", "True", 5.0);
    pzq_timed("PossibleZeroQ[" RESID("-2", "7") "]", "True", 5.0);
}

/* --- True: Erf / Erfi derivative identities (Gaussian appears explicitly) --- */
static void test_erf_derivative_identities_true(void) {
    pzq_timed("PossibleZeroQ[D[Erf[x], x] - 2/Sqrt[Pi] E^(-x^2)]", "True", 5.0);
    pzq_timed("PossibleZeroQ[D[Erf[2 x], x] - 4/Sqrt[Pi] E^(-4 x^2)]", "True", 5.0);
    pzq_timed("PossibleZeroQ[D[Erf[3 x], x] - 6/Sqrt[Pi] E^(-9 x^2)]", "True", 5.0);
    pzq_timed("PossibleZeroQ[D[Erf[x], {x, 2}] + 4 x/Sqrt[Pi] E^(-x^2)]", "True", 5.0);
    pzq_timed("PossibleZeroQ[D[Erfi[x], x] - 2/Sqrt[Pi] E^(x^2)]", "True", 5.0);
    pzq_timed("PossibleZeroQ[D[Erfi[2 x], x] - 4/Sqrt[Pi] E^(4 x^2)]", "True", 5.0);
}

/* --- True: Gaussian / product derivative identities --- */
static void test_gaussian_product_identities_true(void) {
    pzq_timed("PossibleZeroQ[D[E^(-x^2), x] + 2 x E^(-x^2)]", "True", 5.0);
    pzq_timed("PossibleZeroQ[D[E^(-x^2/2), x] + x E^(-x^2/2)]", "True", 5.0);
    pzq_timed("PossibleZeroQ[D[x E^(-x^2), x] - (E^(-x^2) - 2 x^2 E^(-x^2))]", "True", 5.0);
    pzq_timed("PossibleZeroQ[D[E^(-x^2) Erf[x], x] "
              "- (2/Sqrt[Pi] E^(-2 x^2) - 2 x E^(-x^2) Erf[x])]", "True", 5.0);
    pzq_timed("PossibleZeroQ[D[E^(x^2) Erfi[x], x] "
              "- (2 x E^(x^2) Erfi[x] + 2/Sqrt[Pi] E^(2 x^2))]", "True", 5.0);
}

/* --- False: a clean, non-overflowing term added to a zero residual --- */
static void test_matched_nonidentities_false(void) {
    pzq("PossibleZeroQ[" RESID("C[1]", "C[2]") " + E^(-x^2/2)]", "False");
    pzq("PossibleZeroQ[" RESID("C[1]", "C[2]") " + 1]", "False");
    pzq("PossibleZeroQ[" RESID("C[1]", "C[2]") " + x]", "False");
    pzq("PossibleZeroQ[" RESID("C[1]", "C[2]") " + x^2]", "False");
    pzq("PossibleZeroQ[" RESID("C[1]", "C[2]") " + Cos[x]]", "False");
    pzq("PossibleZeroQ[D[Erf[x], x] - 2/Sqrt[Pi] E^(-x^2) + 1]", "False");
    pzq("PossibleZeroQ[D[Erf[x], x] - 2/Sqrt[Pi] E^(-x^2) + E^(-x^2/2)]", "False");
    pzq("PossibleZeroQ[D[Erfi[x], x] - 2/Sqrt[Pi] E^(x^2) + Cos[x]]", "False");
    pzq("PossibleZeroQ[E^(-x^2) Erf[x] - E^(-x^2)]", "False");
}

/* --- Verdict preservation: ExpandAll is value-preserving --- */
static void test_verdict_preservation(void) {
    pzq("PossibleZeroQ[E^(2 x) - E^x E^x]", "True");     /* affine exp: gate skips */
    pzq("PossibleZeroQ[E^(x + y) - E^x E^y]", "True");   /* affine exp: gate skips */
    pzq("PossibleZeroQ[E^x - E^(2 x)]", "False");
    pzq("PossibleZeroQ[E^(-x^2) E^(x^2) - 1]", "True");  /* non-linear: gate fires */
    pzq("PossibleZeroQ[E^(I x) + E^(-I x) - 2 Cos[x]]", "True");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_ode_residuals_true);
    TEST(test_erf_derivative_identities_true);
    TEST(test_gaussian_product_identities_true);
    TEST(test_matched_nonidentities_false);
    TEST(test_verdict_preservation);

    printf("\nAll PossibleZeroQ exp-combining stress tests passed. (31 cases)\n");
    return 0;
}
