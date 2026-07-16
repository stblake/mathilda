/* test_risch_hypertangent.c — the hypertangent case (Bronstein §5.10).
 *
 * Verifies:
 *   - Risch`PolynomialReduce (§5.4): reducing a polynomial modulo derivatives in
 *     a nonlinear monomial to a remainder of degree < deg_t(Dt);
 *   - Risch`IntegrateHypertangentPolynomial (§5.10): the polynomial-part tangent
 *     integrator, returning {q, c} with p - D[q] - c D(t^2+1)/(t^2+1) in k.
 *
 * The monomial t = tan(x) has Dt = 1 + t^2 (a = 1); t = tan(2x) has Dt = 2(1+t^2)
 * (a = 2).  Worked example 5.10.1: integral of (tan^2 x + x tan x + 1).
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    }
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

static void run_zero(const char* diff) {
    char buf[2048];
    snprintf(buf, sizeof buf, "Expand[Together[%s]]", diff);
    run_test(buf, "0");
}

/* t = tan(x): Dt = 1 + t^2. */
#define HT "{x -> 1, t -> 1 + t^2}"

/* ---- PolynomialReduce (§5.4) ------------------------------------------ */
static void test_polynomial_reduce(void) {
    /* Example 5.10.1: reduce t^2 + x t + 1 -> q = t, r = x t. */
    run_zero("Risch`PolynomialReduce[t^2 + x t + 1, t, " HT "][[1]] - t");
    run_zero("Risch`PolynomialReduce[t^2 + x t + 1, t, " HT "][[2]] - x t");
    /* Reconstruction p = D[q] + r for a higher-degree p. */
    run_zero("Risch`Derivation[Risch`PolynomialReduce[t^4 + t + x, t, " HT "][[1]], " HT "] "
             "+ Risch`PolynomialReduce[t^4 + t + x, t, " HT "][[2]] - (t^4 + t + x)");
    /* The remainder has degree < deg_t(Dt) = 2. */
    run_zero("Coefficient[Risch`PolynomialReduce[t^4 + t + x, t, " HT "][[2]], t, 2]");
    run_zero("Coefficient[Risch`PolynomialReduce[t^4 + t + x, t, " HT "][[2]], t, 3]");
}

/* ---- IntegrateHypertangentPolynomial (§5.10) -------------------------- */
static void test_integrate_hypertangent_poly(void) {
#define IHP "Risch`IntegrateHypertangentPolynomial"
    /* Example 5.10.1: {q, c} = {t, x/2}. */
    run_zero(IHP "[t^2 + x t + 1, t, " HT "][[1]] - t");
    run_zero(IHP "[t^2 + x t + 1, t, " HT "][[2]] - x/2");
    /* D[c] = 1/2 != 0, so (the x tan x part of) the integral is not elementary. */
    run_test("D[" IHP "[t^2 + x t + 1, t, " HT "][[2]], x]", "Rational[1, 2]");
    /* The §5.10 certificate: p - D[q] - c D[t^2+1]/(t^2+1) is in k (here 0). */
    run_zero("(t^2 + x t + 1) - Risch`Derivation[" IHP "[t^2 + x t + 1, t, " HT "][[1]], " HT "] "
             "- " IHP "[t^2 + x t + 1, t, " HT "][[2]] Risch`Derivation[t^2 + 1, " HT "]/(t^2 + 1)");

    /* Integral of tan(x): {0, 1/2}; D[c] = 0, so it IS elementary
     * (∫ tan x = (1/2) Log[tan^2 x + 1] = -Log Cos x). */
    run_test(IHP "[t, t, " HT "]", "List[0, Rational[1, 2]]");
    run_test("D[" IHP "[t, t, " HT "][[2]], x]", "0");

    /* tan(2x): Dt = 2(1 + t^2), a = 2, so c = 1/(2 a) = 1/4. */
    run_test(IHP "[t, t, {x -> 1, t -> 2 (1 + t^2)}]", "List[0, Rational[1, 4]]");
#undef IHP
}

/* ---- The §5.10 certificate p - D[q] - c D(t^2+1)/(t^2+1) in k --------- */
/* IntegrateHypertangentPolynomial's output {q, c} satisfies (*): the residual
 * is an element of k, i.e. FREE OF t (it need NOT be zero — a nonzero constant
 * is the leftover base-field integrand ∫(element of k)). */
static void run_certificate(const char* p, const char* deriv) {
    char buf[4096];
    snprintf(buf, sizeof buf,
        "FreeQ[Expand[Together["
        "(%s) - Risch`Derivation[Risch`IntegrateHypertangentPolynomial[%s, t, %s][[1]], %s] "
        "- Risch`IntegrateHypertangentPolynomial[%s, t, %s][[2]] "
        "Risch`Derivation[t^2 + 1, %s]/(t^2 + 1)]], t]",
        p, p, deriv, deriv, p, deriv, deriv);
    run_test(buf, "True");
}

/* ---- Higher-degree p and the Dc = 0 / Dc != 0 boundary ---------------- */
static void test_hypertangent_higher_degree(void) {
    /* The §5.10 certificate holds for a range of pure-t polynomials. */
    run_certificate("t^3", HT);
    run_certificate("t^4", HT);
    run_certificate("t^5", HT);
    run_certificate("t^6", HT);
    run_certificate("t^4 + t^2 + 1", HT);
    run_certificate("t^3 + 2 t^2 + t", HT);

    /* Concrete values: ∫ t^3 -> {t^2/2, -1/2} (c constant: Dc = 0, elementary). */
    run_zero("Risch`IntegrateHypertangentPolynomial[t^3, t, " HT "][[1]] - t^2/2");
    run_zero("Risch`IntegrateHypertangentPolynomial[t^3, t, " HT "][[2]] + 1/2");
    run_test("D[Risch`IntegrateHypertangentPolynomial[t^3, t, " HT "][[2]], x]", "0");
    /* ∫ t^4 -> {-t + t^3/3, 0}. */
    run_zero("Risch`IntegrateHypertangentPolynomial[t^4, t, " HT "][[1]] - (-t + t^3/3)");
    run_test("Risch`IntegrateHypertangentPolynomial[t^4, t, " HT "][[2]]", "0");

    /* Elementary: ∫ x t^2 -> {t x, -1/2}; the log coefficient is constant, Dc = 0. */
    run_test("D[Risch`IntegrateHypertangentPolynomial[x t^2, t, " HT "][[2]], x]", "0");
    /* Non-elementary boundary: the LINEAR-in-t term x t forces c = (x-1)/2, so
     * Dc = 1/2 != 0 (the coefficient of the log term must be constant for an
     * elementary integral — Liouville). */
    run_test("D[Risch`IntegrateHypertangentPolynomial[x t^2 + x t, t, " HT "][[2]], x]", "Rational[1, 2]");
    run_certificate("x t^2 + t^2 + x t + 1", HT);
}

/* ---- Scaling a: tan(a x), Dt = a (1 + t^2) ---------------------------- */
static void test_hypertangent_scaling(void) {
    /* ∫ t with Dt = a(1+t^2) -> c = 1/(2a): a=1 -> 1/2, a=2 -> 1/4, a=3 -> 1/6. */
    run_test("Risch`IntegrateHypertangentPolynomial[t, t, {x -> 1, t -> 1 + t^2}]",
             "List[0, Rational[1, 2]]");
    run_test("Risch`IntegrateHypertangentPolynomial[t, t, {x -> 1, t -> 2 (1 + t^2)}]",
             "List[0, Rational[1, 4]]");
    run_test("Risch`IntegrateHypertangentPolynomial[t, t, {x -> 1, t -> 3 (1 + t^2)}]",
             "List[0, Rational[1, 6]]");
    /* Certificate still holds for a = 3. */
    run_certificate("t^3 + t", "{x -> 1, t -> 3 (1 + t^2)}");
}

/* ---- Robustness: non-hypertangent monomial + malformed input ---------- */
static void test_hypertangent_robustness(void) {
    /* Exponential monomial is NOT hypertangent (Dt = t, not a(t^2+1)). */
    run_test("Head[Risch`IntegrateHypertangentPolynomial[t^2, t, {x -> 1, t -> t}]] "
             "=== Risch`IntegrateHypertangentPolynomial", "True");
    /* Log monomial is not hypertangent either. */
    run_test("Head[Risch`IntegrateHypertangentPolynomial[t, t, {x -> 1, t -> 1/x}]] "
             "=== Risch`IntegrateHypertangentPolynomial", "True");
    /* Wrong arity / malformed derivation. */
    run_test("Head[Risch`IntegrateHypertangentPolynomial[t, t]] "
             "=== Risch`IntegrateHypertangentPolynomial", "True");
    run_test("Head[Risch`IntegrateHypertangentPolynomial[t, t, {x, t}]] "
             "=== Risch`IntegrateHypertangentPolynomial", "True");
}

/* t = tan(2x): Dt = 2(1+t^2);  t = tan(x/2): Dt = (1+t^2)/2. */
#define HT2 "{x -> 1, t -> 2 (1 + t^2)}"
#define HALF "{x -> 1, t -> (1 + t^2)/2}"
#define RR "Risch`ResidueReduce"
#define IHT "Risch`IntegrateHypertangent"

/* ---- ResidueReduce (§5.6, the residue criterion) --------------------- */
static void test_residue_reduce(void) {
    /* Simple normal pole with a constant Rothstein-Trager residue:
     * h = 1/(t-2), residue 1/(1+2^2) = 1/5, log arg t-2. */
    run_zero(RR "[1/(t - 2), t, " HT "][[1]] - 1/5 Log[t - 2]");
    run_test(RR "[1/(t - 2), t, " HT "][[2]]", "True");
    /* No simple poles (h is a polynomial in t) -> {0, True}. */
    run_test(RR "[t^2 + 1, t, " HT "]", "List[0, True]");
    run_test(RR "[t^3 + x, t, " HT "]", "List[0, True]");
    run_test(RR "[0, t, " HT "]", "List[0, True]");
    /* Non-constant residue (depends on x) -> {0, False}: h = 1/(t-x) has residue
     * 1/x^2, so it has no elementary integral over k(t) (Thm 5.6.1). */
    run_test(RR "[1/(t - x), t, " HT "][[2]]", "False");
    run_zero(RR "[1/(t - x), t, " HT "][[1]]");
}

/* h - D[g2] must lie in k[t] (its denominator is free of t) when beta = True. */
static void run_rr_reduces(const char* h, const char* deriv) {
    char buf[4096];
    snprintf(buf, sizeof buf,
        "With[{s = " RR "[%s, t, %s]}, s[[2]] === True && "
        "FreeQ[Denominator[Together[(%s) - Risch`Derivation[s[[1]], %s]]], t]]",
        h, deriv, h, deriv);
    run_test(buf, "True");
}

static void test_residue_reduce_diffback(void) {
    run_rr_reduces("1/(t - 2)", HT);
    run_rr_reduces("3/(t - 2)", HT);
    run_rr_reduces("1/((t - 2)(t - 4))", HT);
    run_rr_reduces("1/(t - 2) + 1/(t + 3)", HT);
    run_rr_reduces("(2 t + 1)/((t - 1)(t - 5))", HT);
    run_rr_reduces("1/(t - 2)", HT2);          /* scaled monomial */
}

/* ---- IntegrateHypertangent (§5.10, the full driver): book examples --- */
static void test_iht_book(void) {
    /* Example 5.10.1: ∫(tan^2 + x tan + 1) = tan(x) + ∫ x tan(x) (non-elem),
     * so g = t and beta = False. */
    run_test(IHT "[t^2 + x t + 1, t, " HT "]", "List[t, False]");
    /* Example 5.10.2: ∫ sin(x)/x, t = tan(x/2): g = 0, beta = False. */
    run_test(IHT "[2 t/(x (t^2 + 1)), t, " HALF "]", "List[0, False]");
    /* Pure log: ∫ D[Log(tan-2)] -> {Log[t-2], True}. */
    run_zero(IHT "[(1 + t^2)/(t - 2), t, " HT "][[1]] - Log[t - 2]");
    run_test(IHT "[(1 + t^2)/(t - 2), t, " HT "][[2]]", "True");
    /* ∫ x tan(x) is not elementary: g = 0, beta = False. */
    run_test(IHT "[x t, t, " HT "]", "List[0, False]");
}

/* When beta = True, f - D[g] must lie in k (free of t). */
static void run_iht_reduces(const char* f, const char* deriv) {
    char buf[4096];
    snprintf(buf, sizeof buf,
        "With[{s = " IHT "[%s, t, %s]}, s[[2]] === True && "
        "FreeQ[Together[(%s) - Risch`Derivation[s[[1]], %s]], t]]",
        f, deriv, f, deriv);
    run_test(buf, "True");
}

/* Constructed: f = D[g0] + base with g0 in k(t) and base in k (free of t), so f
 * is elementary and reduces to the base field — the driver must return beta =
 * True with f - D[g] free of t. */
static void run_iht_constructed(const char* g0, const char* base, const char* deriv) {
    char buf[6000];
    snprintf(buf, sizeof buf,
        "With[{f = Risch`Derivation[%s, %s] + (%s)}, "
        "With[{s = " IHT "[f, t, %s]}, s[[2]] === True && "
        "FreeQ[Together[f - Risch`Derivation[s[[1]], %s]], t]]]",
        g0, deriv, base, deriv, deriv);
    run_test(buf, "True");
}

static void test_iht_constructed(void) {
    /* Rational (no new logs): pure polynomial, reduced (special) poles, mixed. */
    run_iht_constructed("t", "0", HT);
    run_iht_constructed("t^2", "x", HT);
    run_iht_constructed("t^3", "0", HT);                 /* ∫ tan^3, elementary */
    run_iht_constructed("1/(t^2 + 1)", "1", HT);
    run_iht_constructed("t/(t^2 + 1)", "x^2", HT);
    run_iht_constructed("x t", "1/x", HT);
    run_iht_constructed("x^2/(t^2 + 1)^2", "x", HT);     /* repeated special pole */
    run_iht_constructed("t^2 + t/(t^2 + 1)", "x^3", HT);
    /* Example 5.10.3 (elementary reduced element). */
    run_iht_reduces("(t^5 + t^3 + x^2 t + 1)/(t^2 + 1)^3", HT);
    /* Genuine constant-residue logarithms (exercise ResidueReduce). */
    run_iht_reduces("(1 + t^2)/(t - 2)", HT);
    run_iht_reduces("(1 + t^2)/(t - 2) + (1 + t^2)/(t + 3)", HT);
}

static void test_iht_scaling(void) {
    /* Scaled monomial tan(2x). */
    run_iht_constructed("t", "x", HT2);
    run_iht_constructed("t^2/(t^2 + 1)", "1", HT2);
    run_iht_constructed("t^3 + t", "x^2", HT2);
    /* Half-angle monomial tan(x/2). */
    run_iht_constructed("t/(t^2 + 1)", "x", HALF);
    run_iht_constructed("t^2", "1", HALF);
    run_iht_reduces("(1 + t^2)/(2 (t - 3))", HALF);
}

/* ---- Robustness: non-hypertangent + malformed input ------------------ */
static void test_iht_robustness(void) {
    /* IntegrateHypertangent declines non-hypertangent monomials (exp, log). */
    run_test("Head[" IHT "[1/(t^2 + 1), t, {x -> 1, t -> t}]] === " IHT, "True");
    run_test("Head[" IHT "[t, t, {x -> 1, t -> 1/x}]] === " IHT, "True");
    /* Wrong arity / malformed derivation. */
    run_test("Head[" IHT "[t, t]] === " IHT, "True");
    run_test("Head[" IHT "[t, t, {x, t}]] === " IHT, "True");
    /* ResidueReduce is a general (§5.6) algorithm — it declines only on genuinely
     * malformed input (non-symbol monomial, wrong arity), not on non-hypertangent
     * monomials. */
    run_test("Head[" RR "[1/(t - 2), t]] === " RR, "True");
    run_test("Head[" RR "[1/(t - 2), 5, " HT "]] === " RR, "True");
}

/* ---- IntegrateHypertanh (hyperbolic tangent case, special t^2-1) ------ */
/* t = Tanh(x): Dt = 1 - t^2;  t = Tanh(2x): Dt = 2(1 - t^2). */
#define TD "{x -> 1, t -> 1 - t^2}"
#define TD2 "{x -> 1, t -> 2 (1 - t^2)}"
#define IHTH "Risch`IntegrateHypertanh"

static void test_hypertanh_values(void) {
    /* ∫Tanh = Log[Cosh] : g = -1/2 Log[1-t^2], beta True. */
    run_zero(IHTH "[t, t, " TD "][[1]] - (-1/2 Log[1 - t^2])");
    run_test(IHTH "[t, t, " TD "][[2]]", "True");
    /* ∫Tanh^2 : g = -t (leftover base +1). */
    run_test(IHTH "[t^2, t, " TD "]", "List[Times[-1, t], True]");
    /* ∫Sech^2 = Tanh : Sech^2 = 1 - t^2 -> g = t. */
    run_test(IHTH "[1 - t^2, t, " TD "]", "List[t, True]");
    /* Rejects a circular (t^2+1) derivation — that is the hypertangent case. */
    run_test("Head[" IHTH "[t, t, {x -> 1, t -> 1 + t^2}]] === " IHTH, "True");
    /* Malformed. */
    run_test("Head[" IHTH "[t, t]] === " IHTH, "True");
}

/* beta True and f - D[g] free of t (reduced to k). */
static void run_ihth_reduces(const char* f, const char* deriv) {
    char buf[4096];
    snprintf(buf, sizeof buf,
        "With[{s = " IHTH "[%s, t, %s]}, s[[2]] === True && "
        "FreeQ[Together[(%s) - Risch`Derivation[s[[1]], %s]], t]]",
        f, deriv, f, deriv);
    run_test(buf, "True");
}
/* Constructed f = D[g0] + base (g0 in k(t), base in k) must reduce (beta True). */
static void run_ihth_constructed(const char* g0, const char* base, const char* deriv) {
    char buf[6000];
    snprintf(buf, sizeof buf,
        "With[{f = Risch`Derivation[%s, %s] + (%s)}, "
        "With[{s = " IHTH "[f, t, %s]}, s[[2]] === True && "
        "FreeQ[Together[f - Risch`Derivation[s[[1]], %s]], t]]]",
        g0, deriv, base, deriv, deriv);
    run_test(buf, "True");
}

static void test_hypertanh_reduced_and_constructed(void) {
    /* Reduced (t^2-1 pole) via the decoupled real coupled system. */
    run_ihth_reduces("1/(1 - t^2)", TD);            /* = Cosh^2 */
    run_ihth_reduces("1/(1 - t^2)^2", TD);
    run_ihth_reduces("t/(1 - t^2)^2", TD);
    run_ihth_reduces("(t + x)/(1 - t^2)^2", TD);
    /* Constructed: polynomials, reduced poles, mixed, scaled monomial. */
    run_ihth_constructed("t", "x", TD);
    run_ihth_constructed("t^2", "1", TD);
    run_ihth_constructed("t^3", "0", TD);
    run_ihth_constructed("1/(1 - t^2)", "x^2", TD);
    run_ihth_constructed("t/(1 - t^2)", "1", TD);
    run_ihth_constructed("x t^2 + t", "x^3", TD);
    run_ihth_constructed("t", "x", TD2);            /* Tanh(2x), eta = 2 */
    run_ihth_constructed("t^2/(1 - t^2)", "1", TD2);
    /* Non-elementary: ∫ x Tanh^3 has no elementary integral (the x tanh
     * obstruction lifted) -> beta False. */
    run_test(IHTH "[x t^3, t, " TD "][[2]]", "False");
}

int main(void) {
    core_init();

    TEST(test_polynomial_reduce);
    TEST(test_integrate_hypertangent_poly);
    TEST(test_hypertangent_higher_degree);
    TEST(test_hypertangent_scaling);
    TEST(test_hypertangent_robustness);
    TEST(test_residue_reduce);
    TEST(test_residue_reduce_diffback);
    TEST(test_iht_book);
    TEST(test_iht_constructed);
    TEST(test_iht_scaling);
    TEST(test_iht_robustness);
    TEST(test_hypertanh_values);
    TEST(test_hypertanh_reduced_and_constructed);

    printf("All risch_hypertangent tests passed.\n");
    return 0;
}
