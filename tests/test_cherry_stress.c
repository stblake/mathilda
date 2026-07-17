/* test_cherry_stress.c — broad soundness / robustness stress battery for the
 * Cherry special-function integrators (ei, erf, li, dilog) and the B3/B2/A1
 * additions (the C0 seam, multi-term ei, complex constants).
 *
 * The Cherry engines are correct-by-construction (every emission passes the exact
 * rt_verify_antideriv gate before it is returned), so this suite's job is:
 *   (1) INDEPENDENT soundness re-verification — a large corpus, each closed answer
 *       numerically diff-backed to 0 at several safe points (a hole in the internal
 *       gate would surface here as a nonzero diff-back);
 *   (2) ROBUSTNESS — no crash / hang / wrong form across diverse and adversarial
 *       inputs;
 *   (3) DECLINE-SAFETY — genuinely non-elementary or currently-deferred inputs
 *       return unevaluated, never a wrong closed form;
 *   (4) the DECISION predicates (LiElementaryQ / SigmaDecomposition).
 *
 * Diff-back is checked NUMERICALLY at points all > sqrt(2) (so every Log/li
 * argument in the corpus is real and positive — avoids the branch cut of the
 * PowerExpand li convention the engine emits under), which the exact tower gate is
 * insensitive to but a naive Simplify diff-back is not.
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
#include <stdbool.h>

static Expr* eval_str(const char* s) {
    Expr* e = parse_expression(s);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static bool eval_is(const char* s, const char* expected) {
    Expr* r = eval_str(s);
    char* got = expr_to_string_fullform(r);
    bool ok = strcmp(got, expected) == 0;
    if (!ok) printf("  [%s] -> %s (expected %s)\n", s, got, expected);
    free(got);
    expr_free(r);
    return ok;
}

/* Safe evaluation grid: all > sqrt(2), avoiding the small integers / half-integers
 * that are poles of the corpus denominators. */
#define PTS "{3/2, 7/4, 11/4, 13/6, 5/2}"

/* f closes under the recursive Risch engine AND the closed form diff-backs
 * numerically to 0 (soundness) — the core stress assertion. */
static void closes(const char* f) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "Head[Integrate[%s, x, Method -> \"RischTranscendental\"]] =!= Integrate", f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: expected to close", f);
    snprintf(buf, sizeof(buf),
        "With[{r = Integrate[%s, x, Method -> \"RischTranscendental\"]},"
        " Max[Abs[Table[N[(D[r, x] - (%s)) /. x -> p, 25], {p, %s}]]] < 10^-10]",
        f, f, PTS);
    ASSERT_MSG(eval_is(buf, "True"), "%s: numeric diff-back not ~0 (SOUNDNESS)", f);
}

/* f is left unevaluated by the recursive Risch engine (clean decline). */
static void declines(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "Head[Integrate`RischTranscendental[%s, x]] === Integrate`RischTranscendental", f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: expected a clean decline", f);
}

/* Universal soundness: whatever the engine returns, if it closed it must diff-back
 * to 0.  Applied to adversarial inputs where we do not assert closes-vs-declines. */
static void sound(const char* f) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "With[{r = Integrate[%s, x, Method -> \"RischTranscendental\"]},"
        " Head[r] === Integrate ||"
        " Max[Abs[Table[N[(D[r, x] - (%s)) /. x -> p, 25], {p, %s}]]] < 10^-10]",
        f, f, PTS);
    ASSERT_MSG(eval_is(buf, "True"), "%s: closed to a WRONG form (SOUNDNESS)", f);
}

/* ---------------- ei: linear exponent / linear denominator ---------------- */
static void test_ei_linear(void) {
    closes("E^x/(x + 1)");
    closes("E^(2 x + 1)/(x - 1)");
    closes("E^(-x)/(2 x + 3)");
    closes("E^(3 x - 2)/(x + 5)");
    closes("E^x/(2 x - 1)");
    closes("5 E^(x + 2)/(3 x + 1)");
}

/* ---------------- ei: rational exponents ---------------- */
static void test_ei_ratexp(void) {
    closes("E^(1/x)");                 /* Ex 5.1 */
    closes("E^(1/x)/x");
    closes("E^(1/x)/x^2");
}

/* ---------------- ei: higher-multiplicity poles ---------------- */
static void test_ei_polemult(void) {
    closes("E^x/(x + 1)^2");           /* d8 */
    closes("E^x/(x + 1)^3");
    closes("E^x/(x + 1)^4");
    closes("E^x/(x - 1)^2");
    closes("E^x/(x - 2)^3");
    closes("E^x/x^2");
    closes("E^x/x^3");
}

/* ---------------- ei: multiple distinct rational poles ---------------- */
static void test_ei_multipole(void) {
    closes("E^x/((x - 1)(x - 2))");
    closes("E^x/((x - 1)(x - 2)(x - 3))");
    closes("(x^2 + 3) E^x/(x^2 + 3 x + 2)");   /* d11 */
    closes("E^x/((2 x - 1)(3 x + 2))");
    closes("(x + 1) E^x/((x - 2)(x - 3))");
}

/* ---------------- ei: real algebraic constants (C = C-bar, real) ---------------- */
static void test_ei_realalg(void) {
    closes("E^x/(x^2 - 2)");           /* Cherry p.894, sqrt2 */
    closes("E^x/(x^2 - 5)");           /* sqrt5 */
    closes("E^x/(x^2 - x - 1)");       /* golden ratio */
    closes("E^x/(x^2 - 3 x + 1)");
}

/* ---------------- ei: complex constants (A1) — fields Solve cracks ---------------- */
static void test_ei_complex(void) {
    closes("E^x/(x^2 + 1)");           /* Q(i)     */
    closes("E^x/(x^2 + 9)");           /* Q(i)     */
    closes("E^x/(x^2 + 2)");           /* Q(i sqrt2), no shift */
    closes("E^x/(x^2 + 3)");           /* Q(i sqrt3), no shift */
    closes("E^x/(x^2 + 6)");           /* Q(i sqrt6), no shift */
    closes("E^x/(x^2 + 2 x + 2)");     /* Q(i),   shifted */
    closes("E^x/(x^2 + 2 x + 5)");     /* Q(i),   shifted */
    closes("E^x/(x^2 + x + 1)");       /* Q(i sqrt3), shifted */
    closes("(2 x + 1) E^x/(x^2 + x + 3)");  /* Q(i sqrt11), numerator */
    closes("x E^x/(x^2 + x + 1)");     /* numerator over Q(i sqrt3) */
    closes("(x^2 + 1) E^x/(x^2 + x + 1)");  /* d12 */
    /* Shifted complex centers over Q(i sqrt d): the direct Solve over the number
     * field fails (Together will not cancel the conjugate-linear factors back to the
     * real quadratic), so these go through the {1,chs}-basis NF fallback
     * (rt_cherry_ei_conjpair_nf). */
    closes("E^x/(x^2 + 2 x + 3)");     /* Q(i sqrt2), shifted */
    closes("E^x/(x^2 + 2 x + 7)");     /* Q(i sqrt6), shifted */
    closes("E^x/(x^2 + 4 x + 7)");     /* Q(i sqrt3), shifted center 2 */
    closes("E^x/(x^2 - 2 x + 3)");     /* Q(i sqrt2), positive shift */
    closes("(3 x + 1) E^x/(x^2 + 2 x + 3)");  /* NF fallback with numerator */
}

/* Constant exponent offset: E^(c + h(x)) = E^c E^(h(x)); the constant folds into
 * the cofactor (else it inflates deg(p) and defeats P2 recognition). */
static void test_ei_const_offset(void) {
    closes("E^(1/x + 2)");
    closes("E^((x - 1)/x)");           /* = E^(1 - 1/x) */
    closes("E^(1/x^2 + 3)");           /* offset + erf */
    closes("E^(1/x - 1)");
}

/* ---------------- ei: nonlinear exponents ---------------- */
static void test_ei_nonlinear(void) {
    closes("E^(x^2)/x");
    closes("E^(x^2 + x)");
    closes("E^(2 x^2)/x");
}

/* ---------------- ei: multi-term (B2, Thm 5.4 case b) ---------------- */
static void test_ei_multiterm(void) {
    closes("(E^x + E^(2 x))/(x - 1)");
    closes("(E^x + E^(2 x) + E^(3 x))/(x - 1)");
    closes("E^x (E^x + 1)/((x - 1)(x - 2))");
    closes("(E^x + E^(2 x))/((x - 1)(x - 2))");
    closes("(E^x + E^(2 x))/(x^2 - 2)");
}

/* ---------------- erf (q a perfect square) ---------------- */
static void test_erf(void) {
    closes("E^(1/x^2)");
    closes("(1/x + 1/x^2) E^(1/x^2)");   /* Ex 5.2 */
    closes("E^(1/x^2)/x");
}

/* ---------------- li ---------------- */
static void test_li(void) {
    closes("1/Log[x]");
    closes("x/Log[x]^2");                 /* d1 */
    closes("x/Log[x]");
    closes("x^2/Log[x]");
    closes("x^2/Log[x + 1]");             /* d3, multi-li */
    closes("x^3/Log[x^2 - 1]");           /* Ex 5.1, two-log tower */
    closes("1/(Log[x] + 3)");             /* d2, rescale */
    closes("(Log[x]^2 + 3)/(Log[x]^2 + 3 Log[x] + 2)");  /* d4 */
}

/* ---------------- dilog ---------------- */
static void test_dilog(void) {
    closes("Log[x]/(1 + x)");
    closes("Log[x]/(1 - x)");
    closes("Log[x]/(x^2 - 1)");
    closes("Log[2 + x]/x");               /* transcendental spacing (monic kernel) */
}

/* ---------------- B3: the ExpIntegralEiResultant primitive ---------------- */
static void test_b3_resultant(void) {
    /* p.894 g1=x^2-2, p=x, q=1  ->  a^2-2 (roots +-sqrt2). */
    ASSERT_MSG(eval_is("Integrate`ExpIntegralEiResultant[x^2 - 2, x, 1, a, x]",
        "Plus[-2, Power[a, 2]]"), "ExpIntegralEiResultant p.894");
    /* d8 g1=(x+1)^2 -> (a-1)^2 (double root a=1 -> ei(x+1)). */
    ASSERT_MSG(eval_is("Factor[Integrate`ExpIntegralEiResultant[(x + 1)^2, x, 1, a, x]]",
        "Power[Plus[-1, a], 2]"), "ExpIntegralEiResultant d8");
}

/* ---------------- decision predicates (A2 / B1) ---------------- */
static void test_decisions(void) {
    /* li-elementarity: Ex 5.1 is li-elementary (True), Ex 5.2 is not (False). */
    ASSERT_MSG(eval_is("Integrate`LiElementaryQ[x^3/Log[x^2 - 1], x]", "True"),
        "LiElementaryQ Ex 5.1 -> True");
    ASSERT_MSG(eval_is("Integrate`LiElementaryQ[x^2/Log[x^2 - 1], x]", "False"),
        "LiElementaryQ Ex 5.2 -> False");
    /* Sigma-decomposition: Ex 5.1 x^2/2 decomposes; Ex 5.2 x/2 does not. */
    ASSERT_MSG(eval_is("Head[Integrate`SigmaDecomposition[x^2/2, {x - 1, x + 1}, x]]"
        " === List", "True"), "SigmaDecomposition Ex 5.1 -> List");
    ASSERT_MSG(eval_is("Integrate`SigmaDecomposition[x/2, {x - 1, x + 1}, x]", "$Failed"),
        "SigmaDecomposition Ex 5.2 -> $Failed");
}

/* ---------------- decline-safety: must return unevaluated, never wrong ---------------- */
static void test_decline_safety(void) {
    /* complex pair mixed with a P2/reciprocal term (deferred over Q(i)) */
    declines("E^(1/x)/(x^2 + 1)");
    /* degree-3 constant tower: only the real root admitted -> clean decline */
    declines("E^x/(x^3 - 2)");
    /* mixed polynomial+rational exponent (genuinely hard / non-elementary shape) */
    declines("E^((x^2 + 1)/x)");
    declines("E^(2/x^2 + 1/x)");
    /* C-i: non-monic dilog kernel — cherry_dilog (Risch path) declines a non-monic
     * Log[a x + b] kernel (and a rational-root kernel generally); the FULL Integrate
     * cascade DOES close Log[2x+1]/(x+1), so users are unaffected.  This asserts the
     * cherry_dilog-internal decline only; update to closes() if that engine grows
     * rational-root interpolants. */
    declines("Log[2 x + 1]/(x + 1)");
}

/* ---------------- universal soundness on adversarial inputs ---------------- */
static void test_soundness_adversarial(void) {
    /* Whatever these do, a closed answer MUST diff-back to 0. */
    sound("E^x/(x^2 + 2 x + 3)");
    sound("E^(1/x + 2)");
    sound("E^((x^2 + 1)/x)");
    sound("(E^x + E^(3 x))/(x^2 - 3)");
    sound("(x^3 + 1) E^x/(x^2 + x + 1)");
    sound("E^x/(x^2 - 7)");
    sound("E^x/((x - 1)^2 (x - 2))");
    sound("E^(1/x)/(x + 1)");
    sound("x^4/Log[x + 1]");
    sound("Log[3 x + 2]/(x - 1)");
}

int main(void) {
    core_init();
    TEST(test_ei_linear);
    TEST(test_ei_ratexp);
    TEST(test_ei_polemult);
    TEST(test_ei_multipole);
    TEST(test_ei_realalg);
    TEST(test_ei_complex);
    TEST(test_ei_const_offset);
    TEST(test_ei_nonlinear);
    TEST(test_ei_multiterm);
    TEST(test_erf);
    TEST(test_li);
    TEST(test_dilog);
    TEST(test_b3_resultant);
    TEST(test_decisions);
    TEST(test_decline_safety);
    TEST(test_soundness_adversarial);
    printf("All Cherry stress tests passed.\n");
    return 0;
}
