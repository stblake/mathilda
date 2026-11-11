/* test_integrate_meijerconv.c
 *
 * The MeijerConvolution definite-integration stage
 * (src/internal/MeijerConvolution.m, wired in src/calculus/integrate.c):
 * C x^alpha K1(x) K2(x) on (0, Infinity) by the Meijer G convolution theorem,
 * kernels BesselJ / BesselY / BesselK, E^(-p x), E^(-p x^2), Cos, Sin,
 * (a^2 + x^2)^-r at arbitrary scales.  Reference values are Mathematica
 * NIntegrate at 30 digits (closed forms where Mathematica gives them).
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

static void check_eq(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n",
                input, expected, s);
    }
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(p);
    expr_free(r);
}

static void check_value(const char* integral, const char* ref) {
    char buf[1024];
    snprintf(buf, sizeof buf, "Chop[N[(%s) - (%s)], 10^-9]", integral, ref);
    check_eq(buf, "0");
}

static void test_two_kernels(void) {
    check_value("Integrate[x^(3/5) E^(-13/10 x) BesselJ[7/10, 9/10 x], {x, 0, Infinity}]",
                "0.2592479223723994");
    check_value("Integrate[x^2 BesselK[0, x] BesselK[1, 2 x], {x, 0, Infinity}]",
                "0.1413987469155312");
    /* different frequencies, amplitude ~1/x: conditionally convergent; = 0 */
    check_value("Integrate[BesselJ[1, x] BesselY[1, 3 x], {x, 0, Infinity}]", "0");
    /* x J1(x) Y1(3x) has a non-decaying amplitude: not convergent, must decline */
    check_eq("Head[Integrate[x BesselJ[1, x] BesselY[1, 3 x], {x, 0, Infinity}, "
             "Method -> \"MeijerConvolution\"]]", "Integrate");
    check_value("Integrate[x BesselJ[0, x] E^(-x^2/2), {x, 0, Infinity}]", "E^(-1/2)");
    check_value("Integrate[BesselJ[0, 2 x] Cos[x], {x, 0, Infinity}]", "1/Sqrt[3]");
    check_value("Integrate[x BesselJ[0, x]/(1 + x^2), {x, 0, Infinity}]", "BesselK[0, 1]");
}

static void test_symbolic(void) {
    /* Laplace transform of J_nu: a^nu/(Sqrt[a^2+p^2] (p + Sqrt[a^2+p^2])^nu) */
    check_eq("Chop[N[(Integrate[BesselJ[nu, a x] E^(-p x), {x, 0, Infinity}, "
             "Assumptions -> a > 0 && p > 0 && nu > 0] - a^nu/(Sqrt[a^2 + p^2] "
             "(p + Sqrt[a^2 + p^2])^nu)) /. {a -> 2/3, p -> 5/4, nu -> 3/7}], 10^-9]", "0");
}

static void test_declines(void) {
    /* divergent at 0: x^-2 J0 */
    check_eq("Head[Integrate[BesselJ[0, x] BesselJ[0, 2 x]/x^2, {x, 0, Infinity}, "
             "Method -> \"MeijerConvolution\"]]", "Integrate");
    /* equal-frequency product without enough decay */
    check_eq("Head[Integrate[BesselJ[0, x]^2, {x, 0, Infinity}, "
             "Method -> \"MeijerConvolution\"]]", "Integrate");
    /* (x^2 - b^2)^-1 has real poles: the (a + x^2)^-r kernel needs Re[a] > 0,
       and the refuted side condition must decline (it returned a wrong
       finite value for these principal-value integrals before) */
    check_eq("Head[Integrate[Sin[a x]/(x^2 - b^2), {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"MeijerConvolution\"]]", "Integrate");
    check_eq("Head[Integrate[x Cos[a x]/(x^2 - b^2), {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"MeijerConvolution\"]]", "Integrate");
    /* ... while Re[a] > 0 still closes */
    check_eq("Chop[N[(Integrate[Cos[a x]/(x^2 + b^2), {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"MeijerConvolution\"] - "
             "Pi E^(-a b)/(2 b)) /. {a -> 2, b -> 3}, 20], 10^-12]", "0");
}

/* The stage answers in MeijerG; Integrate's final pass
 * (fe_reduce_definite_result) rewrites it to named functions. The targets are
 * the reference system's values; the 30-digit checks confirm the MeijerG the
 * stage produced really is that function (arb numerics on the left via an
 * MPFR argument, so no exact reduction is involved). */
static void test_named_results(void) {
    check_eq("Integrate[E^-x Sin[x]/x, {x, 0, Infinity}]", "1/4 Pi");
    check_eq("Integrate[x BesselJ[0, a x]/(x^2 + c^2), {x, 0, Infinity}, "
             "Assumptions -> a > 0 && c > 0]", "BesselK[0, a c]");
    check_eq("Integrate[BesselJ[0, a x] Sin[b x]/x, {x, 0, Infinity}, "
             "Assumptions -> b > a > 0]", "1/2 Pi");
    check_eq("Abs[N[1/2 MeijerG[{{1/2}, {1}}, {{0, 1/2}, {}}, N[1, 50]] - Pi/4, 40]] < 10^-30", "True");
    check_eq("Abs[N[(2 MeijerG[{{0, 0}, {}}, {{}, {}}, N[4/(a^2 c^2), 50]])/(a^2 c^2) - "
             "BesselK[0, a c] /. {a -> 3/2, c -> 5/7}, 40]] < 10^-30", "True");
    check_eq("Abs[N[1/2 Sqrt[Pi] MeijerG[{{1}, {1}}, {{1/2}, {0}}, N[b^2/a^2, 50]] - Pi/2 /. "
             "{a -> 2/3, b -> 5/4}, 40]] < 10^-30", "True");
    /* b < a is the other branch of the same Piecewise: ArcSin[b/a] */
    check_eq("Integrate[BesselJ[0, a x] Sin[b x]/x, {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0 && b < a]", "ArcSin[b/a]");
}

/* Two oscillating kernels whose frequencies cross (Weber-Schafheitlin):
 * the q < 1 and q > 1 closed forms do not cover q == 1 itself. There the
 * integral converges iff the non-oscillating tail x^(alpha+e1+e2) decays
 * faster than 1/x, and is then continuous with the theorem's value at
 * argument 1; otherwise q == 1 must be excluded by the condition. It used
 * to fall through to Piecewise's implicit default 0 (G&R 3.741.3, 3.828.2,
 * 6.693.1, 6.693.2; issue #132). References: mpmath quadosc, 20 digits. */
static void test_frequency_crossing(void) {
    /* G&R 3.741.3: Pi min(a, b)/2, the boundary a == b included */
    check_value("Integrate[Sin[a x] Sin[b x]/x^2, {x, 0, Infinity}, "
                "Assumptions -> 0 < a <= b] /. {a -> 7/4, b -> 7/4}", "7 Pi/8");
    check_value("Integrate[Sin[a x] Sin[b x]/x^2, {x, 0, Infinity}, "
                "Assumptions -> 0 < a <= b] /. {a -> 5/8, b -> 7/4}", "5 Pi/16");
    /* G&R 3.828.2: side undecided -- both branches plus the boundary */
    check_value("Integrate[Sin[a x] Sin[b x]/x^2, {x, 0, Infinity}, "
                "Assumptions -> a > 0 && b > 0] /. {a -> 5/8, b -> 5/8}", "5 Pi/16");
    check_value("Integrate[Sin[a x] Sin[b x]/x^2, {x, 0, Infinity}, "
                "Assumptions -> a > 0 && b > 0] /. {a -> 7/3, b -> 1/5}", "Pi/10");
    /* G&R 6.693.1: Sin[nu Pi/2]/nu at beta == alpha */
    check_value("Integrate[BesselJ[nu, alpha x] Sin[beta x]/x, {x, 0, Infinity}, "
                "Assumptions -> beta >= alpha && Re[nu] > -1 && alpha > 0] /. "
                "{alpha -> 5/6, beta -> 5/6, nu -> 5/4}", "Sin[5 Pi/8]/(5/4)");
    check_value("Integrate[BesselJ[nu, alpha x] Sin[beta x]/x, {x, 0, Infinity}, "
                "Assumptions -> beta <= alpha && Re[nu] > -1 && alpha > 0 && beta > 0] /. "
                "{alpha -> 5/4, beta -> 5/4, nu -> 3/8}", "Sin[3 Pi/16]/(3/8)");
    /* G&R 6.693.2: Cos[nu Pi/2]/nu at beta == alpha */
    check_value("Integrate[BesselJ[nu, alpha x] Cos[beta x]/x, {x, 0, Infinity}, "
                "Assumptions -> beta <= alpha && Re[nu] > 0 && alpha > 0 && beta > 0] /. "
                "{alpha -> 5/6, beta -> 5/6, nu -> 5/4}", "Cos[5 Pi/8]/(5/4)");
    check_value("Integrate[BesselJ[1, a x] BesselJ[1, b x]/x, {x, 0, Infinity}, "
                "Assumptions -> a > 0 && b > 0] /. {a -> 7/4, b -> 7/4}", "1/2");
    /* divergent at the crossing (Sin^2/x, J0 Sin ~ 1/Sqrt[b^2 - a^2]): the
       boundary is excluded, never a finite number */
    check_eq("Integrate[Sin[a x] Sin[b x]/x, {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0] /. {a -> 1, b -> 1}", "Undefined");
    check_eq("Integrate[BesselJ[0, a x] Sin[b x], {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0] /. {a -> 1, b -> 1}", "Undefined");
    check_value("Integrate[BesselJ[0, a x] Sin[b x], {x, 0, Infinity}, "
                "Assumptions -> a > 0 && b > 0] /. {a -> 5/8, b -> 7/4}",
                "1/Sqrt[(7/4)^2 - (5/8)^2]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_two_kernels);
    TEST(test_symbolic);
    TEST(test_declines);
    TEST(test_named_results);
    TEST(test_frequency_crossing);

    printf("All Integrate MeijerConvolution tests passed!\n");
    return 0;
}
