/* test_integrate_roach.c
 *
 * Integrate[..., Method -> "Roach"] (src/internal/RoachIntegrate.m): definite
 * integrals by Meijer G functions -- every factor a G function (PBM vol. 3
 * sec. 8.4), the one- and two-G integration theorems under the PBM 2.24.1
 * conditions, G -> named functions by Roach's Formula algorithm.
 *
 * Cases are the integrals of Kelly Roach, "Solving Integrals with the Quantum
 * Computer Algebra System" (2008), and the ISSAC '97 gallery. Each answer is
 * checked numerically at a sample point against a split-interval NIntegrate,
 * and some are also pinned to Roach's printed form.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* split NIntegrate reference, defined once */
static void define_ref(void) {
    Expr* e = parse_expression(
        "RoachTest`splitN[f_, {x_, a_, b_}] := Module[{p = If[b === Infinity,"
        " Join[{a}, a + {1, 4, 16, 64, 256}, {Infinity}], Subdivide[a, b, 4]]},"
        " Total[Table[NIntegrate[f, {x, p[[i]], p[[i + 1]]}], {i, Length[p] - 1}]]]");
    Expr* r = evaluate(e);
    expr_free(e); expr_free(r);
}

/* Integrate[f, {x, a, b}, Assumptions -> as, Method -> "Roach"] closes (no
 * Integrate left) and agrees with the split NIntegrate at `sub` to 1e-7. */
static void check_roach(const char* f, const char* rng, const char* as, const char* sub) {
    char buf[4096];
    snprintf(buf, sizeof buf,
        "With[{r = Integrate[%s, %s, Assumptions -> %s, Method -> \"Roach\"]},"
        " FreeQ[r, Integrate] && With[{v = N[(r /. ConditionalExpression[w_, _] :> w) /. %s],"
        " ref = RoachTest`splitN[Evaluate[(%s) /. %s], Evaluate[%s /. %s]]},"
        " Abs[v - ref] < 10^-7 (1 + Abs[ref])]]",
        f, rng, as, sub, f, sub, rng, sub);
    assert_eval_eq(buf, "True", 0);
}

/* Roach 2008, p. 3-4 */
static void test_paper(void) {
    check_roach("Exp[-al x] BesselJ[0, be x]^2", "{x, 0, Infinity}", "al > 0 && be > 0", "{al -> 2, be -> 1/3}");
    check_roach("(a^2 - x^2)^(n - 1/2)", "{x, 0, a}", "a > 0 && n > 0", "{a -> 2, n -> 1/3}");
    check_roach("Exp[-b x - x^2/(4 a)]", "{x, 0, Infinity}", "a > 0 && b > 0", "{a -> 2, b -> 1/3}");
    check_roach("x Exp[-x n]/(Exp[x] + 1)", "{x, 0, Infinity}", "n > 0", "{n -> 2/3}");
    check_roach("ArcSin[x] x^(2 n - 1)", "{x, 0, 1}", "n > 0", "{n -> 2/3}");
    check_roach("x^(a - 1)/Sqrt[x^b + 1]", "{x, 0, Infinity}", "b > 0 && a > 0 && a < b/2", "{a -> 1/3, b -> 3}");
    check_roach("x Exp[-2 x n]/(Exp[x] + 1)", "{x, 0, Infinity}", "n > 0", "{n -> 2/3}");
    check_roach("x Exp[-3 x]/(Exp[-x] + 1)", "{x, 0, Infinity}", "True", "{}");
    check_roach("x^(b - 1) (a x + 1)^(-c)", "{x, 0, Infinity}", "a > 0 && b > 0 && c > b", "{a -> 2, b -> 1/3, c -> 7/5}");
    check_roach("(SinIntegral[a x] - Pi/2) Exp[-b x]", "{x, 0, Infinity}", "a > 0 && b > 0", "{a -> 2, b -> 1/3}");
    check_roach("x^(a - 1) (b + x)^(-n)", "{x, 0, Infinity}", "b > 0 && a > 0 && n > a", "{a -> 1/3, b -> 2, n -> 7/5}");
    check_roach("BesselJ[1, a x] Exp[-x^2 b]", "{x, 0, Infinity}", "a > 0 && b > 0", "{a -> 2, b -> 1/3}");
    check_roach("Exp[-a x] x^n/(b + x)", "{x, 0, Infinity}", "a > 0 && b > 0 && n > -1", "{a -> 2, b -> 1/3, n -> 2/3}");
    check_roach("Log[x + 1] Cos[a x]/x", "{x, 0, Infinity}", "a > 0", "{a -> 2}");
    check_roach("x Erfc[a x] Exp[-x^2 b]", "{x, 0, Infinity}", "a > 0 && b > 0", "{a -> 2, b -> 1/3}");
    check_roach("x/(Exp[3 x] - 1)^(1/3)", "{x, 0, Infinity}", "True", "{}");
    check_roach("x Exp[-x]/Sqrt[Exp[x] - 1]", "{x, 0, Infinity}", "True", "{}");
    check_roach("x^2 Csch[x] Exp[-2 x n]", "{x, 0, Infinity}", "n > 0", "{n -> 2/3}");
}

/* t = E^-x with Log[t]^k by differentiation in the exponent, and the Lerch
 * term for 1/(1 -+ t^m): Roach's printed values */
static void test_exponential(void) {
    assert_eval_eq("Integrate[x/(Exp[3 x] - 1)^(1/3), {x, 0, Infinity}, Method -> \"Roach\"] ==="
                   " Pi (Pi + 3 Sqrt[3] Log[3])/27", "True", 0);
    assert_eval_eq("Integrate[x Exp[-x]/Sqrt[Exp[x] - 1], {x, 0, Infinity}, Method -> \"Roach\"] ==="
                   " Pi Log[2] - Pi/2", "True", 0);
    assert_eval_eq("Integrate[x^2 Csch[x] Exp[-2 x n], {x, 0, Infinity}, Assumptions -> n > 0,"
                   " Method -> \"Roach\"] === Zeta[3, n + 1/2]/2", "True", 0);
}

/* Roach's printed forms */
static void test_forms(void) {
    assert_eval_eq("Integrate[(SinIntegral[a x] - Pi/2) Exp[-b x], {x, 0, Infinity},"
                   " Assumptions -> a > 0 && b > 0, Method -> \"Roach\"] === (ArcTan[a/b] - Pi/2)/b", "True", 0);
    assert_eval_eq("Integrate[BesselJ[1, a x] Exp[-x^2 b], {x, 0, Infinity},"
                   " Assumptions -> a > 0 && b > 0, Method -> \"Roach\"] === (1 - E^(-a^2/(4 b)))/a", "True", 0);
    assert_eval_eq("Integrate[x^(b - 1) (a x + 1)^(-c), {x, 0, Infinity},"
                   " Assumptions -> a > 0 && b > 0 && c > b, Method -> \"Roach\"] ==="
                   " Gamma[b] Gamma[c - b]/(a^b Gamma[c])", "True", 0);
}

/* a symbolic finite upper limit z: the rescale x -> z y gives Log[z y], and a
 * G-form carrying theta(1 - z y) must not stand in for theta(1 - y) -- it
 * returned the value at z = 1 (Mathar M6460), which the numeric check passed
 * because FindInstance chose z -> 1 */
static void test_finite_symbolic_limit(void) {
    check_roach("Log[t]/(1 + a^2 t^2)", "{t, 0, z}", "a > 0 && z > 0", "{a -> 7/3, z -> 7/4}");
    check_roach("Log[t]/(1 + a^2 t^2)", "{t, 0, z}", "a > 0 && z > 0", "{a -> 7/4, z -> 5/4}");
    assert_eval_eq("FreeQ[Integrate[Log[t]/(1 + a^2 t^2), {t, 0, z},"
                   " Assumptions -> a > 0 && z > 0, Method -> \"Roach\"], z]", "False", 0);
}

/* Roach 2008 #18: the convolution is G^{3,4}_{4,4} at 1 in the log case,
 * reduced by the residue sum of its trigonometric Mellin-Barnes integrand
 * (Hrep`unitG). #24: a growing exponential against Ei(-x^2), by analytic
 * continuation from the decaying case (Roach`int2grow). */
static void test_unit_and_growing(void) {
    check_roach("ArcCot[x]/(x + 1)", "{x, 0, Infinity}", "True", "{}");
    check_roach("ArcTan[x]/(x (x + 1))", "{x, 0, Infinity}", "True", "{}");
    check_roach("x ArcCot[x]/((1 + x) (1 + x^2))", "{x, 0, Infinity}", "True", "{}");
    assert_eval_eq("Integrate[ArcCot[x]/(x + 1), {x, 0, Infinity}, Method -> \"Roach\"] ==="
                   " Catalan + Pi Log[2]/4", "True", 0);
    assert_eval_eq("Integrate[ArcCot[x]/(x + 1), {x, 0, Infinity}] === Catalan + Pi Log[2]/4", "True", 0);
    check_roach("ExpIntegralEi[-x^2] Exp[a^2 x^2]", "{x, 0, Infinity}", "0 < a < 1", "{a -> 1/3}");
    check_roach("ExpIntegralEi[-x] Exp[b x] x", "{x, 0, Infinity}", "0 < b < 1", "{b -> 1/2}");
    check_roach("Erfc[x] Exp[b x^2]", "{x, 0, Infinity}", "0 < b < 1", "{b -> 1/2}");
    assert_eval_eq("Integrate[ExpIntegralEi[-x^2] Exp[a^2 x^2], {x, 0, Infinity},"
                   " Assumptions -> 0 < a < 1, Method -> \"Roach\"] === -Sqrt[Pi] ArcSin[a]/a", "True", 0);
    assert_eval_eq("Integrate[ExpIntegralEi[-x^2] Exp[a^2 x^2], {x, 0, Infinity},"
                   " Assumptions -> 0 < a < 1] === -Sqrt[Pi] ArcSin[a]/a", "True", 0);
    /* divergent: the exponential grows faster than Ei(-x^2) decays */
    assert_eval_eq("Head[Integrate[ExpIntegralEi[-x^2] Exp[a^2 x^2], {x, 0, Infinity},"
                   " Assumptions -> a > 1, Method -> \"Roach\"]]", "Integrate", 0);
}

/* the stage declines (stays unevaluated) rather than guessing */
static void test_declines(void) {
    assert_eval_eq("Head[Integrate[Sin[x^x], {x, 0, 1}, Method -> \"Roach\"]]", "Integrate", 0);
    /* divergent: x^(b-1) (1 + x)^-c with c < b */
    assert_eval_eq("Head[Integrate[x^(b - 1) (x + 1)^(-c), {x, 0, Infinity},"
                   " Assumptions -> b > 0 && c > 0 && c < b, Method -> \"Roach\"]]", "Integrate", 0);
}

/* the stage leaves no global state behind (pattern capture, bindings) */
static void test_hygiene(void) {
    assert_eval_eq("Integrate[Exp[-al x] BesselJ[0, be x]^2, {x, 0, Infinity},"
                   " Assumptions -> al > 0 && be > 0, Method -> \"Roach\"]; {ValueQ[x], ValueQ[al], ValueQ[be]}",
                   "{False, False, False}", 0);
    assert_eval_eq("Check[Integrate[x^(b - 1) (a x + 1)^(-c), {x, 0, Infinity},"
                   " Assumptions -> a > 0 && b > 0 && c > b, Method -> \"Roach\"]; ok, failed]", "ok", 0);
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();
    define_ref();
    TEST(test_paper);
    TEST(test_forms);
    TEST(test_exponential);
    TEST(test_finite_symbolic_limit);
    TEST(test_unit_and_growing);
    TEST(test_declines);
    TEST(test_hygiene);
    printf("All Integrate Roach tests passed!\n");
    return 0;
}
