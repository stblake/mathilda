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
    TEST(test_declines);
    TEST(test_hygiene);
    printf("All Integrate Roach tests passed!\n");
    return 0;
}
