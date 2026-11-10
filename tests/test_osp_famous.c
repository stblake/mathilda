/* Regression tests for the fixes driven by the math.osp.no "famous" problems
 * (benchmarks/osp_famous).  Each fix is a general algorithm, so each is tested
 * on the famous instance AND on siblings the famous problem does not name, and
 * every identity is checked NUMERICALLY against an independent value (a direct
 * series sum, NIntegrate, or a hand-derived closed form), never against
 * Mathilda's own rendering of the answer.
 *
 *   Solve      a {eqn} list in one variable no longer mutates its input
 *              (Solve[{x E^x == 1}, {x}] printed as Solve[{Null}, {x}]);
 *              Lambert-W inversion of k (a x + b) B^(p x + q) + c == 0 and of
 *              the sum form L(x) + k B^M(x) + c == 0, with the W_n branch family
 *              over Complexes and the real-branch count over Reals.
 *   DSolve     a condition at a regular singular point of the general solution
 *              (y(0) of Cos[x]/x, u'(0) of Log[r]) is fitted from the Laurent
 *              expansion: Lane-Emden n=0 and n=1 (n=1 SEGFAULTED: an infinite
 *              Equal was handed to Solve), Poiseuille flow (returned {}).
 *   Integrate  the residue Gaussian/affine families and the Ramanujan Gaussian
 *              family no longer read x^2 + x^-2 as the polynomial x^2
 *              (Integrate[E^(-x^2 - 1/x^2), {x, 0, Infinity}] was Sqrt[Pi]/2);
 *              a factored exponent E^(-(x - 1/x)^2) is expanded and its
 *              constant split off before the definite cascade gives up.
 *   Sum        a hypergeometric series whose leading terms vanish starts at its
 *              first nonzero term (Sum[k^10/k!, {k, 0, Infinity}] was 0).
 *   HypergeometricPFQ  an upper parameter exceeding a lower one by a positive
 *              integer n reduces through P(theta), theta = z d/dz.
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* `cond` must evaluate to True. */
static void assert_true(const char* cond) { assert_eval_eq(cond, "True", 0); }

/* N[a - b, 40] vanishes to 1e-25 relative: an exact identity checked numerically. */
static void assert_num_eq(const char* a, const char* b) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "With[{mrD = Quiet[N[(%s) - (%s), 40]], mrS = Quiet[N[%s, 40]]}, "
             "NumericQ[mrD] && Abs[mrD] <= 10^-25 Max[1, Abs[mrS]]]", a, b, b);
    assert_true(buf);
}

/* The expression does not evaluate to anything containing `head`. */
static void assert_free_of(const char* expr, const char* head) {
    char buf[4096];
    snprintf(buf, sizeof buf, "FreeQ[%s, %s]", expr, head);
    assert_true(buf);
}

/* ---- Solve: list-form input is not mutated ------------------------------ */

void test_solve_list_form_input_intact() {
    /* An equation nothing can solve, posed as {eqn}: the unevaluated result
     * must still carry the equation, not Null. */
    assert_true("FreeQ[Solve[{Sin[x] x^x == Cos[x]^Sqrt[x] + 1}, {x}], Null]");
    assert_true("FreeQ[Solve[{Sin[x] x^x == Cos[x]^Sqrt[x] + 1}, x], Null]");
    /* the bare and list forms agree */
    assert_true("Solve[{Log[y] == 3}, y] === Solve[Log[y] == 3, y]");
}

/* ---- Solve: Lambert W ---------------------------------------------------- */

void test_solve_lambert_product_form() {
    /* the omega constant, every branch over Complexes */
    assert_true("(" "Solve[{x Exp[x] == 1}, {x}, Complexes]" ") === (" "{{x -> ConditionalExpression[ProductLog[C[1], 1], Element[C[1], Integers]]}}" ")");
    /* principal branch by default, and Omega = 0.5671432904097838... */
    assert_true("(" "Quiet[Solve[x Exp[x] == 1, x]]" ") === (" "{{x -> ProductLog[1]}}" ")");
    assert_num_eq("x /. Quiet[Solve[x Exp[x] == 1, x]][[1]]", "0.56714329040978387299996866221035554975381578718651250813513");
    /* a base other than E */
    assert_true("(" "Solve[x 2^x == 3, x, Reals]" ") === (" "{{x -> ProductLog[3 Log[2]]/Log[2]}}" ")");
    /* Reals: two branches on [-1/e, 0), none below */
    assert_true("(" "Solve[x Exp[x] == -1/4, x, Reals]" ") === (" "{{x -> ProductLog[-1/4]}, {x -> ProductLog[-1, -1/4]}}" ")");
    assert_true("(" "Solve[x Exp[x] == -1, x, Reals]" ") === (" "{}" ")");
    /* a general linear factor and exponent: every returned root satisfies it */
    assert_true("With[{s = Quiet[Solve[(2 x + 1) Exp[3 x - 1] == 5, x]]}, "
                "Length[s] == 1 && Abs[N[(2 x + 1) Exp[3 x - 1] - 5 /. N[s[[1]], 40]]] < 10^-20]");
    assert_true("With[{s = Solve[(2 x + 1) Exp[3 x - 1] == 5, x, Complexes] /. C[1] -> -1}, "
                "Abs[N[(2 x + 1) Exp[3 x - 1] - 5 /. N[s[[1]] /. ConditionalExpression[v_, _] :> v, 40]]] < 10^-20]");
}

void test_solve_lambert_sum_form() {
    /* x == e^-x is x e^x == 1 */
    assert_true("(" "Quiet[Solve[x == Exp[-x], x]]" ") === (" "{{x -> ProductLog[1]}}" ")");
    assert_true("With[{s = Solve[x + 2^x == 0, x, Reals]}, "
                "Length[s] == 1 && Abs[N[x + 2^x /. N[s[[1]], 40]]] < 10^-20]");
    assert_true("With[{s = Solve[3 x - 1 == 2 Exp[x/2], x, Complexes]}, "
                "And @@ Table[Abs[N[3 x - 1 - 2 Exp[x/2] /. N[(s[[1]] /. ConditionalExpression[v_, _] :> v) "
                "/. C[1] -> n, 40]]] < 10^-20, {n, -2, 2}]]");
    /* not the Lambert shape: left alone */
    assert_true("Head[Quiet[Solve[x^2 Exp[x] + Sin[x] == 1, x]]] === Solve");
}

/* ---- DSolve: conditions at a regular singular point ----------------------- */

void test_dsolve_singular_point_conditions() {
    /* Lane-Emden n = 1: y = Sin[x]/x (used to segfault) */
    assert_num_eq("y[7/10] /. DSolve[{y''[x] + y[x] + 2 y'[x]/x == 0, y[0] == 1, y'[0] == 0}, y, x][[1]]",
                  "Sin[7/10]/(7/10)");
    /* Lane-Emden n = 0: y = 1 - x^2/6 */
    assert_true("DSolve[{y''[x] + 2 y'[x]/x + 1 == 0, y[0] == 1, y'[0] == 0}, y[x], x] === "
                "{{y[x] -> 1 - x^2/6}}");
    /* Poiseuille: u = G (R^2 - r^2)/(4 mu) (returned {}) */
    assert_true("Simplify[(u[r] /. DSolve[{u''[r] + u'[r]/r == -G/mu, u[R] == 0, u'[0] == 0}, u[r], r][[1]]) "
                "- G (R^2 - r^2)/(4 mu)] === 0");
    /* Euler equation: the x^-2 member has a pole at 0, so y(0) == 0 forces its
     * constant to 0 and y(1) == 5 fixes the other */
    assert_true("DSolve[{x^2 y''[x] + x y'[x] - 4 y[x] == 0, y[0] == 0, y[1] == 5}, y[x], x] === "
                "{{y[x] -> 5 x^2}}");
    /* a regular point is unaffected */
    assert_true("(" "DSolve[{y''[x] + y[x] == 0, y[0] == 1, y'[0] == 0}, y[x], x]" ") === (" "{{y[x] -> Cos[x]}}" ")");
}

/* ---- Integrate: essential singularity in the exponent --------------------- */

void test_integrate_cauchy_schlomilch() {
    /* the famous problem: Sqrt[Pi]/2 (was Sqrt[Pi] E^2/2) */
    assert_num_eq("Integrate[Exp[-(x - 1/x)^2], {x, 0, Infinity}]", "Sqrt[Pi]/2");
    assert_num_eq("Integrate[Exp[-(x - 2/x)^2], {x, 0, Infinity}]", "Sqrt[Pi]/2");
    /* the expanded form: Sqrt[Pi]/(2 E^2) (was Sqrt[Pi]/2) */
    assert_num_eq("Integrate[Exp[-x^2 - 1/x^2], {x, 0, Infinity}]", "Sqrt[Pi]/(2 E^2)");
    /* a, b general: (1/2) Sqrt[Pi/a] E^(-2 Sqrt[a b]) */
    assert_num_eq("Integrate[Exp[-2 x^2 - 3/x^2], {x, 0, Infinity}]", "Sqrt[Pi/2]/2 Exp[-2 Sqrt[6]]");
    /* with a power of x: K_1(2) (was 1/2) */
    assert_num_eq("Integrate[x Exp[-x^2 - 1/x^2], {x, 0, Infinity}]", "BesselK[1, 2]");
    /* every answer agrees with quadrature */
    assert_num_eq("Integrate[Exp[-x^2 - 1/x^2], {x, 0, Infinity}]",
                  "NIntegrate[Exp[-x^2 - 1/x^2], {x, 0, Infinity}, WorkingPrecision -> 45]");
    /* the genuine Gaussian is unaffected */
    assert_true("(" "Integrate[Exp[-x^2], {x, -Infinity, Infinity}]" ") === (" "Sqrt[Pi]" ")");
    assert_true("(" "Integrate[Exp[-x^2 + 2 x], {x, -Infinity, Infinity}]" ") === (" "E Sqrt[Pi]" ")");
}

/* ---- Sum: vanishing leading terms ----------------------------------------- */

void test_sum_vanishing_leading_terms() {
    /* Dobinski: Sum k^n/k! = Bell(n) e (was 0) */
    assert_true("(" "Sum[k^10/k!, {k, 0, Infinity}]" ") === (" "115975 E" ")");
    assert_true("(" "Sum[k^3/k!, {k, 0, Infinity}]" ") === (" "5 E" ")");
    assert_true("(" "Sum[k^2/k!, {k, 0, Infinity}]" ") === (" "2 E" ")");
    /* Touchard: Sum k^3 x^k/k! = e^x (x + 3x^2 + x^3) */
    assert_true("Expand[Sum[k^3 x^k/k!, {k, 0, Infinity}] - E^x (x + 3 x^2 + x^3)] === 0");
    /* the series value, independently */
    assert_num_eq("Sum[k^5 (1/3)^k/k!, {k, 0, Infinity}]", "Sum[N[k^5 (1/3)^k/k!, 50], {k, 0, 200}]");
    /* a term vanishing at the first TWO indices */
    assert_num_eq("Sum[k (k - 1) 2^k/k!, {k, 0, Infinity}]", "4 E^2");
}

/* ---- HypergeometricPFQ: integer-shifted parameter pairs ------------------- */

void test_pfq_shift_pair() {
    assert_true("(" "HypergeometricPFQ[{2}, {1}, z]" ") === (" "E^z (1 + z)" ")");
    assert_true("Expand[HypergeometricPFQ[{2, 2}, {1, 1}, z] - E^z (1 + 3 z + z^2)] === 0");
    assert_true("Expand[HypergeometricPFQ[{3}, {1}, z] - E^z (2 + 4 z + z^2)/2] === 0");
    /* a 2F1 with c+1 over c: (1 - z)^(-a-1) (1 - z + a z/c) */
    assert_true("Simplify[HypergeometricPFQ[{a, c + 1}, {c}, z] - (1 - z)^(-a - 1) (1 - z + a z/c)] === 0");
    /* numerically against the series at a point */
    assert_num_eq("HypergeometricPFQ[{7/2, 3}, {5/2}, 1/5]",
                  "Sum[N[Pochhammer[7/2, k] Pochhammer[3, k]/(Pochhammer[5/2, k] k!) (1/5)^k, 50], {k, 0, 400}]");
    /* no closed form for the remainder: stays */
    assert_true("Head[HypergeometricPFQ[{2, b}, {1, c}, z]] === HypergeometricPFQ");
    /* a lower non-positive integer is a pole, never reduced */
    assert_free_of("HypergeometricPFQ[{1/2}, {1/3}, z]", "Hyp`t");
}

int main() {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(test_solve_list_form_input_intact);
    TEST(test_solve_lambert_product_form);
    TEST(test_solve_lambert_sum_form);
    TEST(test_dsolve_singular_point_conditions);
    TEST(test_integrate_cauchy_schlomilch);
    TEST(test_sum_vanishing_leading_terms);
    TEST(test_pfq_shift_pair);

    printf("All OSP famous-problem regression tests passed.\n");
    return 0;
}
