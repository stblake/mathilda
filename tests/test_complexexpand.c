#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

/*
 * Tests for ComplexExpand.
 *
 * Two verification standards are used (see the ComplexExpand plan):
 *   - assert_eval_eq  : exact-string, pinned to Mathilda's canonical
 *                       (Orderless-sorted) output.  Used where the result is
 *                       a clean closed form.
 *   - assert_cx_zero  : surface-form-independent correctness.  Substitutes
 *                       real x, y (and complex z) and checks that
 *                       ComplexExpand[e] - e chops to 0 numerically.  Used for
 *                       the messy inverse-function / Log / Conjugate-target
 *                       results whose printed form differs from (but is equal
 *                       to) the Wolfram reference.
 */

/* Numeric zero-difference: (diff_expr) must vanish at x=3/10, y=7/10,
 * z=3/10+7/10 I. */
static void assert_cx_zero(const char* diff_expr) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
             "Chop[N[(%s) /. {x -> 3/10, y -> 7/10, z -> 3/10 + 7/10 I}, 25]]",
             diff_expr);
    assert_eval_eq(buf, "0", 0);
}

/* ------------------------------------------------------------------ */
/* Circular and hyperbolic functions                                  */
/* ------------------------------------------------------------------ */
void test_cx_trig(void) {
    assert_eval_eq("ComplexExpand[Sin[x + I y]]",
                   "Sin[x] Cosh[y] + I Cos[x] Sinh[y]", 0);
    assert_eval_eq("ComplexExpand[Cos[x + I y]]",
                   "Cos[x] Cosh[y] - I Sin[x] Sinh[y]", 0);
    assert_eval_eq("ComplexExpand[Sinh[x + I y]]",
                   "I Cosh[x] Sin[y] + Sinh[x] Cos[y]", 0);
    assert_eval_eq("ComplexExpand[Tan[x + I y]]",
                   "Sin[2 x]/(Cos[2 x] + Cosh[2 y]) + I Sinh[2 y]/(Cos[2 x] + Cosh[2 y])", 0);
    /* A complex variable given via the second argument. */
    assert_eval_eq("ComplexExpand[Cosh[x + I y] + Tanh[z], {z}]",
                   "Cosh[x] Cos[y] + I Sinh[x] Sin[y] + I Sin[2 Im[z]]/(Cos[2 Im[z]] + Cosh[2 Re[z]]) + Sinh[2 Re[z]]/(Cos[2 Im[z]] + Cosh[2 Re[z]])", 0);
    /* Single complex variable given bare (not in a list). */
    assert_eval_eq("ComplexExpand[Sin[x], x]",
                   "Cosh[Im[x]] Sin[Re[x]] + I Cos[Re[x]] Sinh[Im[x]]", 0);
    /* Remaining circular/hyperbolic + reciprocals: correctness only. */
    assert_cx_zero("ComplexExpand[Cot[x + I y]] - Cot[x + I y]");
    assert_cx_zero("ComplexExpand[Sec[x + I y]] - Sec[x + I y]");
    assert_cx_zero("ComplexExpand[Csc[x + I y]] - Csc[x + I y]");
    assert_cx_zero("ComplexExpand[Coth[x + I y]] - Coth[x + I y]");
    assert_cx_zero("ComplexExpand[Sech[x + I y]] - Sech[x + I y]");
    assert_cx_zero("ComplexExpand[Csch[x + I y]] - Csch[x + I y]");
}

/* ------------------------------------------------------------------ */
/* Polynomials (complex variables via second argument)                */
/* ------------------------------------------------------------------ */
void test_cx_polynomial(void) {
    assert_eval_eq("ComplexExpand[Re[z^2], {z}]", "-Im[z]^2 + Re[z]^2", 0);
    assert_eval_eq("ComplexExpand[Re[z^5 - 2 z^3 - z + 1], z]",
                   "1 + Re[z]^5 - 2 Re[z]^3 - Re[z] - 10 Im[z]^2 Re[z]^3 + 5 Im[z]^4 Re[z] + 6 Im[z]^2 Re[z]", 0);
    assert_eval_eq("ComplexExpand[Re[2 z^3 - z + 1], z]",
                   "1 - Re[z] + 2 Re[z]^3 - 6 Im[z]^2 Re[z]", 0);
    assert_eval_eq("ComplexExpand[a + x^2, {x}]",
                   "a - Im[x]^2 + Re[x]^2 + (2*I) Im[x] Re[x]", 0);
    /* Negative integer power: correctness. */
    assert_cx_zero("ComplexExpand[z^-2, {z}] - z^-2");
}

/* ------------------------------------------------------------------ */
/* Exponential, logarithm, powers of numbers                          */
/* ------------------------------------------------------------------ */
void test_cx_exp_log(void) {
    assert_eval_eq("ComplexExpand[3^(I x)]", "Cos[Log[3] x] + I Sin[Log[3] x]", 0);
    assert_eval_eq("ComplexExpand[(-1)^(4 I/5)]", "E^(-4/5 Pi)", 0);
    assert_eval_eq("ComplexExpand[(-3)^(4 I/5)]",
                   "Cos[4/5 Log[3]] E^(-4/5 Pi) + I E^(-4/5 Pi) Sin[4/5 Log[3]]", 0);
    /* (-1)^(4/5): exact radicals from Mathilda's trig tables. */
    assert_eval_eq("ComplexExpand[(-1)^(4/5)]",
                   "-1/4 - 1/4 Sqrt[5] + (1/4*I) Sqrt[10 - 2 Sqrt[5]]", 0);
    /* ReIm of a composed exponential (documented example). */
    assert_eval_eq("ComplexExpand[ReIm[Cos[I^x]]]",
                   "{Cos[Cos[1/2 Pi x]] Cosh[Sin[1/2 Pi x]], -Sin[Cos[1/2 Pi x]] Sinh[Sin[1/2 Pi x]]}", 0);
    /* Log and a nested composition: correctness. */
    assert_cx_zero("ComplexExpand[Log[x + I y]] - Log[x + I y]");
    assert_cx_zero("ComplexExpand[(I - 1)^(2/5)] - (I - 1)^(2/5)");
    assert_cx_zero("ComplexExpand[Re[Log[Sin[Exp[x + I y]^2]]]] - Re[Log[Sin[Exp[x + I y]^2]]]");
}

/* ------------------------------------------------------------------ */
/* TargetFunctions                                                    */
/* ------------------------------------------------------------------ */
void test_cx_targetfunctions(void) {
    assert_eval_eq("ComplexExpand[Re[Tan[z]], z]",
                   "Sin[2 Re[z]]/(Cosh[2 Im[z]] + Cos[2 Re[z]])", 0);
    assert_eval_eq("ComplexExpand[Re[z^2], {z}, TargetFunctions -> {Abs, Arg}]",
                   "Abs[z]^2 Cos[Arg[z]]^2 - Abs[z]^2 Sin[Arg[z]]^2", 0);
    assert_eval_eq("ComplexExpand[Re[Tan[z]], z, TargetFunctions -> {Abs, Arg}]",
                   "Sin[2 Abs[z] Cos[Arg[z]]]/(Cos[2 Abs[z] Cos[Arg[z]]] + Cosh[2 Abs[z] Sin[Arg[z]]])", 0);
    assert_eval_eq("ComplexExpand[Re[z^2], {z}, TargetFunctions -> Conjugate]",
                   "1/2 (z^2 + Conjugate[z]^2)", 0);
    /* Conjugate target on transcendental heads: correctness. */
    assert_cx_zero("ComplexExpand[Re[Sin[z]], z, TargetFunctions -> Conjugate] - Re[Sin[z]]");
    assert_cx_zero("ComplexExpand[Re[Tan[z]], z, TargetFunctions -> Conjugate] - Re[Tan[z]]");
    /* Abs/Arg target correctness. */
    assert_cx_zero("ComplexExpand[Abs[z^2], {z}, TargetFunctions -> {Abs, Arg}] - Abs[z^2]");
}

/* ------------------------------------------------------------------ */
/* Wrapper heads: Re, Im, Abs, Arg, Conjugate, Sign, ReIm             */
/* ------------------------------------------------------------------ */
void test_cx_wrappers(void) {
    assert_eval_eq("ComplexExpand[Conjugate[z], {z}]", "-I Im[z] + Re[z]", 0);
    assert_eval_eq("ComplexExpand[Abs[z], {z}]", "Sqrt[Im[z]^2 + Re[z]^2]", 0);
    /* Sign of a complex variable -> z/Abs[z]. */
    assert_eval_eq("ComplexExpand[Sign[x], x]",
                   "Re[x]/Sqrt[Im[x]^2 + Re[x]^2] + I Im[x]/Sqrt[Im[x]^2 + Re[x]^2]", 0);
    /* Real variables pass through the real-valued wrappers. */
    assert_eval_eq("ComplexExpand[Re[x + y]]", "x + y", 0);
    assert_eval_eq("ComplexExpand[Im[x + y]]", "0", 0);
    assert_eval_eq("ComplexExpand[Conjugate[x + y]]", "x + y", 0);
}

/* ------------------------------------------------------------------ */
/* Inverse circular / hyperbolic functions                            */
/* ------------------------------------------------------------------ */
void test_cx_inverse(void) {
    /* ArcTan matches the Wolfram oracle exactly. */
    assert_eval_eq("ComplexExpand[{Re[ArcTan[x + I y]], Im[ArcTan[x + I y]]}]",
                   "{-1/2 Arg[1 - I (x + I y)] + 1/2 Arg[1 + I (x + I y)], -1/4 Log[1 + x^2 - 2 y + y^2] + 1/4 Log[1 + x^2 + 2 y + y^2]}", 0);
    /* Documented composite example. */
    assert_eval_eq("ComplexExpand[Re[ArcCot[x + I y]] + Im[ArcSinh[x - I y]]]",
                   "Arg[x - I y + Sqrt[1 + (x - I y)^2]] - 1/2 Arg[1 - I/(x + I y)] + 1/2 Arg[1 + I/(x + I y)]", 0);
    /* Every inverse head: correctness against the true function value. */
    assert_cx_zero("ComplexExpand[ArcSin[x + I y]] - ArcSin[x + I y]");
    assert_cx_zero("ComplexExpand[ArcCos[x + I y]] - ArcCos[x + I y]");
    assert_cx_zero("ComplexExpand[ArcTan[x + I y]] - ArcTan[x + I y]");
    assert_cx_zero("ComplexExpand[ArcCot[x + I y]] - ArcCot[x + I y]");
    assert_cx_zero("ComplexExpand[ArcSec[x + I y]] - ArcSec[x + I y]");
    assert_cx_zero("ComplexExpand[ArcCsc[x + I y]] - ArcCsc[x + I y]");
    assert_cx_zero("ComplexExpand[ArcSinh[x + I y]] - ArcSinh[x + I y]");
    assert_cx_zero("ComplexExpand[ArcCosh[x + I y]] - ArcCosh[x + I y]");
    assert_cx_zero("ComplexExpand[ArcTanh[x + I y]] - ArcTanh[x + I y]");
    assert_cx_zero("ComplexExpand[ArcCoth[x + I y]] - ArcCoth[x + I y]");
    assert_cx_zero("ComplexExpand[ArcSech[x + I y]] - ArcSech[x + I y]");
    assert_cx_zero("ComplexExpand[ArcCsch[x + I y]] - ArcCsch[x + I y]");
}

/* ------------------------------------------------------------------ */
/* Threading and identities                                           */
/* ------------------------------------------------------------------ */
void test_cx_threading(void) {
    assert_eval_eq("ComplexExpand[{Sin[x], Cos[x + I y]}]",
                   "{Sin[x], Cos[x] Cosh[y] - I Sin[x] Sinh[y]}", 0);
    assert_eval_eq("ComplexExpand[x < y && z == 1]", "x < y && z == 1", 0);
    /* Verified complex identities collapse to True. */
    assert_eval_eq("ComplexExpand[Re[z] == (z + Conjugate[z])/2, z]", "True", 0);
    assert_eval_eq("ComplexExpand[Exp[I z] == Cos[z] + I Sin[z]]", "True", 0);
}

/* ------------------------------------------------------------------ */
/* Multiple complex variables                                         */
/* ------------------------------------------------------------------ */
void test_cx_multivar(void) {
    assert_eval_eq("ComplexExpand[Sin[x] Exp[y], {x, y}]",
                   "I Cosh[Im[x]] Sin[Re[x]] E^Re[y] Sin[Im[y]] + Cosh[Im[x]] Sin[Re[x]] Cos[Im[y]] E^Re[y] + I Cos[Re[x]] Sinh[Im[x]] Cos[Im[y]] E^Re[y] - Cos[Re[x]] Sinh[Im[x]] E^Re[y] Sin[Im[y]]", 0);
}

/* ------------------------------------------------------------------ */
/* Argument-count handling                                            */
/* ------------------------------------------------------------------ */
void test_cx_argcount(void) {
    /* Wrong arg counts leave the call unevaluated (a General::argct
     * message is emitted on stderr). */
    assert_eval_eq("ComplexExpand[]", "ComplexExpand[]", 0);
    assert_eval_eq("ComplexExpand[x, 2, 3, 4]", "ComplexExpand[x, 2, 3, 4]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_cx_trig);
    TEST(test_cx_polynomial);
    TEST(test_cx_exp_log);
    TEST(test_cx_targetfunctions);
    TEST(test_cx_wrappers);
    TEST(test_cx_inverse);
    TEST(test_cx_threading);
    TEST(test_cx_multivar);
    TEST(test_cx_argcount);

    printf("All ComplexExpand tests passed!\n");
    return 0;
}
