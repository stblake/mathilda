/* test_meijergreduce.c
 *
 * MeijerGReduce[expr, x] (src/special_functions/meijer_g.c, the rewrite in
 * src/internal/RoachIntegrate.m on the Roach`GForm1 table).
 *
 * The contract checked here is the one that matters: Activate of the result is
 * expr again. Every function is checked NUMERICALLY after Activate at two
 * sample points x > 0 (the half line the G table is written for), with
 * generic values for its parameters -- not by comparing printed forms, which
 * would only pin one spelling of the answer. The structural cases pin the
 * shape: Inactive[MeijerG] factors, constants and powers of x outside, sums
 * mapped, unevaluated when no factor has a form.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    expr_free(p);
    expr_free(r);
    return s;
}

static void check_eq(const char* input, const char* expected) {
    char* s = eval_str(input);
    if (strcmp(s, expected) != 0)
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n", input, expected, s);
    ASSERT_STR_EQ(s, expected);
    free(s);
}

/* f is rewritten (the result contains an Inactive MeijerG and no f-only
 * leftovers are required), and Activate[result] == f at x = 7/10 and
 * x = 23/10 with the parameter values below. */
static void check_reduce(const char* f) {
    char buf[2048];
    snprintf(buf, sizeof buf,
        "With[{r = MeijerGReduce[%s, x], pts = {a -> 13/10, b -> 3/2, n -> 1/3, nu -> 2/5}},"
        " {FreeQ[r, MeijerGReduce], !FreeQ[r, Inactive[MeijerG]],"
        "  Max[Abs[N[(Activate[r] - (%s)) /. pts /. x -> #] & /@ {7/10, 23/10}]] < 10^-9}]",
        f, f);
    check_eq(buf, "{True, True, True}");
}

static void test_numeric_round_trip(void) {
    /* elementary */
    check_reduce("Exp[x]");
    check_reduce("Exp[-a x]");
    check_reduce("Exp[-x^2]");
    check_reduce("E^(2 - 3 x)");
    check_reduce("Sin[x]");
    check_reduce("Cos[b x]");
    check_reduce("Sin[2 x^3]");
    check_reduce("Sin[Sqrt[x]]");
    check_reduce("Sinh[x]");
    check_reduce("Cosh[x]");
    check_reduce("Log[1 + x]");
    check_reduce("Log[1 + a x]/x");
    check_reduce("ArcTan[x]");
    check_reduce("ArcCot[x]");
    check_reduce("1/(1 + x)");
    check_reduce("(1 + x)^a");
    check_reduce("(2 + 3 x^2)^(-1/3)");
    check_reduce("Sqrt[1 - x^2]");
    /* the interval indicator, alone and with a (1 - x)^(a-1) or Log factor */
    check_reduce("UnitStep[1 - x]");
    check_reduce("(1 - x)^(a - 1) UnitStep[1 - x]");
    check_reduce("Log[x] UnitStep[1 - x]");
    /* special functions */
    check_reduce("Erf[x]");
    check_reduce("Erfc[x]");
    check_reduce("ExpIntegralEi[-x]");
    check_reduce("ExpIntegralE[2, x]");
    check_reduce("SinIntegral[x]");
    check_reduce("CosIntegral[x]");
    check_reduce("FresnelS[x]");
    check_reduce("FresnelC[x]");
    check_reduce("Gamma[a, x]");
    check_reduce("BesselJ[0, x]");
    check_reduce("BesselJ[nu, b x]");
    check_reduce("BesselY[0, x]");
    check_reduce("BesselY[nu, x]");
    check_reduce("BesselI[0, x]");
    check_reduce("BesselI[nu, b x]");
    check_reduce("BesselK[n, x]");
    check_reduce("StruveH[0, x]");
    check_reduce("EllipticK[x/3]");
    check_reduce("EllipticE[x/3]");
    check_reduce("BesselJ[0, x]^2");
    /* products, powers, sums, constants */
    check_reduce("x^2 BesselJ[0, x]");
    check_reduce("Sin[x] Cos[x]");
    check_reduce("Sin[x]^2");
    check_reduce("Sin[x] + Cos[x]");
    check_reduce("3 Erf[x]");
    check_reduce("x^a Exp[-x]");
    check_reduce("E^(-x) BesselJ[0, x]");
    check_reduce("x^(3/2) Exp[-x^2/2] FresnelC[x]");
    check_reduce("Exp[-x] + x");
}

static void test_shape(void) {
    /* the four-argument form: argument linear in x, r = 1/q */
    check_eq("MeijerGReduce[Sin[x], x]",
             "Sqrt[Pi] Inactive[MeijerG][{{}, {}}, {{1/2}, {0}}, 1/2 x, 1/2]");
    /* q = 1: the three-argument form */
    check_eq("MeijerGReduce[Log[1 + x], x]",
             "Inactive[MeijerG][{{1, 1}, {}}, {{1}, {0}}, x]");
    check_eq("MeijerGReduce[Exp[-a x], x]",
             "Inactive[MeijerG][{{}, {}}, {{0}, {}}, a x]");
    /* a fractional power keeps the argument c x^q */
    check_eq("MeijerGReduce[Sin[Sqrt[x]], x]",
             "Sqrt[Pi] Inactive[MeijerG][{{}, {}}, {{1/2}, {0}}, 1/4 x]");
    /* constants and powers of x stay outside */
    check_eq("MeijerGReduce[x^2 BesselJ[0, x], x]",
             "x^2 Inactive[MeijerG][{{}, {}}, {{0}, {0}}, 1/2 x, 1/2]");
    check_eq("MeijerGReduce[E^(2 - 3 x), x]",
             "E^2 Inactive[MeijerG][{{}, {}}, {{0}, {}}, 3 x]");
    /* (1 - x)^(a-1) UnitStep[1 - x] is a single G */
    check_eq("MeijerGReduce[(1 - x)^(a - 1) UnitStep[1 - x], x]",
             "Gamma[a] Inactive[MeijerG][{{}, {a}}, {{0}, {}}, x]");
    /* no factor with a form: unevaluated; nothing to convert: itself */
    check_eq("MeijerGReduce[foo[x], x]", "MeijerGReduce[foo[x], x]");
    check_eq("MeijerGReduce[x^3, x]", "x^3");
    check_eq("MeijerGReduce[5, x]", "5");
    /* the unconvertible factor is kept in a product */
    check_eq("FreeQ[MeijerGReduce[foo[x] Sin[x], x], foo[x]]", "False");
    /* a non-symbol variable or the wrong arity: unevaluated */
    check_eq("MeijerGReduce[Sin[x], 3]", "MeijerGReduce[Sin[x], 3]");
    check_eq("MeijerGReduce[Sin[x]]", "MeijerGReduce[Sin[x]]");
    check_eq("Attributes[MeijerGReduce]", "{Protected}");
}

int main(void) {
    symtab_init();
    core_init();
    test_mute_messages();
    TEST(test_numeric_round_trip);
    TEST(test_shape);
    printf("All MeijerGReduce tests passed.\n");
    return 0;
}
