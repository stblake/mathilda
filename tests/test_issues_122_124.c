/* test_issues_122_124.c
 *
 * GitHub #122: Simplify takes a numeric content out of a sum when the
 * factored form is no more complex -- 2 + 2 x -> 2 (1 + x) -- as it already
 * did for a symbolic one (a + a x -> a (1 + x)); a content whose lift costs
 * more stays (3 + 6 x), and an approximate sum is left as written.
 * GitHub #123: OptionsPattern[] under an Orderless head.  The LHS sorts to
 * h[OptionsPattern[], a_] and the arguments to h[1, o -> v], so the options
 * no longer lead; the slot now takes them wherever they sit.
 * GitHub #124: Factor (and the other heads that round-trip inexact input
 * through exact arithmetic) keeps the input's arbitrary precision, and an
 * approximate factorisation is written monic, c (r + x).
 * Expected forms are the reference kernel's (InputForm strings; Mathilda
 * prints a machine 2. as 2.0).
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void eq_if(const char* e, const char* expected) {
    char buf[2048], quoted[2048];
    snprintf(buf, sizeof buf, "ToString[%s, InputForm]", e);
    snprintf(quoted, sizeof quoted, "\"%s\"", expected);   /* a String prints quoted */
    assert_eval_eq(buf, quoted, 0);
}

static void test_simplify_numeric_content(void) {
    /* #122, verbatim */
    eq_if("Simplify[2 + 2 x]", "2*(1 + x)");
    eq_if("Simplify[2 + 2 Sin[x]]", "2*(1 + Sin[x])");
    eq_if("Simplify[6 + 6 x^2]", "6*(1 + x^2)");
    eq_if("Simplify[2 Sqrt[2] + 2 x]", "2*(Sqrt[2] + x)");
    eq_if("Simplify[2 x + 4 y^2]", "2*(x + 2*y^2)");
    eq_if("Simplify[-2 + 2 x]", "2*(-1 + x)");
    eq_if("Simplify[2 + 2 x^3]", "2*(1 + x^3)");
    eq_if("Simplify[2 + 2 E^x]", "2*(1 + E^x)");
    eq_if("Simplify[1/2 + x/3]", "(3 + 2*x)/6");
    /* already factored before, and still */
    eq_if("Simplify[a + a x]", "a*(1 + x)");
    eq_if("Simplify[-2 - 2 x]", "-2*(1 + x)");
    eq_if("Simplify[(2 + 2 x)/y]", "(2*(1 + x))/y");
    /* the lift costs more: kept */
    eq_if("Simplify[3 + 6 x]", "3 + 6*x");
    eq_if("Simplify[4 + 6 x]", "4 + 6*x");
    /* an approximate sum is left as written */
    eq_if("Simplify[2. + 2. x]", "2.0 + 2.0*x");
    eq_if("Simplify[-2. + 2. x^2]", "-2.0 + 2.0*x^2");
}

static void test_orderless_optionspattern(void) {
    /* #123, verbatim */
    assert_eval_eq("SetAttributes[fun, Orderless]; Options[fun] = {Explicit -> False}; "
                   "fun[a_, b_, c_, OptionsPattern[]] := (MATCH) /; OptionValue[Explicit]; "
                   "fun[a, b, c, Explicit -> True]", "MATCH", 0);
    assert_eval_eq("fun[a, b, c]", "fun[a, b, c]", 0);
    assert_eval_eq("fun[Explicit -> True, c, b, a]", "MATCH", 0);
    /* no Condition, and OptionValue in the body */
    assert_eval_eq("SetAttributes[oh1, Orderless]; Options[oh1] = {e1 -> False}; "
                   "oh1[a_, OptionsPattern[]] := OptionValue[e1]; oh1[1, e1 -> True]", "True", 0);
    assert_eval_eq("oh1[1]", "False", 0);
    assert_eval_eq("oh1[1, {e1 -> 7}]", "7", 0);
    /* an option-shaped argument can still go to a blank */
    assert_eval_eq("SetAttributes[oh2, Orderless]; oh2[a_, OptionsPattern[]] := {a}; "
                   "oh2[u -> 1]", "{u -> 1}", 0);
    /* the named form, which already worked */
    assert_eval_eq("SetAttributes[oh3, Orderless]; Options[oh3] = {e1 -> False}; "
                   "oh3[a_, o:OptionsPattern[]] := {a, {o}}; oh3[1, e1 -> True]",
                   "{1, {e1 -> True}}", 0);
}

static void test_factor_precision(void) {
    /* #124, verbatim: the reference's 49.52; Mathilda tracks no error growth */
    assert_eval_eq("Precision @ D[#, x] & @ Factor[N[Pi, 50] x + 1] > 49", "True", 0);
    /* the digits are right, not just the precision label */
    assert_eval_eq("Abs[Coefficient[Expand[Factor[N[Pi, 50] x + 1]], x] - N[Pi, 50]] < 10^-45",
                   "True", 0);
    assert_eval_eq("Abs[First[Factor[N[Pi, 50] x + 1]] - N[Pi, 50]] < 10^-45", "True", 0);
    assert_eval_eq("Abs[(Factor[N[Pi, 50] x^2 - 1] /. x -> 2) - (4 N[Pi, 50] - 1)] < 10^-45",
                   "True", 0);
    /* the siblings sharing the round trip keep it too */
    assert_eval_eq("Precision[Together[N[Pi, 40] x + 1/x] /. x -> 3] > 39", "True", 0);
    assert_eval_eq("Precision[Simplify[N[Pi, 40] x + x] /. x -> 3] > 39", "True", 0);
    /* an approximate factorisation is written monic: c (r + x) */
    eq_if("Factor[2.5 x + 1]", "2.5*(0.4 + x)");
    eq_if("Factor[3.5 x^2 + 2.1 x]", "3.5*x*(0.6 + x)");
    eq_if("Factor[2.5 x + 3.5 y + 1]", "2.5*(0.4 + x + 1.4*y)");
    eq_if("Factor[(2.5 x + 1)^2]", "6.25*(0.4 + x)^2");
    /* a machine number still round-trips at machine precision */
    assert_eval_eq("Precision[Factor[2.5 x + 1]]", "MachinePrecision", 0);
    assert_eval_eq("Precision[Factor[N[Pi, 50] x + 2.5]]", "MachinePrecision", 0);
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();
    TEST(test_simplify_numeric_content);
    TEST(test_orderless_optionspattern);
    TEST(test_factor_precision);
    printf("All issue #122/#123/#124 tests passed.\n");
    return 0;
}
