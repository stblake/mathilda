/* test_feyncalc_round14.c
 *
 * Kernel fixes found by FeynCalc's last failing unit tests (round 14, v0.6168+):
 *   - Factor of a product with a negative real coefficient and a sum with
 *     Gaussian coefficients decides the sign on the expansion, as the
 *     reference does: Factor[-2 (I x - y)] = 2 (-I x + y), the same as
 *     Factor[-2 I x + 2 y] (DiracTrace in the BMHV scheme left two equal
 *     traces that did not cancel).
 *   - InputForm keeps the parentheses of a held nested product or sum
 *     (HoldForm[a*b*(c*d)], HoldForm[a - (b - c)] -- the latter used to
 *     print as a - b - c, which reads back as a different value) and groups a
 *     complex coefficient with its numerator ((I*x)/y); a lone quotient
 *     numerator stays bare (HoldForm[a/b/c]).
 *   - Block localises attributes and built-in rules: Block[{f}, ...] sees no
 *     attributes and no C builtin (FeynCalc's TeX preview runs //. on inert
 *     Plus and Times inside Block[{Plus, Times}, ...]).
 *   - Replace / ReplaceAt / ReplaceList match whole expressions only: the
 *     Flat-subset rewrite (1 + x + y /. x + y -> z) is ReplaceAll's and
 *     ReplaceRepeated's alone (Replace[a + b + c, a + b -> x] was c + x).
 * Expected values computed by the reference kernel.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* lhs and rhs evaluate to the same expression */
static void same(const char* lhs, const char* rhs) {
    char buf[4096];
    snprintf(buf, sizeof buf, "(%s) === (%s)", lhs, rhs);
    assert_eval_eq(buf, "True", 0);
}

static void test_factor_gaussian_sign_on_expansion(void) {
    same("Factor[-2 (I x - y)]", "2*(-I*x + y)");
    same("Factor[-2 I x + 2 y]", "2*(-I*x + y)");
    same("Factor[-4 (I x - y) z]", "4*(-I*x + y)*z");
    same("Factor[-2 (I x - y) (a + b)]", "2*(a + b)*(-I*x + y)");
    same("Factor[-2 (I a b - c d) e]", "2*(-I*a*b + c*d)*e");
    /* the two BMHV traces of DiracTrace-ID46, factored, now cancel */
    same("Factor[-2 (I x - y)] - Factor[2 (-I x + y)]", "0");
    /* unchanged: a real or positive coefficient, a power of a sum */
    same("Factor[-2 (x - I y)]", "-2*(x - I*y)");
    same("Factor[(2 + 2 I) (x - y)]", "(2 + 2*I)*(x - y)");
    same("Factor[-I (x - y)]", "-I*(x - y)");
}

static void eq_if(const char* e, const char* expected) {
    char buf[2048], quoted[2048];
    snprintf(buf, sizeof buf, "ToString[%s, InputForm]", e);
    snprintf(quoted, sizeof quoted, "\"%s\"", expected);   /* a String prints quoted */
    assert_eval_eq(buf, quoted, 0);
}

static void test_inputform_nested_parens(void) {
    eq_if("HoldForm[(25*NH)/27]", "HoldForm[(25*NH)/27]");
    eq_if("HoldForm[((c + a)*ep)/M]", "HoldForm[((c + a)*ep)/M]");
    eq_if("HoldForm[a + b + (c + d)]", "HoldForm[a + b + (c + d)]");
    eq_if("HoldForm[a*b*(c*d)]", "HoldForm[a*b*(c*d)]");
    eq_if("HoldForm[a*(b*c)/d]", "HoldForm[a*((b*c)/d)]");
    eq_if("HoldForm[a - (b - c)]", "HoldForm[a - (b - c)]");
    eq_if("HoldForm[-(a b)]", "HoldForm[-(a*b)]");
    eq_if("HoldForm[x (y/z)]", "HoldForm[x*(y/z)]");
    eq_if("HoldForm[(a/b)/c]", "HoldForm[a/b/c]");
    eq_if("I x/y", "(I*x)/y");
    eq_if("(-I) x/y", "((-I)*x)/y");
    eq_if("(2 I) x/y", "((2*I)*x)/y");
    eq_if("I/y", "I/y");
    eq_if("25 NH/27", "(25*NH)/27");
    /* the held text reads back as the same value */
    same("ReleaseHold[ToExpression[ToString[HoldForm[a - (b - c)], InputForm]]]", "a - b + c");
}

static void test_block_localises_attributes_and_builtins(void) {
    same("Block[{Times}, Hold[Evaluate[Attributes[Times]]]]", "Hold[{}]");
    same("SetAttributes[blkf, Orderless]; Block[{blkf}, Hold[Evaluate[{Attributes[blkf], blkf[b, a]}]]]",
         "Hold[{{}, blkf[b, a]}]");
    same("{Attributes[blkf], blkf[b, a]}", "{{Orderless}, blkf[a, b]}");
    same("Block[{Sin}, Hold[Evaluate[Sin[{0, 1}]]]]", "Hold[Sin[{0, 1}]]");
    same("Sin[{0, 1}]", "{0, Sin[1]}");
    same("Attributes[Times]", "{Flat, Listable, NumericFunction, OneIdentity, Orderless, Protected}");
    same("Block[{$RecursionLimit = 50}, $RecursionLimit]", "50");
    same("Block[{x = 2}, x^2]", "4");
}

static void test_replace_whole_expression_only(void) {
    /* Replace, ReplaceAt and ReplaceList match whole expressions; only
     * ReplaceAll / ReplaceRepeated rewrite a subset of a Flat head */
    same("Replace[a + b + c, a + b -> x]", "a + b + c");
    same("Replace[a b c, a b -> x]", "a b c");
    same("Replace[f[a + b + c], a + b -> x, {1}]", "f[a + b + c]");
    same("Replace[f[a + b + c], a + b -> x, All]", "f[a + b + c]");
    same("ReplaceAt[f[a + b + c], a + b -> x, 1]", "f[a + b + c]");
    same("Replace[a + b + c, {a + b -> x, c -> 1}]", "a + b + c");
    same("ReplaceList[a + b + c, a + b -> x]", "{}");
    same("(a + b + c) /. a + b -> x", "c + x");
    same("(a + b + c) //. a + b -> x", "c + x");
    same("Replace[a + b + c, a + b + c -> x]", "x");
    same("Replace[a + b + c, x_ + b -> {x}]", "{a + c}");
}

static void test_numloop_seed_and_collect_head(void) {
    /* a loop does not numericalise a non-numeric value its variable holds
     * (user rules on machine numbers fired: FeynCalc TID-ID50) */
    same("nlhit = 0; nlh[x_Real, y_] := (nlhit++; nlhh); nlres = nlh[0, k - p]; "
         "Do[nlres = nlfoo[nlres], {i, 1}]; {nlhit, nlres}",
         "{0, nlfoo[nlh[0, k - p]]}");
    /* Collect's head runs at the caller's level: Together answers factored */
    same("Collect[c1 X + c2, {c1, c2}, Together[(a b)/(x0 (1 + x0))] &]",
         "(a b c1)/(x0 (1 + x0)) + (a b c2)/(x0 (1 + x0))");
    same("Collect[(a b)/(x0 + x0^2) c1, c1, Together]", "(a b c1)/(x0 (1 + x0))");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_factor_gaussian_sign_on_expansion);
    TEST(test_inputform_nested_parens);
    TEST(test_block_localises_attributes_and_builtins);
    TEST(test_replace_whole_expression_only);
    TEST(test_numloop_seed_and_collect_head);
    printf("All FeynCalc round-14 tests passed.\n");
    return 0;
}
