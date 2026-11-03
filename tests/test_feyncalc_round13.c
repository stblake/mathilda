/* test_feyncalc_round13.c
 *
 * Kernel fixes for FeynCalc's Feynman-parameter prefactors (v0.6146):
 *   - Power writes an exact rational base under a symbolic exponent in the
 *     reference's normal form: (1/4)^e = 4^(-e), 4^(-D/2) = 2^(-D),
 *     12^(D/2) = 2^D 3^(D/2), (-1)^(x/2) = I^x;
 *   - Simplify moves the integer content of a symbolic exponent into an exact
 *     base when that is no more complex (4^(-2 + 2 e) -> 16^(-1 + e));
 *   - the special-function engine no longer trades a Gamma ratio for
 *     Legendre's duplication when the score rises and no family disappears;
 *   - Together splits the integer part off B^(k + s) for a base whose normal
 *     form is a negated associate ((b - a)^(-1 + n) = -(-a + b)^n/(a - b)),
 *     and keeps a sum's symbolic power an opaque factor;
 *   - InputForm puts a power with a negative symbolic exponent in the
 *     denominator (c/x^n, 3^(x/2)/2^x).
 * Expected values computed by the reference kernel (InputForm strings).
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

static void eq_if(const char* e, const char* expected) {
    char buf[2048], quoted[2048];
    snprintf(buf, sizeof buf, "ToString[%s, InputForm]", e);
    snprintf(quoted, sizeof quoted, "\"%s\"", expected);   /* a String prints quoted */
    assert_eval_eq(buf, quoted, 0);
}

static void test_power_rational_base_symbolic_exponent(void) {
    eq_if("4^(-D/2)", "2^(-D)");
    eq_if("12^(D/2)", "2^D*3^(D/2)");
    eq_if("(1/4)^(2 - 2 e)", "4^(-2 + 2*e)");
    eq_if("4^(1 + x/2)", "4^(1 + x/2)");
    eq_if("6^(x/2)", "6^(x/2)");
    eq_if("8^(2 x/3)", "2^(2*x)");
    eq_if("64^(x/4)", "2^((3*x)/2)");
    eq_if("1000^(x/2)", "10^((3*x)/2)");
    eq_if("72^(x/2)", "2^((3*x)/2)*3^x");
    eq_if("(-3)^(x/2)", "I^x*3^(x/2)");
    eq_if("(-1)^(x/2)", "I^x");
    eq_if("(-8)^(x/3)", "(-1)^(x/3)*2^x");
    eq_if("(-6)^(x/3)", "(-6)^(x/3)");
    eq_if("(2/3)^(-x)", "(3/2)^x");
    eq_if("(4/3)^(-x)", "(4/3)^(-x)");
    eq_if("(2^70)^(x/2)", "34359738368^x");
    eq_if("2^(D/2) 4^(-D/2)", "2^(-1/2*D)");
}

static void test_simplify_exponent_content(void) {
    eq_if("Simplify[4^(-2 + 2 e)]", "16^(-1 + e)");
    eq_if("Simplify[2^(4 x)]", "16^x");
    eq_if("Simplify[2^(2 x) 3^(2 x)]", "36^x");
    eq_if("Simplify[2^(2 x + 2)]", "4^(1 + x)");
    eq_if("Simplify[4^(3 e - 3)]", "64^(-1 + e)");
}

static void test_simplify_gamma_ratio_kept(void) {
    eq_if("Simplify[Gamma[-2 + 2 e]/Gamma[-2 + e]]", "Gamma[-2 + 2*e]/Gamma[-2 + e]");
    eq_if("Simplify[4^(-2 + 2 e) Gamma[-2 + 2 e]/Gamma[-2 + e]]", "(16^(-1 + e)*Gamma[-2 + 2*e])/Gamma[-2 + e]");
    eq_if("Simplify[((-1)^e Gamma[e + (-4 + 2 e)/2])/Gamma[-2 + e]]", "((-1)^e*Gamma[-2 + 2*e])/Gamma[-2 + e]");
}

static void test_together_signed_base_integer_part(void) {
    eq_if("Together[(b - a)^(-1 + n)]", "-((-a + b)^n/(a - b))");
    eq_if("Together[(a + b)^(-1 + n)]", "(a + b)^(-1 + n)");
    eq_if("Together[(1 - x)^(-1 + n)]", "-((1 - x)^n/(-1 + x))");
    eq_if("Together[(b - a)^(-1 + n) + 1]", "(a - b - (-a + b)^n)/(a - b)");
    eq_if("Together[(b - a)^(n - 2)]", "(-a + b)^n/(a - b)^2");
    eq_if("Together[1/((a - b) (c + d)^n)]", "1/((a - b)*(c + d)^n)");
    eq_if("Together[(2 - 2 x)^(-1 + n)]", "-1/2*(2 - 2*x)^n/(-1 + x)");
    eq_if("Together[(b - a)^(-1 + n) (c - d)^(-1 + m)]", "-(((-a + b)^n*(c - d)^(-1 + m))/(a - b))");
}

static void test_inputform_negative_symbolic_exponent(void) {
    eq_if("c x^(-n)", "c/x^n");
    eq_if("3^(x/2) 2^(-x)", "3^(x/2)/2^x");
    eq_if("a E^(-x)", "a/E^x");
    eq_if("-x^(-n)", "-x^(-n)");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_power_rational_base_symbolic_exponent);
    TEST(test_simplify_exponent_content);
    TEST(test_simplify_gamma_ratio_kept);
    TEST(test_together_signed_base_integer_part);
    TEST(test_inputform_negative_symbolic_exponent);
    printf("All FeynCalc round-13 tests passed.\n");
    return 0;
}
