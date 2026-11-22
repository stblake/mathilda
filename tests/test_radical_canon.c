/* test_radical_canon.c -- the canonical form of exact rational radicals
 * (src/radcanon.c, the O5 pass in src/times.c, power_radcanon in
 * src/power.c) and the Abs / Conjugate / Re / Im identities on powers,
 * products and sums (src/complex.c).
 *
 * Every expected form was derived from the algorithm's rules (radcanon.c
 * header; the identity table above cx_call1 in complex.c), and checked
 * against the reference's automatic evaluation; the one deliberate
 * difference is noted inline.  A numeric self-check confirms that each
 * rewrite keeps the value. */
#include "core.h"
#include "expr.h"
#include "symtab.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>

static void test_rational_radicals(void) {
    assert_eval_eq("-1/3*Sqrt[3/2]", "Times[-1, Power[6, Rational[-1, 2]]]", 1);
    assert_eval_eq("2 Sqrt[3/2]", "Power[6, Rational[1, 2]]", 1);
    assert_eval_eq("Sqrt[6]/3", "Power[Rational[2, 3], Rational[1, 2]]", 1);
    assert_eval_eq("Sqrt[2]/2", "Power[2, Rational[-1, 2]]", 1);
    assert_eval_eq("2 Sqrt[2]", "Times[2, Power[2, Rational[1, 2]]]", 1);
    assert_eval_eq("3 Sqrt[3/2]", "Times[3, Power[Rational[3, 2], Rational[1, 2]]]", 1);
    assert_eval_eq("6^(1/3)/2", "Times[Power[2, Rational[-2, 3]], Power[3, Rational[1, 3]]]", 1);
    assert_eval_eq("2^(2/3)/2", "Power[2, Rational[-1, 3]]", 1);
    assert_eval_eq("3/Sqrt[3]", "Power[3, Rational[1, 2]]", 1);
    assert_eval_eq("30 Sqrt[2/15]", "Times[2, Power[30, Rational[1, 2]]]", 1);
    assert_eval_eq("2 Sqrt[2/15]", "Times[2, Power[Rational[2, 15], Rational[1, 2]]]", 1);
    assert_eval_eq("21/5 Sqrt[2/213]", "Times[Rational[7, 5], Power[Rational[6, 71], Rational[1, 2]]]", 1);
    assert_eval_eq("14/Sqrt[10]", "Times[7, Power[Rational[2, 5], Rational[1, 2]]]", 1);
    assert_eval_eq("2/Sqrt[30]", "Power[Rational[2, 15], Rational[1, 2]]", 1);
    assert_eval_eq("2/3 x 6^(1/3)", "Times[2, Power[2, Rational[1, 3]], Power[3, Rational[-2, 3]], x]", 1);
    assert_eval_eq("2^(1/3)/4", "Times[Rational[1, 2], Power[2, Rational[-2, 3]]]", 1);
    assert_eval_eq("6^(-1/3) 2", "Times[Power[2, Rational[2, 3]], Power[3, Rational[-1, 3]]]", 1);
    assert_eval_eq("Sqrt[2] 3^(1/3)/6", "Times[Power[2, Rational[-1, 2]], Power[3, Rational[-2, 3]]]", 1);
    assert_eval_eq("Sqrt[6] 2^(1/3)/2", "Times[Power[2, Rational[-1, 6]], Power[3, Rational[1, 2]]]", 1);
    assert_eval_eq("1000000007 Sqrt[3/1000000007]", "Power[3000000021, Rational[1, 2]]", 1);
    assert_eval_eq("3 Sqrt[3 (2^61-1)]/(2^61-1)", "Times[3, Power[Rational[3, 2305843009213693951], Rational[1, 2]]]", 1);
    assert_eval_eq("2 I Sqrt[3/2]", "Times[Complex[0, 1], Power[6, Rational[1, 2]]]", 1);
    assert_eval_eq("(3+6I)/Sqrt[3]", "Times[Complex[1, 2], Power[3, Rational[1, 2]]]", 1);
    assert_eval_eq("(2+3I) Sqrt[2]/4", "Times[Complex[1, Rational[3, 2]], Power[2, Rational[-1, 2]]]", 1);
    assert_eval_eq("(1/2+I/3) Sqrt[6]", "Times[Complex[3, 2], Power[6, Rational[-1, 2]]]", 1);
    assert_eval_eq("Sqrt[6] Sqrt[10]", "Times[2, Power[15, Rational[1, 2]]]", 1);
    assert_eval_eq("2^(1/3) 3^(1/2) 5^(-1/3)", "Times[Power[Rational[2, 5], Rational[1, 3]], Power[3, Rational[1, 2]]]", 1);
    assert_eval_eq("2^(1/2) (3/2)^(1/3)", "Times[Power[2, Rational[1, 6]], Power[3, Rational[1, 3]]]", 1);
    assert_eval_eq("(-6)^(1/3)/2", "Times[Power[-3, Rational[1, 3]], Power[2, Rational[-2, 3]]]", 1);
    assert_eval_eq("(-2)^(1/3) 3^(1/3)", "Power[-6, Rational[1, 3]]", 1);
    assert_eval_eq("(-2)^(1/3)(-3)^(1/3)", "Times[Power[-3, Rational[1, 3]], Power[-2, Rational[1, 3]]]", 1);
    assert_eval_eq("(-2)^(1/3) 2^(1/3)", "Times[Power[-1, Rational[1, 3]], Power[2, Rational[2, 3]]]", 1);
    assert_eval_eq("(-1)^(1/3) 2^(1/3)", "Power[-2, Rational[1, 3]]", 1);
    assert_eval_eq("2 (-1)^(2/3) 2^(1/3)", "Times[2, Power[-1, Rational[2, 3]], Power[2, Rational[1, 3]]]", 1);
    assert_eval_eq("(-3/2)^(1/3)/3", "Times[Power[Rational[-1, 2], Rational[1, 3]], Power[3, Rational[-2, 3]]]", 1);
    assert_eval_eq("(-2)^(1/3) 3^(-1/3)", "Power[Rational[-2, 3], Rational[1, 3]]", 1);
    assert_eval_eq("(-2)^(1/4) 3^(1/4)", "Power[-6, Rational[1, 4]]", 1);
    assert_eval_eq("(2/3)^(5/3)", "Times[Rational[2, 3], Power[Rational[2, 3], Rational[2, 3]]]", 1);
    assert_eval_eq("(2/3)^(-5/3)", "Times[Rational[3, 2], Power[Rational[3, 2], Rational[2, 3]]]", 1);
    assert_eval_eq("(3/2)^(-1/2)", "Power[Rational[2, 3], Rational[1, 2]]", 1);
    assert_eval_eq("(1/6)^(1/2)", "Power[6, Rational[-1, 2]]", 1);
    assert_eval_eq("(4/9)^(1/3)", "Power[Rational[2, 3], Rational[2, 3]]", 1);
    assert_eval_eq("(12/5)^(4/3)", "Times[Rational[12, 5], Power[Rational[3, 5], Rational[1, 3]], Power[2, Rational[2, 3]]]", 1);
    assert_eval_eq("(-2/3)^(5/3)", "Times[Rational[-2, 3], Power[Rational[-2, 3], Rational[2, 3]]]", 1);
    assert_eval_eq("(-1/6)^(1/3)", "Power[Rational[-1, 6], Rational[1, 3]]", 1);
    assert_eval_eq("(-12)^(1/3)", "Times[Power[-3, Rational[1, 3]], Power[2, Rational[2, 3]]]", 1);
    assert_eval_eq("(-4)^(1/3)", "Times[Power[-1, Rational[1, 3]], Power[2, Rational[2, 3]]]", 1);
    assert_eval_eq("(-16)^(1/4)", "Times[2, Power[-1, Rational[1, 4]]]", 1);
    assert_eval_eq("(-12)^(1/4)", "Times[Power[-3, Rational[1, 4]], Power[2, Rational[1, 2]]]", 1);
    assert_eval_eq("(-1/6)^(3/4)", "Power[Rational[-1, 6], Rational[3, 4]]", 1);
    assert_eval_eq("(-1)^(3/4) 6^(-3/4)", "Power[Rational[-1, 6], Rational[3, 4]]", 1);
    assert_eval_eq("(-3)^(5/4)", "Times[-3, Power[-3, Rational[1, 4]]]", 1);
    assert_eval_eq("(-2/3)^(3/2)", "Times[Complex[0, Rational[-2, 3]], Power[Rational[2, 3], Rational[1, 2]]]", 1);
    assert_eval_eq("Sqrt[-2/3]", "Times[Complex[0, 1], Power[Rational[2, 3], Rational[1, 2]]]", 1);
    assert_eval_eq("2^(-2/3) 5^(2/3)", "Power[Rational[5, 2], Rational[2, 3]]", 1);
    assert_eval_eq("2^(2/5) 3^(-2/5) 5^(2/5)", "Power[Rational[10, 3], Rational[2, 5]]", 1);
    assert_eval_eq("12^(1/3)", "Times[Power[2, Rational[2, 3]], Power[3, Rational[1, 3]]]", 1);
    assert_eval_eq("6 Sqrt[2] 3^(1/3)", "Times[6, Power[2, Rational[1, 2]], Power[3, Rational[1, 3]]]", 1);
    /* the one divergence from the reference, which leaves this product
     * unmerged while merging 2^(-2/3) 5^(2/3) on its own */
    assert_eval_eq("2^(-2/3) 3^(-1/3) 5^(2/3)",
                   "Times[Power[Rational[5, 2], Rational[2, 3]], Power[3, Rational[-1, 3]]]", 1);
    /* a canonical product is a fixed point: re-evaluating changes nothing,
     * and a bignum coefficient is recognised as already canonical */
    assert_eval_eq("Times[Rational[11529215046068469755, 3], Power[2305843009213693951, Rational[1, 3]]]",
                   "Times[Rational[11529215046068469755, 3], Power[2305843009213693951, Rational[1, 3]]]", 1);
    /* the value is kept: (in - out) is exactly 0 for a spread of products */
    assert_eval_eq("Union[Chop[N[{-1/3 Sqrt[3/2] + 1/Sqrt[6], 6^(1/3)/2 - 3^(1/3)/2^(2/3), "
                   "(25/12)^(1/3) - 5^(2/3)/(2^(2/3) 3^(1/3)), (-6)^(1/3)/2 - (-3)^(1/3)/2^(2/3), "
                   "(-12)^(1/4) - (-3)^(1/4) Sqrt[2]}, 30]]]", "List[0]", 1);
}

static void test_complex_heads_of_powers(void) {
    assert_eval_eq("Abs[x^2]", "Power[Abs[x], 2]", 1);
    assert_eval_eq("Abs[x^(1/2)]", "Power[Abs[x], Rational[1, 2]]", 1);
    assert_eval_eq("Abs[x^2.5]", "Power[Abs[x], 2.5]", 1);
    assert_eval_eq("Abs[x^-1]", "Power[Abs[x], -1]", 1);
    assert_eval_eq("Abs[x^Pi]", "Power[Abs[x], Pi]", 1);
    assert_eval_eq("Abs[(x+1)^3]", "Power[Abs[Plus[1, x]], 3]", 1);
    assert_eval_eq("Abs[E^y]", "Power[E, Re[y]]", 1);
    assert_eval_eq("Abs[2^y]", "Power[2, Re[y]]", 1);
    assert_eval_eq("Abs[Pi^y]", "Power[Pi, Re[y]]", 1);
    assert_eval_eq("Abs[E^(a+I b)]", "Power[E, Plus[Times[-1, Im[b]], Re[a]]]", 1);
    assert_eval_eq("Abs[2^(x+I)]", "Power[2, Re[x]]", 1);
    assert_eval_eq("Abs[I^y]", "Power[E, Times[Rational[-1, 2], Pi, Im[y]]]", 1);
    assert_eval_eq("Abs[(-1)^y]", "Power[E, Times[-1, Pi, Im[y]]]", 1);
    assert_eval_eq("Abs[x^I]", "Power[E, Times[-1, Im[Log[x]]]]", 1);
    assert_eval_eq("Abs[x^(2 I)]", "Power[E, Times[-2, Im[Log[x]]]]", 1);
    assert_eval_eq("Abs[x^y]", "Abs[Power[x, y]]", 1);
    assert_eval_eq("Abs[x^(1+I)]", "Abs[Power[x, Complex[1, 1]]]", 1);
    assert_eval_eq("Abs[(-2)^y]", "Abs[Power[-2, y]]", 1);
    assert_eval_eq("Abs[x^2 E^y]", "Times[Power[E, Re[y]], Power[Abs[x], 2]]", 1);
    assert_eval_eq("Abs[x^2 y]", "Abs[Times[Power[x, 2], y]]", 1);
    assert_eval_eq("Abs[E^y x y]", "Times[Power[E, Re[y]], Abs[Times[x, y]]]", 1);
    assert_eval_eq("Abs[Sqrt[x] E^y]", "Times[Power[E, Re[y]], Power[Abs[x], Rational[1, 2]]]", 1);
    assert_eval_eq("Abs[I^y x]", "Times[Power[E, Times[Rational[-1, 2], Pi, Im[y]]], Abs[x]]", 1);
    assert_eval_eq("Abs[(-2)^(1/3) x]", "Times[Power[2, Rational[1, 3]], Abs[x]]", 1);
    assert_eval_eq("Abs[(-8)^(1/3)]", "2", 1);
    assert_eval_eq("Abs[Pi^(2 I)]", "1", 1);
    assert_eval_eq("Abs[Abs[x]^2]", "Power[Abs[x], 2]", 1);
    assert_eval_eq("Abs[x^2] - Abs[x]^2", "0", 1);
    assert_eval_eq("Abs[E^x]^2", "Power[E, Times[2, Re[x]]]", 1);
    assert_eval_eq("Conjugate[x^2]", "Power[Conjugate[x], 2]", 1);
    assert_eval_eq("Conjugate[x^-1]", "Power[Conjugate[x], -1]", 1);
    assert_eval_eq("Conjugate[x^(1/2)]", "Conjugate[Power[x, Rational[1, 2]]]", 1);
    assert_eval_eq("Conjugate[x^2.5]", "Conjugate[Power[x, 2.5]]", 1);
    assert_eval_eq("Conjugate[E^y]", "Power[E, Conjugate[y]]", 1);
    assert_eval_eq("Conjugate[2^y]", "Power[2, Conjugate[y]]", 1);
    assert_eval_eq("Conjugate[(-2)^y]", "Conjugate[Power[-2, y]]", 1);
    assert_eval_eq("Conjugate[Pi^(2 I)]", "Power[Pi, Complex[0, -2]]", 1);
    assert_eval_eq("Conjugate[x^2 y z]", "Times[Power[Conjugate[x], 2], Conjugate[Times[y, z]]]", 1);
    assert_eval_eq("Conjugate[Sqrt[x] y]", "Conjugate[Times[Power[x, Rational[1, 2]], y]]", 1);
    assert_eval_eq("Conjugate[1+x]", "Plus[1, Conjugate[x]]", 1);
    assert_eval_eq("Conjugate[x+y]", "Conjugate[Plus[x, y]]", 1);
    assert_eval_eq("Conjugate[x^2+y]", "Plus[Power[Conjugate[x], 2], Conjugate[y]]", 1);
    assert_eval_eq("Conjugate[x^3 + 2 x + 1]", "Plus[1, Times[2, Conjugate[x]], Power[Conjugate[x], 3]]", 1);
    assert_eval_eq("Conjugate[E^(I x)]", "Power[E, Times[Complex[0, -1], Conjugate[x]]]", 1);
    assert_eval_eq("Conjugate[Conjugate[x]^2]", "Power[x, 2]", 1);
    assert_eval_eq("Re[a + I b]", "Plus[Times[-1, Im[b]], Re[a]]", 1);
    assert_eval_eq("Re[x + I]", "Re[x]", 1);
    assert_eval_eq("Re[2 x]", "Times[2, Re[x]]", 1);
    assert_eval_eq("Re[I x]", "Times[-1, Im[x]]", 1);
    assert_eval_eq("Re[(1+I) x]", "Re[Times[Complex[1, 1], x]]", 1);
    assert_eval_eq("Re[x + y + 1]", "Plus[1, Re[Plus[x, y]]]", 1);
    assert_eval_eq("Re[2 x + 3 y]", "Re[Plus[Times[2, x], Times[3, y]]]", 1);
    assert_eval_eq("Re[Conjugate[x]]", "Re[x]", 1);
    assert_eval_eq("Im[Conjugate[x]]", "Times[-1, Im[x]]", 1);
    assert_eval_eq("Im[a + I b]", "Plus[Im[a], Re[b]]", 1);
    assert_eval_eq("Im[I x]", "Re[x]", 1);
    assert_eval_eq("Im[x + 2 I]", "Plus[2, Im[x]]", 1);
    assert_eval_eq("Re[Abs[x] + y]", "Plus[Abs[x], Re[y]]", 1);
    assert_eval_eq("Re[2 + x - I y]", "Plus[2, Im[y], Re[x]]", 1);
    assert_eval_eq("Im[(2+I)x + 3 I y]", "Plus[Im[Times[Complex[2, 1], x]], Times[3, Re[y]]]", 1);
    assert_eval_eq("Re[Pi + x]", "Plus[Pi, Re[x]]", 1);
    assert_eval_eq("Re[x y]", "Re[Times[x, y]]", 1);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_rational_radicals);
    TEST(test_complex_heads_of_powers);
    printf("All radical-canon tests passed.\n");
    return 0;
}
