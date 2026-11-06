/* Tests for CountRoots[poly, x] and CountRoots[poly, {x, a, b}].
 *
 * Expected counts are derived by hand from the factored or classical form of
 * each polynomial (the roots are known), or from theorems: Chebyshev T_n has n
 * simple zeros cos((2k - 1) pi/(2n)) in (-1, 1); x^n - 1 has the n-th roots of
 * unity; Wilkinson's perturbed product prod (x - k) - 2^-23 x^19 has exactly
 * 10 real roots (Wilkinson 1959; an exact Sturm count); the Mignotte
 * polynomial x^n - 2 (10 x - 1)^2 has two roots straddling 1/10 at distance
 * ~ sqrt(2) 10^-(n/2 + 1).
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "test_utils.h"

#include <stdio.h>
#include <string.h>

void test_countroots_real_line() {
    assert_eval_eq("CountRoots[x^2 - 2, x]", "2", 0);
    assert_eval_eq("CountRoots[x^2 + 1, x]", "0", 0);
    assert_eval_eq("CountRoots[x^5 - x - 1, x]", "1", 0);
    /* multiplicities count */
    assert_eval_eq("CountRoots[(x - 1)^3 (x + 2), x]", "4", 0);
    assert_eval_eq("CountRoots[(x^2 + 1)^2 (x - 3)^2, x]", "2", 0);
    assert_eval_eq("CountRoots[x^12 - 1, x]", "2", 0);
    assert_eval_eq("CountRoots[ChebyshevT[20, x], x]", "20", 0);
    /* constants and the zero polynomial */
    assert_eval_eq("CountRoots[5, x]", "0", 0);
    assert_eval_eq("CountRoots[0, x]", "0", 0);
    /* rational and approximate coefficients */
    assert_eval_eq("CountRoots[x^3/3 - x/7 + 1/1000, x]", "3", 0);
    assert_eval_eq("CountRoots[x^2 - 2., x]", "2", 0);
    /* Gaussian coefficients: the real roots of gcd(Re, Im) */
    assert_eval_eq("CountRoots[(x - I) (x - 2), x]", "1", 0);
    assert_eval_eq("CountRoots[(x - I) (x - 2 I), x]", "0", 0);
    /* an unexpanded product: Wilkinson's perturbation */
    assert_eval_eq("CountRoots[Product[x - k, {k, 1, 20}] - x^19/2^23, x]", "10", 0);
}

void test_countroots_interval() {
    /* x^3 - x = x (x - 1) (x + 1): closed intervals include their ends */
    assert_eval_eq("CountRoots[x^3 - x, {x, 0, 1}]", "2", 0);
    assert_eval_eq("CountRoots[x^3 - x, {x, -1, 1}]", "3", 0);
    assert_eval_eq("CountRoots[x^3 - x, {x, 1/2, 3/4}]", "0", 0);
    assert_eval_eq("CountRoots[x^3 - x, {x, 1, -1}]", "3", 0);
    assert_eval_eq("CountRoots[(x - 1)^2 (x - 3), {x, 1, 3}]", "3", 0);
    assert_eval_eq("CountRoots[(x - 1)^2 (x - 3), {x, 2, 4}]", "1", 0);
    /* infinite ends */
    assert_eval_eq("CountRoots[x^2 - 2, {x, 0, Infinity}]", "1", 0);
    assert_eval_eq("CountRoots[x^2 - 2, {x, -Infinity, 0}]", "1", 0);
    assert_eval_eq("CountRoots[(x + 5) (x - 1/3) (x - 7), {x, -Infinity, Infinity}]", "3", 0);
    /* Chebyshev zeros all lie in (-1, 1), none in [1, Infinity) */
    assert_eval_eq("CountRoots[ChebyshevT[15, x], {x, -1, 1}]", "15", 0);
    assert_eval_eq("CountRoots[ChebyshevT[15, x], {x, 1, Infinity}]", "0", 0);
    /* ... and 7 of the 15 are positive (cos((2k - 1) pi/30) > 0 for k <= 7) */
    assert_eval_eq("CountRoots[ChebyshevT[15, x], {x, 1/10^6, 1}]", "7", 0);
    /* Mignotte: the two roots near 1/10, 1.4e-21 apart, both in (0, 1/5) */
    assert_eval_eq("CountRoots[x^40 - 2 (-1 + 10 x)^2, {x, 0, 1/5}]", "2", 0);
    assert_eval_eq("CountRoots[x^40 - 2 (-1 + 10 x)^2, {x, 0, 1/10}]", "1", 0);
    /* approximate ends: their exact binary values */
    assert_eval_eq("CountRoots[x^2 - 2, {x, 1.4, 1.5}]", "1", 0);
}

void test_countroots_box() {
    /* x^2 + 1: both roots on the edges of [0, 0] x [-1, 2] -- closed */
    assert_eval_eq("CountRoots[x^2 + 1, {x, -I, 2 I}]", "2", 0);
    assert_eval_eq("CountRoots[x^2 + 1, {x, 1 - I, 2 + 2 I}]", "0", 0);
    /* x^4 - 1: all four roots on the boundary of the square */
    assert_eval_eq("CountRoots[x^4 - 1, {x, -1 - I, 1 + I}]", "4", 0);
    /* the 12th roots of unity in the closed first-quadrant unit square:
     * angles 0 (corner 1), 30, 60 (inside), 90 (corner I) */
    assert_eval_eq("CountRoots[x^12 - 1, {x, 0, 1 + I}]", "4", 0);
    /* interior roots only: (1 + I)/2, double; (1 - I)/2 and 3 outside */
    assert_eval_eq("CountRoots[(2 x^2 - 2 x + 1)^2 (x - 3), {x, 0, 1 + I}]", "2", 0);
    /* a box corner order other than lower-left / upper-right */
    assert_eval_eq("CountRoots[x^12 - 1, {x, 1, I}]", "4", 0);
    /* every root of a degree-9 polynomial lies in a large enough box */
    assert_eval_eq("CountRoots[x^9 - 3 x^4 + 7, {x, -10 - 10 I, 10 + 10 I}]", "9", 0);
}

void test_countroots_unevaluated() {
    /* a second variable: CountRoots::nuniv */
    assert_eval_eq("Check[CountRoots[x^2 y, x], fail]", "fail", 0);
    assert_eval_eq("Quiet[CountRoots[x^2 y - 1, x]]", "CountRoots[-1 + x^2 y, x]", 0);
    /* not a polynomial: CountRoots::nmet */
    assert_eval_eq("Check[CountRoots[Sin[x], x], fail]", "fail", 0);
    assert_eval_eq("Quiet[CountRoots[Sin[x], x]]", "CountRoots[Sin[x], x]", 0);
    assert_eval_eq("CountRoots[p, x]", "CountRoots[p, x]", 0);
}

void test_countroots_attributes() {
    SymbolDef* d = symtab_get_def("CountRoots");
    ASSERT(d && (d->attributes & ATTR_PROTECTED));
    ASSERT_MSG(d->docstring && strstr(d->docstring, "CountRoots[poly, x]"), "docstring");
}

int main() {
    symtab_init();
    core_init();
    TEST(test_countroots_real_line);
    TEST(test_countroots_interval);
    TEST(test_countroots_box);
    TEST(test_countroots_unevaluated);
    TEST(test_countroots_attributes);
    printf("All CountRoots tests passed.\n");
    return 0;
}
