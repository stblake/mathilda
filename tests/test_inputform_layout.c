/* test_inputform_layout.c
 *
 * InputForm of products, quotients and complex numbers, as the reference
 * writes them (ToString[1/100 a + 5 x, InputForm] is "a/100 + 5*x"):
 * explicit "*", a rational coefficient over the denominator, x^(-1) for a
 * lone reciprocal, -(a*b) for a negated product, -1/2*a for a -1/q
 * coefficient, I/2 and (2*I)/3 for rational imaginary parts.  FeynCalc's
 * Write2 round-trips through this text (7 tests).  Every expected string was
 * produced by the reference kernel.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_inputform_reference_layout(void) {
    assert_eval_eq("ToString[1/100 a + 5 x, InputForm]", "\"a/100 + 5*x\"", 0);
    assert_eval_eq("ToString[a b, InputForm]", "\"a*b\"", 0);
    assert_eval_eq("ToString[2 x, InputForm]", "\"2*x\"", 0);
    assert_eval_eq("ToString[-x, InputForm]", "\"-x\"", 0);
    assert_eval_eq("ToString[-2 x, InputForm]", "\"-2*x\"", 0);
    assert_eval_eq("ToString[x/2, InputForm]", "\"x/2\"", 0);
    assert_eval_eq("ToString[3 x/4, InputForm]", "\"(3*x)/4\"", 0);
    assert_eval_eq("ToString[-3 x/4, InputForm]", "\"(-3*x)/4\"", 0);
    assert_eval_eq("ToString[x^2/3, InputForm]", "\"x^2/3\"", 0);
    assert_eval_eq("ToString[a/b, InputForm]", "\"a/b\"", 0);
    assert_eval_eq("ToString[2 a/b, InputForm]", "\"(2*a)/b\"", 0);
    assert_eval_eq("ToString[a b/(c d), InputForm]", "\"(a*b)/(c*d)\"", 0);
    assert_eval_eq("ToString[1/x, InputForm]", "\"x^(-1)\"", 0);
    assert_eval_eq("ToString[-1/x, InputForm]", "\"-x^(-1)\"", 0);
    assert_eval_eq("ToString[2/x, InputForm]", "\"2/x\"", 0);
    assert_eval_eq("ToString[x^-2, InputForm]", "\"x^(-2)\"", 0);
    assert_eval_eq("ToString[Sqrt[x], InputForm]", "\"Sqrt[x]\"", 0);
    assert_eval_eq("ToString[1/Sqrt[x], InputForm]", "\"1/Sqrt[x]\"", 0);
    assert_eval_eq("ToString[x^(1/3), InputForm]", "\"x^(1/3)\"", 0);
    assert_eval_eq("ToString[1/2, InputForm]", "\"1/2\"", 0);
    assert_eval_eq("ToString[-1/2, InputForm]", "\"-1/2\"", 0);
    assert_eval_eq("ToString[2 + 3 I, InputForm]", "\"2 + 3*I\"", 0);
    assert_eval_eq("ToString[I x, InputForm]", "\"I*x\"", 0);
    assert_eval_eq("ToString[1.5 x, InputForm]", "\"1.5*x\"", 0);
    assert_eval_eq("ToString[(a + b)/2, InputForm]", "\"(a + b)/2\"", 0);
    assert_eval_eq("ToString[(a + b) c, InputForm]", "\"(a + b)*c\"", 0);
    assert_eval_eq("ToString[-(a + b), InputForm]", "\"-a - b\"", 0);
    assert_eval_eq("ToString[2 (a + b), InputForm]", "\"2*(a + b)\"", 0);
    assert_eval_eq("ToString[a - b, InputForm]", "\"a - b\"", 0);
    assert_eval_eq("ToString[f[x]^2, InputForm]", "\"f[x]^2\"", 0);
    assert_eval_eq("ToString[(-1)^n, InputForm]", "\"(-1)^n\"", 0);
    assert_eval_eq("ToString[2^n, InputForm]", "\"2^n\"", 0);
    assert_eval_eq("ToString[x^y^z, InputForm]", "\"x^y^z\"", 0);
    assert_eval_eq("ToString[(x^y)^z, InputForm]", "\"(x^y)^z\"", 0);
    assert_eval_eq("ToString[Pi/4, InputForm]", "\"Pi/4\"", 0);
    assert_eval_eq("ToString[Sqrt[2]/2, InputForm]", "\"1/Sqrt[2]\"", 0);
    assert_eval_eq("ToString[1/(2 Sqrt[3]), InputForm]", "\"1/(2*Sqrt[3])\"", 0);
    assert_eval_eq("ToString[E^x, InputForm]", "\"E^x\"", 0);
    assert_eval_eq("ToString[a.b, InputForm]", "\"a . b\"", 0);
    assert_eval_eq("ToString[-a b, InputForm]", "\"-(a*b)\"", 0);
    assert_eval_eq("ToString[x y/2, InputForm]", "\"(x*y)/2\"", 0);
    assert_eval_eq("ToString[(1/2) Sqrt[x], InputForm]", "\"Sqrt[x]/2\"", 0);
    assert_eval_eq("ToString[4 I Pi, InputForm]", "\"(4*I)*Pi\"", 0);
    assert_eval_eq("ToString[-I, InputForm]", "\"-I\"", 0);
    assert_eval_eq("ToString[Log[x]/x, InputForm]", "\"Log[x]/x\"", 0);
    assert_eval_eq("ToString[(1 + x)^(-1/2), InputForm]", "\"1/Sqrt[1 + x]\"", 0);
    assert_eval_eq("ToString[a^-b, InputForm]", "\"a^(-b)\"", 0);
    assert_eval_eq("ToString[a/(b + c), InputForm]", "\"a/(b + c)\"", 0);
    assert_eval_eq("ToString[{1/2 x}, InputForm]", "\"{x/2}\"", 0);
    assert_eval_eq("ToString[a - 2 b, InputForm]", "\"a - 2*b\"", 0);
    assert_eval_eq("ToString[a - b/2, InputForm]", "\"a - b/2\"", 0);
    assert_eval_eq("ToString[-a/2 - b, InputForm]", "\"-1/2*a - b\"", 0);
    assert_eval_eq("ToString[x^2 y^3, InputForm]", "\"x^2*y^3\"", 0);
    assert_eval_eq("ToString[(x + 1)^2/(x - 1)^3, InputForm]", "\"(1 + x)^2/(-1 + x)^3\"", 0);
    assert_eval_eq("ToString[1/(x y), InputForm]", "\"1/(x*y)\"", 0);
    assert_eval_eq("ToString[-1/(x y), InputForm]", "\"-(1/(x*y))\"", 0);
    assert_eval_eq("ToString[2 x^-1 y, InputForm]", "\"(2*y)/x\"", 0);
    assert_eval_eq("ToString[x/(2 y), InputForm]", "\"x/(2*y)\"", 0);
    assert_eval_eq("ToString[(2 x)^3, InputForm]", "\"8*x^3\"", 0);
    assert_eval_eq("ToString[(x/2)^3, InputForm]", "\"x^3/8\"", 0);
    assert_eval_eq("ToString[(-x)^3, InputForm]", "\"-x^3\"", 0);
    assert_eval_eq("ToString[f[x/2, -y], InputForm]", "\"f[x/2, -y]\"", 0);
    assert_eval_eq("ToString[a -> b/3, InputForm]", "\"a -> b/3\"", 0);
    assert_eval_eq("ToString[{1, 1/2, -x/3}, InputForm]", "\"{1, 1/2, -1/3*x}\"", 0);
    assert_eval_eq("ToString[3.5 x^2, InputForm]", "\"3.5*x^2\"", 0);
    assert_eval_eq("ToString[-2.5 x, InputForm]", "\"-2.5*x\"", 0);
    assert_eval_eq("ToString[x^(-1/2) y, InputForm]", "\"y/Sqrt[x]\"", 0);
    assert_eval_eq("ToString[Sqrt[x] y^-1, InputForm]", "\"Sqrt[x]/y\"", 0);
    assert_eval_eq("ToString[1/(1 + x)^2, InputForm]", "\"(1 + x)^(-2)\"", 0);
    assert_eval_eq("ToString[(1/2) (a + b) c, InputForm]", "\"((a + b)*c)/2\"", 0);
    assert_eval_eq("ToString[-1/2 (a - b), InputForm]", "\"(-a + b)/2\"", 0);
    assert_eval_eq("ToString[I/2, InputForm]", "\"I/2\"", 0);
    assert_eval_eq("ToString[(2 I) x/3, InputForm]", "\"((2*I)/3)*x\"", 0);
    assert_eval_eq("ToString[x^-1 y^-1 z, InputForm]", "\"z/(x*y)\"", 0);
    assert_eval_eq("ToString[2 Pi I, InputForm]", "\"(2*I)*Pi\"", 0);
    assert_eval_eq("ToString[E^(-x/2), InputForm]", "\"E^(-1/2*x)\"", 0);
    assert_eval_eq("ToString[Sin[x]/2 + Cos[x]^2/3, InputForm]", "\"Cos[x]^2/3 + Sin[x]/2\"", 0);
    assert_eval_eq("ToString[a/(b c) + 1, InputForm]", "\"1 + a/(b*c)\"", 0);
    assert_eval_eq("ToString[D[f[x], x]/x, InputForm]", "\"Derivative[1][f][x]/x\"", 0);
    assert_eval_eq("ToString[(-1)^(1/3) x, InputForm]", "\"(-1)^(1/3)*x\"", 0);
    assert_eval_eq("ToString[1/(2 x^2), InputForm]", "\"1/(2*x^2)\"", 0);
    assert_eval_eq("ToString[-3 x^-2, InputForm]", "\"-3/x^2\"", 0);
    assert_eval_eq("ToString[Gamma[1/2]/x, InputForm]", "\"Sqrt[Pi]/x\"", 0);
    assert_eval_eq("ToString[-1/(2 x^2), InputForm]", "\"-1/2*1/x^2\"", 0);
    assert_eval_eq("ToString[-1/Sqrt[x], InputForm]", "\"-(1/Sqrt[x])\"", 0);
    assert_eval_eq("ToString[-1/(2 x), InputForm]", "\"-1/2*1/x\"", 0);
    assert_eval_eq("ToString[-2/(x y), InputForm]", "\"-2/(x*y)\"", 0);
    assert_eval_eq("ToString[-3/(2 x), InputForm]", "\"-3/(2*x)\"", 0);
    assert_eval_eq("ToString[-1/2 a b, InputForm]", "\"-1/2*(a*b)\"", 0);
    assert_eval_eq("ToString[-1/2 x/y, InputForm]", "\"-1/2*x/y\"", 0);
    assert_eval_eq("ToString[-3/2 x/y, InputForm]", "\"(-3*x)/(2*y)\"", 0);
    assert_eval_eq("ToString[I/2, InputForm]", "\"I/2\"", 0);
    assert_eval_eq("ToString[2 I/3, InputForm]", "\"(2*I)/3\"", 0);
    assert_eval_eq("ToString[-I/2, InputForm]", "\"-1/2*I\"", 0);
    assert_eval_eq("ToString[1/2 + I/3, InputForm]", "\"1/2 + I/3\"", 0);
    assert_eval_eq("ToString[3 + I/2, InputForm]", "\"3 + I/2\"", 0);
    assert_eval_eq("ToString[(1 + 2 I) x, InputForm]", "\"(1 + 2*I)*x\"", 0);
    assert_eval_eq("ToString[-2/3 I, InputForm]", "\"(-2*I)/3\"", 0);
    assert_eval_eq("ToString[I x/2, InputForm]", "\"(I/2)*x\"", 0);
    assert_eval_eq("ToString[-I x, InputForm]", "\"(-I)*x\"", 0);
    assert_eval_eq("ToString[-I x/2, InputForm]", "\"(-1/2*I)*x\"", 0);
    assert_eval_eq("ToString[Complex[1/2, -1/2], InputForm]", "\"1/2 - I/2\"", 0);
    assert_eval_eq("ToString[2/3 a b c, InputForm]", "\"(2*a*b*c)/3\"", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_inputform_reference_layout);
    printf("All InputForm layout tests passed.\n");
    return 0;
}
