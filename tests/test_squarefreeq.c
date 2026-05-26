#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) {
        printf("Failed to parse: %s\n", input);
        ASSERT(0);
    }
    Expr* res = evaluate(e);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, res_str);
        ASSERT(0);
    } else {
        printf("PASS: %s -> %s\n", input, res_str);
    }
    free(res_str);
    expr_free(res);
    expr_free(e);
}

/* Integers: trivial cases. */
static void test_squarefreeq_trivial_integers(void) {
    run_test("SquareFreeQ[0]",  "False");
    run_test("SquareFreeQ[1]",  "True");
    run_test("SquareFreeQ[-1]", "True");
    run_test("SquareFreeQ[2]",  "True");
    run_test("SquareFreeQ[-2]", "True");
}

/* Small positive integers: textbook cases. */
static void test_squarefreeq_small_integers(void) {
    run_test("SquareFreeQ[6]",   "True");   /* 6 = 2 * 3 */
    run_test("SquareFreeQ[10]",  "True");   /* spec example */
    run_test("SquareFreeQ[30]",  "True");   /* 2 * 3 * 5 */
    run_test("SquareFreeQ[4]",   "False");  /* 2^2 */
    run_test("SquareFreeQ[8]",   "False");  /* 2^3 */
    run_test("SquareFreeQ[9]",   "False");  /* 3^2 */
    run_test("SquareFreeQ[12]",  "False");  /* 2^2 * 3 */
    run_test("SquareFreeQ[20]",  "False");  /* 2^2 * 5  -- spec example */
    run_test("SquareFreeQ[25]",  "False");  /* 5^2 */
    run_test("SquareFreeQ[100]", "False");  /* (2*5)^2 */
}

/* Negative integers: sign is a unit, ignored. */
static void test_squarefreeq_negative_integers(void) {
    run_test("SquareFreeQ[-10]", "True");
    run_test("SquareFreeQ[-4]",  "False");
    run_test("SquareFreeQ[-20]", "False");
}

/* Large bigint cases (>= 10^18). */
static void test_squarefreeq_large_integers(void) {
    /* 10^70 + 3 is prime (Mathematica-confirmed); trivially square-free. */
    run_test("SquareFreeQ[10^70 + 3]", "True");
    /* (2^61 - 1) is the Mersenne prime M61: square-free. */
    run_test("SquareFreeQ[2^61 - 1]",  "True");
    /* (2^31 - 1)^2 = (M31)^2: not square-free. */
    run_test("SquareFreeQ[(2^31 - 1)^2]", "False");
    /* 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 (primorial(9)): square-free. */
    run_test("SquareFreeQ[223092870]", "True");
    /* 2^2 * 3^2 * 5: not square-free. */
    run_test("SquareFreeQ[180]", "False");
}

/* Rational numbers: True iff both numerator and denominator are square-free. */
static void test_squarefreeq_rationals(void) {
    run_test("SquareFreeQ[2/3]",  "True");   /* spec example */
    run_test("SquareFreeQ[6/5]",  "True");
    run_test("SquareFreeQ[1/4]",  "False");  /* denominator 4 = 2^2 */
    run_test("SquareFreeQ[4/3]",  "False");  /* numerator 4 = 2^2 */
    run_test("SquareFreeQ[-2/3]", "True");
    run_test("SquareFreeQ[-8/9]", "False");
}

/* Reals are never "manifestly" square-free per Mathematica. */
static void test_squarefreeq_reals(void) {
    run_test("SquareFreeQ[1.5]", "False");
    run_test("SquareFreeQ[3.0]", "False");
}

/* Gaussian integers with the Automatic mode (Complex literal triggers Z[i]). */
static void test_squarefreeq_gaussian_auto(void) {
    run_test("SquareFreeQ[3 + 2 I]", "True");   /* N=13 prime ≡ 1 (mod 4): split, e=1 */
    run_test("SquareFreeQ[1 + I]",   "True");   /* N=2: (1+i) is Gaussian prime, e=1 */
    run_test("SquareFreeQ[1 - I]",   "True");   /* associate of (1+i) */
    run_test("SquareFreeQ[2 + I]",   "True");   /* N=5 prime ≡ 1 (mod 4), e=1 */
    run_test("SquareFreeQ[3 I]",     "True");   /* 3 ≡ 3 mod 4 Gaussian prime, mult 1 */
    /* 7 + 7 I = 7 * (1 + I); norm = 7^2 * 2 = 98.  7 is itself a Gaussian
     * prime (7 ≡ 3 mod 4) with mult 1; (1 + I) is the Gaussian prime above
     * 2 with mult 1 — square-free. */
    run_test("SquareFreeQ[7 + 7 I]", "True");
}

/* Fix the borderline case in the previous test. */
static void test_squarefreeq_gaussian_squared_factor(void) {
    /* 9 + 9 I = 9 * (1 + I); 9 = 3^2 and 3 is a Gaussian prime → not sqfree. */
    run_test("SquareFreeQ[9 + 9 I]",   "False");
    /* (1 + I)^2 = 2 I: norm 4, mult of (1+I) is 2 → not sqfree. */
    run_test("SquareFreeQ[2 I]",       "False");
    /* (2 + I)^2 = 3 + 4 I: norm 25 = 5^2; pi=(2+i) repeats → not sqfree. */
    run_test("SquareFreeQ[3 + 4 I]",   "False");
    /* (1 + 2 I) * (2 + I) = 5 I: norm 25 = 5^2 but distinct conjugate primes → sqfree. */
    run_test("SquareFreeQ[5 I]",       "True");
}

/* GaussianIntegers option: explicit True overrides default integer behaviour. */
static void test_squarefreeq_gaussian_option(void) {
    /* Default: 2 in Z is prime, square-free. */
    run_test("SquareFreeQ[2]", "True");
    /* With GaussianIntegers -> True: 2 = -i*(1+i)^2 in Z[i] → not sqfree. */
    run_test("SquareFreeQ[2, GaussianIntegers -> True]", "False");
    /* 5 in Z is prime, 5 in Z[i] = (2+i)(2-i) split (two distinct primes). */
    run_test("SquareFreeQ[5]", "True");
    run_test("SquareFreeQ[5, GaussianIntegers -> True]", "True");
    /* GaussianIntegers -> False forces Z (default behaviour). */
    run_test("SquareFreeQ[3 + 2 I, GaussianIntegers -> False]", "False");
    /* Same with Automatic (default). */
    run_test("SquareFreeQ[3 + 2 I, GaussianIntegers -> Automatic]", "True");
}

/* Univariate polynomials. */
static void test_squarefreeq_univariate(void) {
    run_test("SquareFreeQ[6 + 6 x + x^2]",   "True");   /* spec example */
    run_test("SquareFreeQ[x]",               "True");
    run_test("SquareFreeQ[x + 1]",           "True");
    run_test("SquareFreeQ[x^2 - 1]",         "True");   /* (x-1)(x+1) */
    run_test("SquareFreeQ[x^2 + 2 x + 1]",   "False");  /* (x+1)^2 */
    run_test("SquareFreeQ[x^3 - 3 x + 2]",   "False");  /* (x-1)^2 (x+2) */
    run_test("SquareFreeQ[(x - 1)*(x - 2)*(x - 3)]", "True");
    run_test("SquareFreeQ[(x + 1)^2 (x - 1)]", "False");
}

/* Multivariate polynomials. */
static void test_squarefreeq_multivariate(void) {
    run_test("SquareFreeQ[x^3 - x^2 y]",  "False");  /* spec example: x^2 (x - y) */
    run_test("SquareFreeQ[x*y + x + y + 1]", "True"); /* (x+1)(y+1) */
    run_test("SquareFreeQ[(x + y)^2]", "False");
    run_test("SquareFreeQ[(x + y) (x - y)]", "True"); /* x^2 - y^2 */
    run_test("SquareFreeQ[(x + y) (x + y + 1)]", "True");
}

/* Polynomial with explicit vars argument restricting the indeterminates. */
static void test_squarefreeq_with_vars(void) {
    /* x*y^2: viewed as a poly in x alone (with y as a constant in the
     * coefficient ring), it's linear in x → square-free. */
    run_test("SquareFreeQ[x y^2, x]", "True");
    /* Same expression in y: y^2 * (x as constant) → not square-free in y. */
    run_test("SquareFreeQ[x y^2, y]", "False");
    /* Explicit var list: both x and y considered → not square-free. */
    run_test("SquareFreeQ[x y^2, {x, y}]", "False");
    /* (x+1)^2 in y (no y dependence): vacuously sqfree wrt y. */
    run_test("SquareFreeQ[(x + 1)^2, y]", "True");
    /* (x+1)^2 in x: not sqfree. */
    run_test("SquareFreeQ[(x + 1)^2, x]", "False");
}

/* Symbolic expressions that aren't polynomials over the (auto-collected)
 * vars are reported as not-manifestly-square-free → False. */
static void test_squarefreeq_symbolic(void) {
    run_test("SquareFreeQ[Sqrt[2]]",  "False");
    run_test("SquareFreeQ[Sin[x]]",   "False");
    run_test("SquareFreeQ[Log[x]]",   "False");
    run_test("SquareFreeQ[Exp[2 Pi I / 3]]", "False");
    run_test("SquareFreeQ[Pi]",       "False");
}

/* Strings and arbitrary expressions: stay False (never unevaluated). */
static void test_squarefreeq_garbage(void) {
    run_test("SquareFreeQ[\"hello\"]", "False");
    run_test("SquareFreeQ[{1, 2, 3}]", "False");  /* Not Listable. */
}

/* SquareFreeQ is Protected, not Listable. */
static void test_squarefreeq_attributes(void) {
    /* Attributes should be exactly {Protected}. */
    run_test("Attributes[SquareFreeQ]", "List[Protected]");
}

/* Argument-error diagnostics. SquareFreeQ[] returns unevaluated (NULL +
 * stderr diagnostic). Verify the call shape is preserved. */
static void test_squarefreeq_argb(void) {
    fprintf(stderr, "--- next line(s) are an EXPECTED SquareFreeQ::argb diagnostic ---\n");
    run_test("SquareFreeQ[]", "SquareFreeQ[]");
}

/* Bad positional / option arg returns unevaluated. */
static void test_squarefreeq_nonopt(void) {
    fprintf(stderr, "--- next line(s) are an EXPECTED SquareFreeQ::nonopt diagnostic ---\n");
    run_test("SquareFreeQ[1, 2, 3]", "SquareFreeQ[1, 2, 3]");
    fprintf(stderr, "--- next line(s) are an EXPECTED SquareFreeQ::nonopt diagnostic ---\n");
    run_test("SquareFreeQ[4, GaussianIntegers -> Foo]",
             "SquareFreeQ[4, Rule[GaussianIntegers, Foo]]");
}

/* The Modulus option is accepted only for Modulus -> 0 (the default integer
 * ring).  Any other value -- non-zero integer, non-integer, or symbolic --
 * is not yet implemented and must emit SquareFreeQ::modnotimpl while leaving
 * the call unevaluated. */
static void test_squarefreeq_modulus_option(void) {
    run_test("SquareFreeQ[10, Modulus -> 0]", "True");
    fprintf(stderr, "--- next line(s) are an EXPECTED SquareFreeQ::modnotimpl diagnostic ---\n");
    run_test("SquareFreeQ[10, Modulus -> 7]",
             "SquareFreeQ[10, Rule[Modulus, 7]]");
    fprintf(stderr, "--- next line(s) are an EXPECTED SquareFreeQ::modnotimpl diagnostic ---\n");
    run_test("SquareFreeQ[x^2 + 1, Modulus -> 2]",
             "SquareFreeQ[Plus[1, Power[x, 2]], Rule[Modulus, 2]]");
}

/* Regression: the dispatcher must not double-free vars on an early
 * nonopt return. */
static void test_squarefreeq_vars_then_nonopt(void) {
    fprintf(stderr, "--- next line(s) are an EXPECTED SquareFreeQ::nonopt diagnostic ---\n");
    run_test("SquareFreeQ[x y^2, x, 42]",
             "SquareFreeQ[Times[x, Power[y, 2]], x, 42]");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running SquareFreeQ tests...\n");
    TEST(test_squarefreeq_trivial_integers);
    TEST(test_squarefreeq_small_integers);
    TEST(test_squarefreeq_negative_integers);
    TEST(test_squarefreeq_large_integers);
    TEST(test_squarefreeq_rationals);
    TEST(test_squarefreeq_reals);
    TEST(test_squarefreeq_gaussian_auto);
    TEST(test_squarefreeq_gaussian_squared_factor);
    TEST(test_squarefreeq_gaussian_option);
    TEST(test_squarefreeq_univariate);
    TEST(test_squarefreeq_multivariate);
    TEST(test_squarefreeq_with_vars);
    TEST(test_squarefreeq_symbolic);
    TEST(test_squarefreeq_garbage);
    TEST(test_squarefreeq_attributes);
    TEST(test_squarefreeq_argb);
    TEST(test_squarefreeq_nonopt);
    TEST(test_squarefreeq_modulus_option);
    TEST(test_squarefreeq_vars_then_nonopt);
    printf("All SquareFreeQ tests passed!\n");
    return 0;
}
