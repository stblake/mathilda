/*
 * test_solve.c
 *
 * Unit tests for `Solve` (src/solve.c) and the polynomial-equality
 * specialist `Solve`SolvePolynomialEquality` (src/solvepoly.c).
 *
 * Outputs are compared against FullForm strings, so the canonical
 * representation of a Solve result (List[List[Rule[var, val]], ...])
 * is asserted exactly.  Order of solutions matches the per-degree
 * solver -- when the order changes intentionally, update the strings.
 */

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
        printf("FAIL: failed to parse: %s\n", input);
        ASSERT(0);
        return;
    }
    Expr* res = evaluate(e);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        free(res_str);
        expr_free(res);
        expr_free(e);
        ASSERT(0);
        return;
    }
    printf("PASS: %s -> %s\n", input, res_str);
    free(res_str);
    expr_free(res);
    expr_free(e);
}

/* Linear: a x + b == 0  →  -b/a. */
static void test_linear(void) {
    run_test("Solve[2 x + 3 == 0, x]",
             "List[List[Rule[x, Rational[-3, 2]]]]");
    run_test("Solve[2 x + 3 == 1, x]",
             "List[List[Rule[x, -1]]]");
    run_test("Solve[x == 7, x]",
             "List[List[Rule[x, 7]]]");
}

/* Quadratic with real roots (no Reals constraint needed). */
static void test_quadratic_real(void) {
    run_test("Solve[x^2 - 5 x + 6 == 0, x]",
             "List[List[Rule[x, 2]], List[Rule[x, 3]]]");
}

/* Quadratic with imaginary roots; Reals filter prunes to {}. */
static void test_quadratic_complex_and_reals(void) {
    run_test("Solve[x^2 + 1 == 0, x]",
             "List[List[Rule[x, Complex[0, -1]]], "
                  "List[Rule[x, Complex[0, 1]]]]");
    run_test("Solve[x^2 + 1 == 0, x, Reals]",
             "List[]");
}

/* Multiplicity: (x-1)^2 == 0 emits x -> 1 twice. */
static void test_multiplicity(void) {
    run_test("Solve[(x-1)^2 == 0, x]",
             "List[List[Rule[x, 1]], List[Rule[x, 1]]]");
}

/* Pre-factored cubic.  Mathilda's canonical Times order reverses the
 * factors, so roots emerge as 3, 2, 1. */
static void test_pre_factored_cubic(void) {
    run_test("Solve[(x-1)(x-2)(x-3) == 0, x]",
             "List[List[Rule[x, 3]], List[Rule[x, 2]], "
                  "List[Rule[x, 1]]]");
}

/* Binomial x^4 - 1 = 0 over Complexes vs Reals. */
static void test_binomial(void) {
    run_test("Solve[x^4 - 1 == 0, x]",
             "List[List[Rule[x, 1]], List[Rule[x, Complex[0, 1]]], "
                  "List[Rule[x, -1]], List[Rule[x, Complex[0, -1]]]]");
    run_test("Solve[x^4 - 1 == 0, x, Reals]",
             "List[List[Rule[x, -1]], List[Rule[x, 1]]]");
}

/* Root-of-unity binomials x^n + 1 = 0.  The k-th root is (-1)^((2k+1)/n);
 * the polynomial solver constructs each root as the principal n-th root
 * times Power[-1, 2k/n] (rather than Exp[2 k Pi I / n], which would
 * expand via Euler's formula into trigonometric sums for non-special
 * angles).  Power then canonicalises each Power[-1, p/q] to the standard
 * `sign * (-1)^(r/q)` form, matching Mathematica's output. */
static void test_root_of_unity_odd(void) {
    /* x^5 + 1 = 0  ->  (-1)^((2k+1)/5) for k = 0..4. */
    run_test("Solve[x^5 + 1 == 0, x]",
             "List[List[Rule[x, Power[-1, Rational[1, 5]]]], "
                  "List[Rule[x, Power[-1, Rational[3, 5]]]], "
                  "List[Rule[x, -1]], "
                  "List[Rule[x, Times[-1, Power[-1, Rational[2, 5]]]]], "
                  "List[Rule[x, Times[-1, Power[-1, Rational[4, 5]]]]]]");
    /* x^3 + 1 = 0 -- the simplest odd-degree case. */
    run_test("Solve[x^3 + 1 == 0, x]",
             "List[List[Rule[x, Power[-1, Rational[1, 3]]]], "
                  "List[Rule[x, -1]], "
                  "List[Rule[x, Times[-1, Power[-1, Rational[2, 3]]]]]]");
}

/* x^5 - 2 = 0: positive constant.  Roots are 2^(1/5) (-1)^(2k/5) and
 * Power's even-q/odd-q canonicalisation produces the Mathematica form
 * (sign times (-1)^(p/5) times 2^(1/5)) for each k. */
static void test_binomial_positive_constant(void) {
    run_test("Solve[x^5 - 2 == 0, x]",
             "List[List[Rule[x, Power[2, Rational[1, 5]]]], "
                  "List[Rule[x, Times[Power[-1, Rational[2, 5]], "
                                     "Power[2, Rational[1, 5]]]]], "
                  "List[Rule[x, Times[Power[-1, Rational[4, 5]], "
                                     "Power[2, Rational[1, 5]]]]], "
                  "List[Rule[x, Times[-1, Power[-1, Rational[1, 5]], "
                                     "Power[2, Rational[1, 5]]]]], "
                  "List[Rule[x, Times[-1, Power[-1, Rational[3, 5]], "
                                     "Power[2, Rational[1, 5]]]]]]");
}

/* Biquadratic x^4 - 5 x^2 + 4 = 0  --  four real radical roots, no
 * Quartics option required. */
static void test_biquadratic(void) {
    run_test("Solve[x^4 - 5 x^2 + 4 == 0, x]",
             "List[List[Rule[x, 1]], List[Rule[x, -1]], "
                  "List[Rule[x, 2]], List[Rule[x, -2]]]");
    /* Same result regardless of Quartics -> False. */
    run_test("Solve[x^4 - 5 x^2 + 4 == 0, x, Quartics -> False]",
             "List[List[Rule[x, 1]], List[Rule[x, -1]], "
                  "List[Rule[x, 2]], List[Rule[x, -2]]]");
}

/* Cubic default: held Root[] objects.  Cardano kicks in only when
 * Cubics -> True (covered separately). */
static void test_cubic_root_default(void) {
    run_test("Solve[x^3 + x + 1 == 0, x]",
             "List[List[Rule[x, "
               "Root[Function[Plus[1, Slot[1], Power[Slot[1], 3]]], 1]]], "
                  "List[Rule[x, "
               "Root[Function[Plus[1, Slot[1], Power[Slot[1], 3]]], 2]]], "
                  "List[Rule[x, "
               "Root[Function[Plus[1, Slot[1], Power[Slot[1], 3]]], 3]]]]");
}

/* Quintic: 5 held Root[] objects. */
static void test_quintic(void) {
    run_test("Solve[x^5 - x - 1 == 0, x]",
             "List[List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
               "Power[Slot[1], 5]]], 1]]], "
                  "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
               "Power[Slot[1], 5]]], 2]]], "
                  "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
               "Power[Slot[1], 5]]], 3]]], "
                  "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
               "Power[Slot[1], 5]]], 4]]], "
                  "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
               "Power[Slot[1], 5]]], 5]]]]");
}

/* Mixed factor list: one linear factor (root x=1) plus an
 * irreducible quintic (5 Root[] objects). */
static void test_mixed(void) {
    run_test("Solve[(x-1)(x^5 - x - 1) == 0, x]",
             "List[List[Rule[x, 1]], "
                  "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
               "Power[Slot[1], 5]]], 1]]], "
                  "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
               "Power[Slot[1], 5]]], 2]]], "
                  "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
               "Power[Slot[1], 5]]], 3]]], "
                  "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
               "Power[Slot[1], 5]]], 4]]], "
                  "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
               "Power[Slot[1], 5]]], 5]]]]");
}

/* Trivial / tautology / free-of-var cases that simplify before
 * reaching the polynomial dispatcher. */
static void test_trivial(void) {
    run_test("Solve[1 == 0, x]",  "List[]");
    run_test("Solve[0 == 0, x]",  "List[List[]]");
    /* Treat unknown constant y as a parameter (constant in x);
     * y is not provably zero, so no solutions in x. */
    run_test("Solve[y == 0, x]",  "List[]");
}

/* Non-polynomial input: Solve should leave the call unevaluated. */
static void test_non_polynomial(void) {
    run_test("Solve[Sin[x] == 0, x]",
             "Solve[Equal[Sin[x], 0], x]");
}

/* Multivariate vars list (system in two variables) -- the initial cut
 * declines to dispatch and the call stays unevaluated. */
static void test_multivariate_declined(void) {
    run_test("Solve[x + 1 == 0, {x, y}]",
             "Solve[Equal[Plus[x, 1], 0], List[x, y]]");
}

/* Cardano via Cubics -> True for the binomial x^3 - 1 = 0 (the
 * binomial fast-path actually handles this without invoking Cardano,
 * so we test x^3 + x + 1 == 0 with Cubics -> True to exercise the
 * full general-cubic branch). */
static void test_cubic_radical(void) {
    /* Just count solutions for the Cardano output: pre-computing the
     * exact FullForm is brittle (algebraic simplification rounds it
     * out across versions).  We check that the result is a length-3
     * outer List of length-1 inner Lists of Rules. */
    Expr* e = parse_expression("Solve[x^3 + x + 1 == 0, x, Cubics -> True]");
    ASSERT(e);
    Expr* res = evaluate(e);
    ASSERT(res);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.arg_count == 3);
    for (size_t i = 0; i < 3; i++) {
        Expr* inner = res->data.function.args[i];
        ASSERT(inner->type == EXPR_FUNCTION);
        ASSERT(inner->data.function.arg_count == 1);
        Expr* rule = inner->data.function.args[0];
        ASSERT(rule->type == EXPR_FUNCTION);
        ASSERT(rule->data.function.head->type == EXPR_SYMBOL);
        ASSERT(strcmp(rule->data.function.head->data.symbol, "Rule") == 0);
        /* RHS should NOT be a held Root[] -- Cubics -> True must emit
         * a radical expression. */
        Expr* rhs = rule->data.function.args[1];
        if (rhs->type == EXPR_FUNCTION
            && rhs->data.function.head->type == EXPR_SYMBOL
            && strcmp(rhs->data.function.head->data.symbol, "Root") == 0) {
            printf("FAIL: Cubics -> True still emitted Root[] for "
                   "x^3 + x + 1 == 0\n");
            ASSERT(0);
        }
    }
    printf("PASS: Solve[x^3 + x + 1 == 0, x, Cubics -> True] -- 3 radical rules\n");
    expr_free(res);
    expr_free(e);
}

/* n-quadratic n=3: x^6 - 9 x^3 + 8 = 0  -- u = x^3 gives u in {1, 8},
 * then cube roots.  Just check we get 6 rules. */
static void test_nquadratic_n3(void) {
    Expr* e = parse_expression("Solve[x^6 - 9 x^3 + 8 == 0, x]");
    ASSERT(e);
    Expr* res = evaluate(e);
    ASSERT(res);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.head->data.symbol, "List") == 0);
    if (res->data.function.arg_count != 6) {
        printf("FAIL: Solve[x^6 - 9 x^3 + 8 == 0, x] -- expected 6 rules, got %zu\n",
               res->data.function.arg_count);
        ASSERT(0);
    } else {
        printf("PASS: Solve[x^6 - 9 x^3 + 8 == 0, x] -- 6 rules\n");
    }
    expr_free(res);
    expr_free(e);
}

/* ------------------------------------------------------------------ *
 *  Rational equalities (Together + cross-multiplication path).        *
 *                                                                    *
 *  Each side of the input is brought to a single fraction via         *
 *  Together; the equation N1/D1 == N2/D2 is then cleared to           *
 *  N1*D2 - N2*D1 == 0 and solved polynomially.  Spurious roots        *
 *  introduced at the poles of the original denominators are dropped   *
 *  by the extraneous-root filter (only candidates that *provably*     *
 *  zero D1*D2 are removed; symbolic / undetermined values are kept).  *
 * ------------------------------------------------------------------ */

/* Parametric: a/x + b == 0  →  x = -a/b. */
static void test_rational_parametric(void) {
    run_test("Solve[a/x + b == 0, x]",
             "List[List[Rule[x, Times[-1, a, Power[b, -1]]]]]");
    run_test("Solve[a/x + b == 0, x, Reals]",
             "List[List[Rule[x, Times[-1, a, Power[b, -1]]]]]");
}

/* Concrete one-sided fraction: 1/(x-1) == 2  →  x = 3/2. */
static void test_rational_one_sided(void) {
    run_test("Solve[1/(x-1) == 2, x]",
             "List[List[Rule[x, Rational[3, 2]]]]");
    /* (x+1)/(x-1) == 3  →  3(x-1) = x+1  →  x = 2 */
    run_test("Solve[(x+1)/(x-1) == 3, x]",
             "List[List[Rule[x, 2]]]");
    /* (x-2)/(x-1) == 0  →  x = 2 (no extraneous root at x=1). */
    run_test("Solve[(x-2)/(x-1) == 0, x]",
             "List[List[Rule[x, 2]]]");
}

/* Two-sided cross-multiplication: 1/x == 1/(x+1)  →  no solution. */
static void test_rational_cross_mult(void) {
    /* x+1 = x  →  1 = 0, no soln. */
    run_test("Solve[1/x == 1/(x+1), x]",
             "List[]");
    /* 1/(x-1) - 1/(x+1) == 0  →  2 = 0, no soln. */
    run_test("Solve[1/(x-1) - 1/(x+1) == 0, x]",
             "List[]");
    /* 1/(x-1) + 1/(x+1) == 0  →  2x = 0  →  x = 0. */
    run_test("Solve[1/(x-1) + 1/(x+1) == 0, x]",
             "List[List[Rule[x, 0]]]");
    /* (x-1)/(x+1) == (x+2)/(x-2)  →  x = 0. */
    run_test("Solve[(x-1)/(x+1) == (x+2)/(x-2), x]",
             "List[List[Rule[x, 0]]]");
    /* Parametric cross-mult: a/(x-1) == b/(x+1)  →  x = -(a+b)/(a-b). */
    Expr* e = parse_expression("Solve[a/(x-1) == b/(x+1), x]");
    ASSERT(e);
    Expr* res = evaluate(e);
    ASSERT(res && res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.arg_count == 1);
    printf("PASS: Solve[a/(x-1) == b/(x+1), x] -- 1 parametric rule\n");
    expr_free(res);
    expr_free(e);
}

/* Extraneous-root filtering: solutions that make a cleared
 * denominator vanish must be dropped. */
static void test_rational_extraneous(void) {
    /* x(x-1) = 2(x-1)  →  x ∈ {1, 2}, but x=1 makes (x-1)=0
     * extraneous.  Keep only x=2. */
    run_test("Solve[x/(x-1) == 2/(x-1), x]",
             "List[List[Rule[x, 2]]]");
    /* 2/(x^2-1) == 1/(x-1).  Cross-mult: 2(x-1) = x^2-1
     * = (x-1)(x+1)  ⇒  (x-1)((x+1) - 2) = (x-1)(x-1) = 0.
     * x=1 is extraneous on both sides; result is {}. */
    run_test("Solve[2/(x^2-1) == 1/(x-1), x]",
             "List[]");
    /* (x^2 - 1)/(x - 1) == 0.  Together cancels to (x+1),
     * yielding x = -1 directly (no x=1 to filter). */
    run_test("Solve[(x^2-1)/(x-1) == 0, x]",
             "List[List[Rule[x, -1]]]");
}

/* Mixed polynomial+rational: x + 1/x == 5/2  →  2x^2 - 5x + 2 = 0. */
static void test_rational_mixed_poly(void) {
    run_test("Solve[x + 1/x == 5/2, x]",
             "List[List[Rule[x, Rational[1, 2]]], "
                  "List[Rule[x, 2]]]");
    /* x/(x+1) + (x+1)/x == 2.  Expands to 1 = 0, no soln. */
    run_test("Solve[x/(x+1) + (x+1)/x == 2, x]",
             "List[]");
}

/* Pure-denominator equations: 1/x^2 == 4  →  x^2 = 1/4. */
static void test_rational_denominator_only(void) {
    run_test("Solve[1/x^2 == 4, x]",
             "List[List[Rule[x, Rational[1, 2]]], "
                  "List[Rule[x, Rational[-1, 2]]]]");
    run_test("Solve[1/x^2 == 4, x, Reals]",
             "List[List[Rule[x, Rational[1, 2]]], "
                  "List[Rule[x, Rational[-1, 2]]]]");
    /* 1/x^2 == -4 has no real solution (Reals filter). */
    run_test("Solve[1/x^2 == -4, x, Reals]",
             "List[]");
}

/* Regression: rational equation that cross-multiplies to a quadratic
 * with zero constant term, e.g. a/x + b == c/x  →  b x^2 + (a - c) x == 0.
 * Must factor out x rather than apply the quadratic formula; otherwise
 * a spurious Sqrt[(a-c)^2] term appears that cannot simplify for
 * symbolic a, c. */
static void test_rational_quadratic_zero_const(void) {
    run_test("Solve[a/x + b == c/x, x]",
             "List[List[Rule[x, Times[-1, Power[b, -1], "
                                        "Plus[a, Times[-1, c]]]]]]");
}

/* Sum-of-fractions reduction: 1/(x-1) - 2/(x+1) == 0  →  -x + 3 = 0. */
static void test_rational_sum_of_fractions(void) {
    run_test("Solve[1/(x-1) - 2/(x+1) == 0, x]",
             "List[List[Rule[x, 3]]]");
    /* 1/x + 1/(x-1) == 0  →  (x-1)+x = 0  →  x = 1/2. */
    run_test("Solve[1/x + 1/(x-1) == 0, x]",
             "List[List[Rule[x, Rational[1, 2]]]]");
    /* 1/x^2 - 1/(x+1)^2 == 0  →  (x+1)^2 - x^2 = 2x+1 = 0  →  x = -1/2. */
    run_test("Solve[1/x^2 - 1/(x+1)^2 == 0, x]",
             "List[List[Rule[x, Rational[-1, 2]]]]");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running solve tests...\n");
    TEST(test_linear);
    TEST(test_quadratic_real);
    TEST(test_quadratic_complex_and_reals);
    TEST(test_multiplicity);
    TEST(test_pre_factored_cubic);
    TEST(test_binomial);
    TEST(test_root_of_unity_odd);
    TEST(test_binomial_positive_constant);
    TEST(test_biquadratic);
    TEST(test_cubic_root_default);
    TEST(test_quintic);
    TEST(test_mixed);
    TEST(test_trivial);
    TEST(test_non_polynomial);
    TEST(test_multivariate_declined);
    TEST(test_cubic_radical);
    TEST(test_nquadratic_n3);
    TEST(test_rational_parametric);
    TEST(test_rational_one_sided);
    TEST(test_rational_cross_mult);
    TEST(test_rational_extraneous);
    TEST(test_rational_mixed_poly);
    TEST(test_rational_denominator_only);
    TEST(test_rational_quadratic_zero_const);
    TEST(test_rational_sum_of_fractions);
    printf("All solve tests passed!\n");
    return 0;
}
