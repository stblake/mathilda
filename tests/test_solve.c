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

/* Pre-factored cubic.  The canonical post-solve sort orders real roots
 * by numerical value, so the result emerges as 1, 2, 3 regardless of
 * the factor order in the input. */
static void test_pre_factored_cubic(void) {
    run_test("Solve[(x-1)(x-2)(x-3) == 0, x]",
             "List[List[Rule[x, 1]], List[Rule[x, 2]], "
                  "List[Rule[x, 3]]]");
}

/* Binomial x^4 - 1 = 0 over Complexes vs Reals.  Real roots (-1, 1)
 * come first sorted numerically; complex roots follow with
 * Complex[0, -1] (= -I) before Complex[0, 1] (= I) by canonical
 * argument order. */
static void test_binomial(void) {
    run_test("Solve[x^4 - 1 == 0, x]",
             "List[List[Rule[x, -1]], List[Rule[x, 1]], "
                  "List[Rule[x, Complex[0, -1]]], "
                  "List[Rule[x, Complex[0, 1]]]]");
    run_test("Solve[x^4 - 1 == 0, x, Reals]",
             "List[List[Rule[x, -1]], List[Rule[x, 1]]]");
}

/* Root-of-unity binomials x^n + 1 = 0.  The k-th root is (-1)^((2k+1)/n);
 * the polynomial solver constructs each root as the principal n-th root
 * times Power[-1, 2k/n] (rather than Exp[2 k Pi I / n], which would
 * expand via Euler's formula into trigonometric sums for non-special
 * angles).  Power then canonicalises each Power[-1, p/q] to the standard
 * `sign * (-1)^(r/q)` form.  The canonical Solve sort then orders the
 * results by (-1)-exponent mod 1, matching Mathematica's output (the
 * real root -1 first, then the primitive 2n-th roots interleaved by
 * canonical exponent 1/n, 2/n, ..., (n-1)/n). */
static void test_root_of_unity_odd(void) {
    /* x^5 + 1 = 0  ->  {-1, (-1)^(1/5), -(-1)^(2/5),
     *                   (-1)^(3/5), -(-1)^(4/5)}. */
    run_test("Solve[x^5 + 1 == 0, x]",
             "List[List[Rule[x, -1]], "
                  "List[Rule[x, Power[-1, Rational[1, 5]]]], "
                  "List[Rule[x, Times[-1, Power[-1, Rational[2, 5]]]]], "
                  "List[Rule[x, Power[-1, Rational[3, 5]]]], "
                  "List[Rule[x, Times[-1, Power[-1, Rational[4, 5]]]]]]");
    /* x^3 + 1 = 0 -- the simplest odd-degree case. */
    run_test("Solve[x^3 + 1 == 0, x]",
             "List[List[Rule[x, -1]], "
                  "List[Rule[x, Power[-1, Rational[1, 3]]]], "
                  "List[Rule[x, Times[-1, Power[-1, Rational[2, 3]]]]]]");
}

/* x^5 - 2 = 0: positive constant.  Roots share the magnitude 2^(1/5);
 * the canonical Solve sort orders them by (-1)-exponent mod 1, i.e.
 * 0, 1/5, 2/5, 3/5, 4/5.  Each non-principal root canonicalises via
 * Power's even-q/odd-q normalisation to `sign * (-1)^(p/5) * 2^(1/5)`. */
static void test_binomial_positive_constant(void) {
    run_test("Solve[x^5 - 2 == 0, x]",
             "List[List[Rule[x, Power[2, Rational[1, 5]]]], "
                  "List[Rule[x, Times[-1, Power[-1, Rational[1, 5]], "
                                     "Power[2, Rational[1, 5]]]]], "
                  "List[Rule[x, Times[Power[-1, Rational[2, 5]], "
                                     "Power[2, Rational[1, 5]]]]], "
                  "List[Rule[x, Times[-1, Power[-1, Rational[3, 5]], "
                                     "Power[2, Rational[1, 5]]]]], "
                  "List[Rule[x, Times[Power[-1, Rational[4, 5]], "
                                     "Power[2, Rational[1, 5]]]]]]");
}

/* Biquadratic x^4 - 5 x^2 + 4 = 0  --  four real radical roots, no
 * Quartics option required.  All concrete reals, sorted numerically. */
static void test_biquadratic(void) {
    run_test("Solve[x^4 - 5 x^2 + 4 == 0, x]",
             "List[List[Rule[x, -2]], List[Rule[x, -1]], "
                  "List[Rule[x, 1]], List[Rule[x, 2]]]");
    /* Same result regardless of Quartics -> False. */
    run_test("Solve[x^4 - 5 x^2 + 4 == 0, x, Quartics -> False]",
             "List[List[Rule[x, -2]], List[Rule[x, -1]], "
                  "List[Rule[x, 1]], List[Rule[x, 2]]]");
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

/* Single equation with multi-variable list -- the linear-system
 * specialist solves for the rightmost variable (here `y` is absent
 * from the equation, so `x` becomes the pivot and `y` stays free).
 * Emits Solve::svars; only the pivot variable's rule appears. */
static void test_multivariate_one_eq_one_pivot(void) {
    run_test("Solve[x + 1 == 0, {x, y}]",
             "List[List[Rule[x, -1]]]");
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
             "List[List[Rule[x, Rational[-1, 2]]], "
                  "List[Rule[x, Rational[1, 2]]]]");
    run_test("Solve[1/x^2 == 4, x, Reals]",
             "List[List[Rule[x, Rational[-1, 2]]], "
                  "List[Rule[x, Rational[1, 2]]]]");
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

/* ------------------------------------------------------------------ *
 *  Reals domain: extended coverage.                                   *
 *                                                                    *
 *  Many cases below are regressions for the per-degree branches.      *
 *  In particular the odd-degree-binomial test for negative bases      *
 *  (x^3 + 1, x^5 + 32, x^7 + 128) exercises the sign-aware real-root  *
 *  selector in solve_binomial -- without it, Mathilda would emit the  *
 *  complex principal root `(-c)^(1/n)` instead of `-c^(1/n)`.         *
 * ------------------------------------------------------------------ */

/* Linear over Reals always has the unique solution -b/a, even when    *
 * that value happens to be rational rather than integer. */
static void test_reals_linear(void) {
    run_test("Solve[2 x + 3 == 0, x, Reals]",
             "List[List[Rule[x, Rational[-3, 2]]]]");
    run_test("Solve[x - 7 == 0, x, Reals]",
             "List[List[Rule[x, 7]]]");
}

/* Quadratic, discriminant-aware: Δ > 0 → 2 real, Δ = 0 → 1 real, Δ < 0
 * → 0 real.  D = 0 case ensures the merged single root path works. */
static void test_reals_quadratic_discriminant(void) {
    /* Δ = 1 - 0 - 0... actually 25 - 24 = 1 > 0 → 2 real roots. */
    run_test("Solve[x^2 - 5 x + 6 == 0, x, Reals]",
             "List[List[Rule[x, 2]], List[Rule[x, 3]]]");
    /* Δ = 0: x^2 + 2 x + 1 = (x+1)^2.  Multiplicity-2 root -1 is
     * emitted twice -- the Reals path preserves multiplicity in step
     * with the default Complexes path. */
    run_test("Solve[x^2 + 2 x + 1 == 0, x, Reals]",
             "List[List[Rule[x, -1]], List[Rule[x, -1]]]");
    /* Δ < 0: no real solutions. */
    run_test("Solve[x^2 + 2 x + 5 == 0, x, Reals]",
             "List[]");
    /* Δ < 0 binomial form. */
    run_test("Solve[x^2 + 2 == 0, x, Reals]",
             "List[]");
    /* Irrational real roots:  x^2 - 2 == 0 → ±Sqrt[2].  Canonical
     * order puts -Sqrt[2] first via concrete-vs-symbolic interleave. */
    run_test("Solve[x^2 - 2 == 0, x, Reals]",
             "List[List[Rule[x, Times[-1, Power[2, Rational[1, 2]]]]], "
                  "List[Rule[x, Power[2, Rational[1, 2]]]]]");
}

/* Even-degree binomial: sign of -b/a decides 0 / 1 / 2 real roots. */
static void test_reals_binomial_even(void) {
    /* x^4 - 1 → ±1 over Reals; bare classifier hits the binomial fast
     * path.  (The existing test_binomial covers this too -- duplicated
     * here for the domain coverage matrix.) */
    run_test("Solve[x^4 - 1 == 0, x, Reals]",
             "List[List[Rule[x, -1]], List[Rule[x, 1]]]");
    /* x^4 + 1 == 0 has no real roots. */
    run_test("Solve[x^4 + 1 == 0, x, Reals]",
             "List[]");
    /* x^2 - 16 ==  0 → ±4. */
    run_test("Solve[x^2 - 16 == 0, x, Reals]",
             "List[List[Rule[x, -4]], List[Rule[x, 4]]]");
}

/* Odd-degree binomial: the unique real root.  The negative-base cases
 * here are the regression for the solve_binomial sign-aware selector. */
static void test_reals_binomial_odd(void) {
    /* x^3 + 1 == 0  →  the real cube root of -1 is -1. */
    run_test("Solve[x^3 + 1 == 0, x, Reals]",
             "List[List[Rule[x, -1]]]");
    /* x^3 - 8 == 0  →  real cube root of 8 is 2 (already worked
     * pre-fix because sign is positive). */
    run_test("Solve[x^3 - 8 == 0, x, Reals]",
             "List[List[Rule[x, 2]]]");
    /* x^5 + 32 == 0  →  real fifth root of -32 is -2. */
    run_test("Solve[x^5 + 32 == 0, x, Reals]",
             "List[List[Rule[x, -2]]]");
    /* x^5 - 32 == 0  →  real fifth root of 32 is 2. */
    run_test("Solve[x^5 - 32 == 0, x, Reals]",
             "List[List[Rule[x, 2]]]");
    /* x^7 + 128 == 0  →  real seventh root of -128 is -2. */
    run_test("Solve[x^7 + 128 == 0, x, Reals]",
             "List[List[Rule[x, -2]]]");
    /* x^3 == 0  →  x = 0 (zero-base edge case). */
    run_test("Solve[x^3 == 0, x, Reals]",
             "List[List[Rule[x, 0]]]");
    /* Rational real root: 8 x^3 - 1 == 0  →  x = 1/2. */
    run_test("Solve[8 x^3 - 1 == 0, x, Reals]",
             "List[List[Rule[x, Rational[1, 2]]]]");
}

/* n-quadratic (u = x^n).  All four sub-cases of inner-quadratic root
 * sign × outer-binomial parity exercised. */
static void test_reals_nquadratic(void) {
    /* u^2 - 5 u + 4 == 0 with u = x^2: u ∈ {1, 4}, both ≥ 0 → x in
     * {-2, -1, 1, 2}. */
    run_test("Solve[x^4 - 5 x^2 + 4 == 0, x, Reals]",
             "List[List[Rule[x, -2]], List[Rule[x, -1]], "
                  "List[Rule[x, 1]], List[Rule[x, 2]]]");
    /* u^2 + 5 u + 4 == 0 with u = x^2: u ∈ {-1, -4}, both < 0 →
     * x has no real solution. */
    run_test("Solve[x^4 + 5 x^2 + 4 == 0, x, Reals]",
             "List[]");
    /* u^2 - 3 u - 4 == 0 with u = x^2: u ∈ {-1, 4}.  u=-1 gives no
     * real x; u=4 gives x ∈ {-2, 2}. */
    run_test("Solve[x^4 - 3 x^2 - 4 == 0, x, Reals]",
             "List[List[Rule[x, -2]], List[Rule[x, 2]]]");
    /* n = 3: u^2 - 9 u + 8 == 0  →  u ∈ {1, 8}, cube roots give
     * x ∈ {1, 2}. */
    run_test("Solve[x^6 - 9 x^3 + 8 == 0, x, Reals]",
             "List[List[Rule[x, 1]], List[Rule[x, 2]]]");
}

/* Slow-path factoring lets a polynomial with all-rational roots
 * decompose into degree-1 factors, even at higher degrees -- so
 * `Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x, Reals]` is fully closed-form. */
static void test_reals_factored_cubic(void) {
    run_test("Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x, Reals]",
             "List[List[Rule[x, 1]], List[Rule[x, 2]], "
                  "List[Rule[x, 3]]]");
    /* Mixed factor list: dropping a x^2 + 1 leaves only the real roots. */
    run_test("Solve[(x-1)(x-2)(x^2+1) == 0, x, Reals]",
             "List[List[Rule[x, 1]], List[Rule[x, 2]]]");
}

/* Parametric / symbolic input over Reals: when neither side is
 * concrete, we keep the symbolic answer (the solver cannot prove the
 * candidate is complex). */
static void test_reals_parametric(void) {
    run_test("Solve[a x + b == 0, x, Reals]",
             "List[List[Rule[x, Times[-1, Power[a, -1], b]]]]");
    run_test("Solve[a/x + b == 0, x, Reals]",
             "List[List[Rule[x, Times[-1, a, Power[b, -1]]]]]");
}

/* Tautology / contradiction / non-polynomial under Reals follow the
 * same routing as the default-domain path. */
static void test_reals_edge_cases(void) {
    run_test("Solve[0 == 0, x, Reals]",  "List[List[]]");
    run_test("Solve[1 == 0, x, Reals]",  "List[]");
    /* Non-polynomial input: still unevaluated, including Reals. */
    run_test("Solve[Sin[x] == 0, x, Reals]",
             "Solve[Equal[Sin[x], 0], x, Reals]");
}

/* ------------------------------------------------------------------ *
 *  Integers domain.                                                   *
 *                                                                    *
 *  The Integers filter is applied as a post-pass on the Reals solver: *
 *  every candidate value is checked to be a *provably* concrete       *
 *  integer (EXPR_INTEGER or EXPR_BIGINT).  Rationals, irrationals     *
 *  (Sqrt[2], etc.), held Root[]'s, and symbolic residues are dropped. *
 * ------------------------------------------------------------------ */

/* Linear: integer iff -b/a collapses to an integer literal. */
static void test_integers_linear(void) {
    run_test("Solve[x - 7 == 0, x, Integers]",
             "List[List[Rule[x, 7]]]");
    /* x = -3/2: not an integer → dropped. */
    run_test("Solve[2 x + 3 == 0, x, Integers]",
             "List[]");
    /* Already integer-valued, but as a Rational the canonical form
     * would still be `7`. */
    run_test("Solve[3 x - 21 == 0, x, Integers]",
             "List[List[Rule[x, 7]]]");
}

/* Quadratic: roots that happen to be integers are kept; rationals,
 * irrationals, and complex roots are dropped. */
static void test_integers_quadratic(void) {
    /* Integer roots: 2 and 3. */
    run_test("Solve[x^2 - 5 x + 6 == 0, x, Integers]",
             "List[List[Rule[x, 2]], List[Rule[x, 3]]]");
    /* x = 1/2 (drop) and x = 3 (keep). */
    run_test("Solve[2 x^2 - 7 x + 3 == 0, x, Integers]",
             "List[List[Rule[x, 3]]]");
    /* Irrational roots ±Sqrt[2] dropped. */
    run_test("Solve[x^2 - 2 == 0, x, Integers]",
             "List[]");
    /* Complex roots dropped (would already be empty under Reals). */
    run_test("Solve[x^2 + 1 == 0, x, Integers]",
             "List[]");
    /* Δ = 0, integer double root: emitted twice (multiplicity
     * preserved through the integer filter). */
    run_test("Solve[x^2 + 2 x + 1 == 0, x, Integers]",
             "List[List[Rule[x, -1]], List[Rule[x, -1]]]");
}

/* Binomial: integer iff -b/a is a perfect n-th power. */
static void test_integers_binomial(void) {
    /* x^2 - 4 = 0  →  ±2, both integers. */
    run_test("Solve[x^2 - 4 == 0, x, Integers]",
             "List[List[Rule[x, -2]], List[Rule[x, 2]]]");
    /* x^2 - 5 = 0  →  ±Sqrt[5], dropped. */
    run_test("Solve[x^2 - 5 == 0, x, Integers]",
             "List[]");
    /* x^4 - 1 = 0 over Reals gives ±1, both integers. */
    run_test("Solve[x^4 - 1 == 0, x, Integers]",
             "List[List[Rule[x, -1]], List[Rule[x, 1]]]");
    /* Odd-degree binomial: integer real root. */
    run_test("Solve[x^3 - 8 == 0, x, Integers]",
             "List[List[Rule[x, 2]]]");
    /* Odd-degree binomial with negative base: relies on the
     * sign-aware real-root selector for the value to surface as
     * `-1` rather than `(-1)^(1/3)` (which would not type-match). */
    run_test("Solve[x^3 + 1 == 0, x, Integers]",
             "List[List[Rule[x, -1]]]");
    run_test("Solve[x^5 + 32 == 0, x, Integers]",
             "List[List[Rule[x, -2]]]");
    run_test("Solve[x^5 - 32 == 0, x, Integers]",
             "List[List[Rule[x, 2]]]");
    /* x^5 ± 1.  Integer roots: -1 and 1. */
    run_test("Solve[x^5 + 1 == 0, x, Integers]",
             "List[List[Rule[x, -1]]]");
    run_test("Solve[x^5 - 1 == 0, x, Integers]",
             "List[List[Rule[x, 1]]]");
    /* 8 x^3 - 1 = 0  →  x = 1/2 (Rational, dropped). */
    run_test("Solve[8 x^3 - 1 == 0, x, Integers]",
             "List[]");
}

/* n-quadratic over Integers: same pipeline as Reals, then filter. */
static void test_integers_nquadratic(void) {
    /* u in {1, 4} gives x in {-2, -1, 1, 2}; all integers. */
    run_test("Solve[x^4 - 5 x^2 + 4 == 0, x, Integers]",
             "List[List[Rule[x, -2]], List[Rule[x, -1]], "
                  "List[Rule[x, 1]], List[Rule[x, 2]]]");
    /* All-complex (Reals → {}). */
    run_test("Solve[x^4 + 5 x^2 + 4 == 0, x, Integers]",
             "List[]");
    /* Mixed: only u = 4 gives real x ∈ {-2, 2}; both integers. */
    run_test("Solve[x^4 - 3 x^2 - 4 == 0, x, Integers]",
             "List[List[Rule[x, -2]], List[Rule[x, 2]]]");
    /* n = 3: u in {1, 8} → x in {1, 2}, both integers. */
    run_test("Solve[x^6 - 9 x^3 + 8 == 0, x, Integers]",
             "List[List[Rule[x, 1]], List[Rule[x, 2]]]");
}

/* Factored cubic with integer roots: slow path produces three linear
 * factors and each linear root is integer-typed.  The default
 * `Cubics -> False` does NOT obstruct the result because the cubic
 * is reducible and never reaches the Root[] branch. */
static void test_integers_factored_cubic(void) {
    run_test("Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x, Integers]",
             "List[List[Rule[x, 1]], List[Rule[x, 2]], "
                  "List[Rule[x, 3]]]");
    /* Pre-factored form: same answer. */
    run_test("Solve[(x-1)(x-2)(x-3) == 0, x, Integers]",
             "List[List[Rule[x, 1]], List[Rule[x, 2]], "
                  "List[Rule[x, 3]]]");
    /* Mixed: irreducible x^2+1 → 0 real roots → filtered out; the
     * two integer linear roots survive. */
    run_test("Solve[(x-1)(x-2)(x^2+1) == 0, x, Integers]",
             "List[List[Rule[x, 1]], List[Rule[x, 2]]]");
}

/* Irreducible cubic of degree 3 with no integer roots: the default
 * Cubics -> False path emits 3 Root[] objects which all fail the
 * integer test → result is {}.  An irreducible cubic that happens to
 * have a rational root (impossible -- it would factor) does not arise
 * here; the rational-root case routes via factoring. */
static void test_integers_irreducible_cubic(void) {
    run_test("Solve[x^3 + x + 1 == 0, x, Integers]",
             "List[]");
}

/* Multiplicity is preserved by the integers filter: each integer copy
 * is kept independently. */
static void test_integers_multiplicity(void) {
    run_test("Solve[(x-1)^2 == 0, x, Integers]",
             "List[List[Rule[x, 1]], List[Rule[x, 1]]]");
    run_test("Solve[(x-1)^2 (x-2) == 0, x, Integers]",
             "List[List[Rule[x, 1]], List[Rule[x, 1]], "
                  "List[Rule[x, 2]]]");
}

/* Parametric over Integers: the symbolic answer `Times[...]` is not
 * EXPR_INTEGER → dropped.  This matches Mathematica's "cannot prove
 * integer-valued" stance. */
static void test_integers_parametric(void) {
    run_test("Solve[a x + b == 0, x, Integers]",
             "List[]");
    run_test("Solve[a/x + b == 0, x, Integers]",
             "List[]");
}

/* Rational equations under Integers: the cross-multiplication
 * pipeline still runs, then both the extraneous-root filter and the
 * integer filter apply. */
static void test_integers_rational(void) {
    /* 1/(x-1) == 2 → x = 3/2.  Not integer. */
    run_test("Solve[1/(x-1) == 2, x, Integers]",
             "List[]");
    /* (x+1)/(x-1) == 3 → x = 2.  Integer. */
    run_test("Solve[(x+1)/(x-1) == 3, x, Integers]",
             "List[List[Rule[x, 2]]]");
    /* 1/x^2 == 4 → x = ±1/2.  Rationals dropped. */
    run_test("Solve[1/x^2 == 4, x, Integers]",
             "List[]");
}

/* Tautology / contradiction / non-polynomial under Integers. */
static void test_integers_edge_cases(void) {
    run_test("Solve[0 == 0, x, Integers]",  "List[List[]]");
    run_test("Solve[1 == 0, x, Integers]",  "List[]");
    /* Non-polynomial: unevaluated. */
    run_test("Solve[Sin[x] == 0, x, Integers]",
             "Solve[Equal[Sin[x], 0], x, Integers]");
}

/* Unsupported domains (Rationals, Algebraics, Booleans, Primes) leave
 * the call unevaluated -- the router declines to dispatch.  This
 * guards against accidentally widening the recognised-domain set. */
static void test_unsupported_domains_unevaluated(void) {
    run_test("Solve[x + 1 == 0, x, Rationals]",
             "Solve[Equal[Plus[x, 1], 0], x, Rationals]");
    run_test("Solve[x + 1 == 0, x, Algebraics]",
             "Solve[Equal[Plus[x, 1], 0], x, Algebraics]");
}

/* ============================================================ */
/* Linear-system specialist (Solve`SolveLinearSystem) tests.    */
/* ============================================================ */

/* Two integer-coefficient equations in two unknowns; unique solution. */
static void test_linsys_two_eq_two_var_integers(void) {
    run_test("Solve[3 x + 2 y == 11 && x + y == 12, {x, y}]",
             "List[List[Rule[x, -13], Rule[y, 25]]]");
}

/* Symbolic-coefficient 2x2 system.  Mathilda's canonical form may
 * differ from Mathematica's print form, but the values are
 * mathematically identical:
 *   x = (1 - c)/a              <=>  Mathematica: -((-1 + c)/a)
 *   y = (-2 a + b - b c)/(a d) <=>  Mathematica: -((2 a - b + b c)/(a d))  */
static void test_linsys_symbolic_coefficients(void) {
    run_test("Solve[a x + c == 1 && b x - d y == 2, {x, y}]",
             "List[List[Rule[x, Times[Power[a, -1], "
                              "Plus[1, Times[-1, c]]]], "
                       "Rule[y, Times[Power[a, -1], "
                              "Plus[Times[-2, a], b, Times[-1, Times[b, c]]], "
                              "Power[d, -1]]]]]");
}

/* Single equation, two-variable list -- underdetermined system.
 * The rightmost variable (y) is the pivot; x is free.
 * Solve::svars warning fires; only the rule for y appears. */
static void test_linsys_underdetermined_emits_svars(void) {
    /* y = 11/2 - 3x/2 = (11 - 3x)/2.  Mathilda's Times canonical:
     *   Rational[1, 2] * Plus[11, -3 x]                                  */
    run_test("Solve[3 x + 2 y == 11, {x, y}]",
             "List[List[Rule[y, Times[Rational[1, 2], "
                              "Plus[11, Times[-3, x]]]]]]");
}

/* Three equations in two unknowns; over-determined inconsistent. */
static void test_linsys_overdetermined_inconsistent(void) {
    run_test("Solve[3 x + 2 y == 11 && x + y == 12 && 3 x + y == 32, {x, y}]",
             "List[]");
}

/* Three equations in two unknowns where the third is implied by the
 * first two -- over-determined but *consistent*.  Result is the same
 * as solving just the first two. */
static void test_linsys_overdetermined_consistent(void) {
    run_test("Solve[3 x + 2 y == 11 && x + y == 12 "
             "      && 4 x + 3 y == 23, {x, y}]",
             "List[List[Rule[x, -13], Rule[y, 25]]]");
}

/* List-form (rather than And-form) input is accepted equally. */
static void test_linsys_list_form_input(void) {
    run_test("Solve[{3 x + 2 y == 11, x + y == 12}, {x, y}]",
             "List[List[Rule[x, -13], Rule[y, 25]]]");
}

/* Three-variable system.  Standard textbook example.
 *    x + y + z = 6
 *    2x - y + z = 3
 *    x + 2y - z = 2     ->   {x, y, z} = {1, 2, 3}                    */
static void test_linsys_three_var(void) {
    run_test("Solve[x + y + z == 6 && 2 x - y + z == 3 "
             "      && x + 2 y - z == 2, {x, y, z}]",
             "List[List[Rule[x, 1], Rule[y, 2], Rule[z, 3]]]");
}

/* Rational-coefficient system. */
static void test_linsys_rational_coefficients(void) {
    /* x/2 + y/3 == 1 && x + y == 5/2
     * From first: 3x + 2y == 6.  Combined with 2x + 2y == 5:
     * x == 1, y == 3/2. */
    run_test("Solve[x/2 + y/3 == 1 && x + y == 5/2, {x, y}]",
             "List[List[Rule[x, 1], Rule[y, Rational[3, 2]]]]");
}

/* Reals domain: pure-real coefficients trivially yield real solutions. */
static void test_linsys_reals_domain(void) {
    run_test("Solve[3 x + 2 y == 11 && x + y == 12, {x, y}, Reals]",
             "List[List[Rule[x, -13], Rule[y, 25]]]");
}

/* Integers domain: unique integer solution survives the filter. */
static void test_linsys_integers_integer_solution(void) {
    run_test("Solve[3 x + 2 y == 11 && x + y == 12, {x, y}, Integers]",
             "List[List[Rule[x, -13], Rule[y, 25]]]");
}

/* Integers domain: rational solution is dropped. */
static void test_linsys_integers_rational_dropped(void) {
    /* 2x + y == 1, x + y == 0  ->  x = 1, y = -1.  Integer, kept. */
    run_test("Solve[2 x + y == 1 && x + y == 0, {x, y}, Integers]",
             "List[List[Rule[x, 1], Rule[y, -1]]]");
    /* 2x + y == 0, 4x + y == 1  ->  x = 1/2, y = -1.  Rational x, dropped. */
    run_test("Solve[2 x + y == 0 && 4 x + y == 1, {x, y}, Integers]",
             "List[]");
}

/* Trivial single Equal expressed via And of one (or via List of one). */
static void test_linsys_single_eq_and_one(void) {
    run_test("Solve[{x + y == 3, x - y == 1}, {x, y}]",
             "List[List[Rule[x, 2], Rule[y, 1]]]");
}

/* Tautological row (0 == 0) inside a conjunction is dropped silently;
 * unique solution is recovered from the remaining equations. */
static void test_linsys_zero_row_consistent(void) {
    run_test("Solve[x + y == 3 && 0 == 0 && x - y == 1, {x, y}]",
             "List[List[Rule[x, 2], Rule[y, 1]]]");
}

/* False conjunct collapses the And to False before the specialist
 * sees the input -- Solve returns {} immediately. */
static void test_linsys_false_conjunct(void) {
    run_test("Solve[x + y == 3 && 1 == 0, {x, y}]",
             "List[]");
}

/* Non-affine system declines to dispatch -- left unevaluated. */
static void test_linsys_nonlinear_declined(void) {
    run_test("Solve[x^2 + y == 1 && x - y == 0, {x, y}]",
             "Solve[And[Equal[Plus[Power[x, 2], y], 1], "
                       "Equal[Plus[x, Times[-1, y]], 0]], List[x, y]]");
    /* Cross-term x*y -- non-affine. */
    run_test("Solve[x y == 1 && x + y == 2, {x, y}]",
             "Solve[And[Equal[Times[x, y], 1], Equal[Plus[x, y], 2]], "
                  "List[x, y]]");
}

/* Direct invocation of the qualified builtin Solve`SolveLinearSystem. */
static void test_linsys_qualified_builtin(void) {
    run_test("Solve`SolveLinearSystem[{x + y == 3, x - y == 1}, {x, y}]",
             "List[List[Rule[x, 2], Rule[y, 1]]]");
    run_test("Solve`SolveLinearSystem[{x + y == 0}, {x, y}]",
             "List[List[Rule[y, Times[-1, x]]]]");
}

/* ------------------------------------------------------------------ *
 *  Approximate-number preprocessing.                                  *
 *                                                                    *
 *  An equation containing inexact (Real) numbers is force-rationalised
 *  through src/common.c, dispatched to the exact specialist, and the
 *  bindings are numericalised on the way out so the caller observes
 *  inexact-in / inexact-out semantics (same contract Integrate uses).
 * ------------------------------------------------------------------ */

/* Linear equation with float coefficients: 1.5 x + 3 == 0 -> x == -2.0.
 * Coefficients are rationalised (3/2, 3), solved exactly to -2, then
 * numericalised. */
static void test_approx_linear(void) {
    run_test("Solve[1.5 x + 3 == 0, x]",
             "List[List[Rule[x, -2.0]]]");
}

/* Quadratic with float constant term: x^2 - 2.25 == 0 -> x == ±1.5. */
static void test_approx_quadratic(void) {
    run_test("Solve[x^2 - 2.25 == 0, x]",
             "List[List[Rule[x, -1.5]], List[Rule[x, 1.5]]]");
}

/* Quadratic with all-float coefficients reducing to rational roots:
 * 0.5 x^2 - x - 1.5 == 0  <=>  x^2 - 2 x - 3 == 0  <=>  x ∈ {-1, 3}. */
static void test_approx_quadratic_rational_roots(void) {
    run_test("Solve[0.5 x^2 - x - 1.5 == 0, x]",
             "List[List[Rule[x, -1.0]], List[Rule[x, 3.0]]]");
}

/* Two-variable linear system with float coefficients.  Exercises the
 * linear-system specialist path through the common preprocessor. */
static void test_approx_linsys(void) {
    run_test("Solve[{1.5 x + y == 4.5, x - y == 0.5}, {x, y}]",
             "List[List[Rule[x, 2.0], Rule[y, 1.5]]]");
}

/* Exact integer input must pass straight through the preprocessor with
 * no numericalisation -- common_scan_inexact reports has_inexact ==
 * false, the post-process is skipped, the result stays exact-integer. */
static void test_approx_pure_exact_untouched(void) {
    run_test("Solve[x^2 - 4 == 0, x]",
             "List[List[Rule[x, -2]], List[Rule[x, 2]]]");
}

#ifdef USE_MPFR
/* Precision propagation: pure-MPFR input at 30 decimal digits flows
 * back out at the same 30-digit (≈ 100-bit) precision.  We assert via
 * `Precision[...]` rather than the FullForm of the binding, because
 * the printer collapses trailing zeros on clean-integer-valued MPFR
 * (e.g. `7.0`) so a string match would hide the precision. */
static void test_approx_mpfr_precision_propagation(void) {
    run_test("Precision[(Solve[N[1, 30] x == 7, x][[1]][[1]])[[2]]]",
             "30.103");
}

/* Mixed Real (53 bits) + MPFR (100 bits): the minimum precision wins,
 * so the output is machine precision -- the Real "infects" the answer
 * with its lower precision, matching standard Mathematica semantics. */
static void test_approx_mixed_real_mpfr_picks_min(void) {
    run_test("Precision[(Solve[(0.5 + N[Pi, 30]) x == 1, x][[1]][[1]])[[2]]]",
             "MachinePrecision");
}
#endif

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
    TEST(test_multivariate_one_eq_one_pivot);
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
    TEST(test_reals_linear);
    TEST(test_reals_quadratic_discriminant);
    TEST(test_reals_binomial_even);
    TEST(test_reals_binomial_odd);
    TEST(test_reals_nquadratic);
    TEST(test_reals_factored_cubic);
    TEST(test_reals_parametric);
    TEST(test_reals_edge_cases);
    TEST(test_integers_linear);
    TEST(test_integers_quadratic);
    TEST(test_integers_binomial);
    TEST(test_integers_nquadratic);
    TEST(test_integers_factored_cubic);
    TEST(test_integers_irreducible_cubic);
    TEST(test_integers_multiplicity);
    TEST(test_integers_parametric);
    TEST(test_integers_rational);
    TEST(test_integers_edge_cases);
    TEST(test_unsupported_domains_unevaluated);

    /* Linear-system specialist tests. */
    TEST(test_linsys_two_eq_two_var_integers);
    TEST(test_linsys_symbolic_coefficients);
    TEST(test_linsys_underdetermined_emits_svars);
    TEST(test_linsys_overdetermined_inconsistent);
    TEST(test_linsys_overdetermined_consistent);
    TEST(test_linsys_list_form_input);
    TEST(test_linsys_three_var);
    TEST(test_linsys_rational_coefficients);
    TEST(test_linsys_reals_domain);
    TEST(test_linsys_integers_integer_solution);
    TEST(test_linsys_integers_rational_dropped);
    TEST(test_linsys_single_eq_and_one);
    TEST(test_linsys_zero_row_consistent);
    TEST(test_linsys_false_conjunct);
    TEST(test_linsys_nonlinear_declined);
    TEST(test_linsys_qualified_builtin);

    /* Approximate-number preprocessing (common.c). */
    TEST(test_approx_linear);
    TEST(test_approx_quadratic);
    TEST(test_approx_quadratic_rational_roots);
    TEST(test_approx_linsys);
    TEST(test_approx_pure_exact_untouched);
#ifdef USE_MPFR
    TEST(test_approx_mpfr_precision_propagation);
    TEST(test_approx_mixed_real_mpfr_picks_min);
#endif

    printf("All solve tests passed!\n");
    return 0;
}
