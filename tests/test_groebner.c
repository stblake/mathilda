/* Unit tests for GroebnerBasis (src/poly/groebner.c +
 * src/poly/groebnerbasis.c).
 *
 * The expected outputs are taken straight from Wolfram Mathematica
 * v13's documentation page for GroebnerBasis; matching them byte-for-
 * byte in FullForm is the strongest correctness signal we have.
 *
 * Run binary directly: ./groebner_tests
 * (per MEMORY.md note: ctest is not configured in tests/CMakeLists.txt). */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"
#include "gb_flint.h"
#include "gb_f4.h"

/* FullForm string match. */
static void check_eq(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, expected) != 0) {
        fprintf(stdout, "FAIL: %s\n  expected: %s\n  got:      %s\n",
                input, expected, got);
        ASSERT_STR_EQ(got, expected);
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* "expr evaluates to True" assertion. */
static void check_true(const char* input) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, "True") != 0) {
        fprintf(stdout, "FAIL: %s\n  expected: True\n  got:      %s\n",
                input, got);
        ASSERT_STR_EQ(got, "True");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Suppress stderr once for nimpl/argt diagnostics. */
static void mute_stderr_once(void) {
    static int done = 0;
    if (!done) { test_mute_messages(); done = 1; }
}

/* ------------------------------------------------------------------ */
/* 1. Headline spec example: bivariate                                 */
/* ------------------------------------------------------------------ */

static void test_bivariate_two_polys(void) {
    /* GroebnerBasis[{x^2-2y^2, x y-3}, {x, y}] -> {-9 + 2 y^4, 3 x - 2 y^3} */
    check_eq("GroebnerBasis[{x^2-2y^2, x y-3}, {x, y}]",
             "List[Plus[-9, Times[2, Power[y, 4]]], "
                  "Plus[Times[3, x], Times[-2, Power[y, 3]]]]");
}

/* ------------------------------------------------------------------ */
/* 2. No common roots -> {1}                                           */
/* ------------------------------------------------------------------ */

static void test_no_common_roots_three(void) {
    check_eq("GroebnerBasis[{x+y, x^2-1, y^2-2 x}, {x, y}]", "List[1]");
}

static void test_no_common_roots_four(void) {
    /* {x^2+y^2+z^2-1, x y-z+2, z^2-3+x, x-y^2+1} -> {1} */
    check_eq("GroebnerBasis[{x^2+y^2+z^2-1, x y-z+2, z^2-3+x, x-y^2+1}, "
             "{x, y, z}]", "List[1]");
}

/* ------------------------------------------------------------------ */
/* 3. Trivariate, finite roots                                         */
/* ------------------------------------------------------------------ */

static void test_trivariate_finite(void) {
    /* Three-poly basis under default Lex; degree-8 in z. */
    check_eq("GroebnerBasis[{x^2+y^2+z^2-1, x y-z+2, z^2-2 x+3 y}, "
             "{x, y, z}]",
             "List["
             "Plus[1024, Times[-832, z], Times[-215, Power[z, 2]], "
                  "Times[156, Power[z, 3]], Times[-25, Power[z, 4]], "
                  "Times[24, Power[z, 5]], Times[13, Power[z, 6]], "
                  "Power[z, 8]], "
             "Plus[-11552, Times[2560, y], Times[2197, z], "
                  "Times[2764, Power[z, 2]], Times[443, Power[z, 3]], "
                  "Times[728, Power[z, 4]], Times[169, Power[z, 5]], "
                  "Times[32, Power[z, 6]], Times[13, Power[z, 7]]], "
             "Plus[-34656, Times[5120, x], Times[6591, z], "
                  "Times[5732, Power[z, 2]], Times[1329, Power[z, 3]], "
                  "Times[2184, Power[z, 4]], Times[507, Power[z, 5]], "
                  "Times[96, Power[z, 6]], Times[39, Power[z, 7]]]]");
}

/* ------------------------------------------------------------------ */
/* 4. Trivariate, infinite roots (4-poly basis)                        */
/* ------------------------------------------------------------------ */

static void test_trivariate_infinite(void) {
    /* {x^2+y^2+z^2-1, x y-z+2} over {x, y, z} -> 4-poly basis. */
    check_eq("GroebnerBasis[{x^2+y^2+z^2-1, x y-z+2}, {x, y, z}]",
             "List["
             "Plus[4, Times[-1, Power[y, 2]], Power[y, 4], "
                  "Times[-4, z], Power[z, 2], "
                  "Times[Power[y, 2], Power[z, 2]]], "
             "Plus[Times[-2, x], Times[-1, y], Power[y, 3], "
                  "Times[x, z], Times[y, Power[z, 2]]], "
             "Plus[2, Times[x, y], Times[-1, z]], "
             "Plus[-1, Power[x, 2], Power[y, 2], Power[z, 2]]]");
}

/* ------------------------------------------------------------------ */
/* 5. Univariate reduces to PolynomialGCD                              */
/* ------------------------------------------------------------------ */

static void test_univariate_as_gcd_list(void) {
    /* GroebnerBasis[{(x-1)(x-2), (x-2)(x-3)}, {x}] == PolynomialGCD(...) */
    check_eq("GroebnerBasis[{(x-1)(x-2), (x-2)(x-3)}, {x}]",
             "List[Plus[-2, x]]");
}

static void test_univariate_as_gcd_symbol(void) {
    /* Single-symbol shorthand: GroebnerBasis[polys, x] form. */
    check_eq("GroebnerBasis[{(x-1)(x-2), (x-2)(x-3)}, x]",
             "List[Plus[-2, x]]");
}

/* ------------------------------------------------------------------ */
/* 6. Linear -> Gaussian elimination                                   */
/* ------------------------------------------------------------------ */

static void test_linear_gauss(void) {
    check_eq("GroebnerBasis[{x+2y+3z, 4x+5y+6z, 7x+8y+9z}, {x, y, z}]",
             "List[Plus[y, Times[2, z]], Plus[x, Times[-1, z]]]");
}

/* ------------------------------------------------------------------ */
/* 7. Elimination (3-arg form)                                         */
/* ------------------------------------------------------------------ */

static void test_eliminate_xy(void) {
    /* {x^2+y^2+z^2-1, x y z-3} eliminate z -> {9 - x^2y^2 + x^4y^2 + x^2y^4}.
     * Print FullForm. The -1*(x^2 y^2) term is a flat Times[-1, x^2, y^2].
     * gb_to_expr grafts Times[-1, ...] onto an already-built Times[x^2, y^2],
     * and Plus's like-term collection splices the two (v0.355); before that
     * it stayed a non-canonical Times-of-Times. */
    check_eq("GroebnerBasis[{x^2+y^2+z^2-1, x y z-3}, {x, y}, {z}]",
             "List[Plus[9, "
                  "Times[-1, Power[x, 2], Power[y, 2]], "
                  "Times[Power[x, 4], Power[y, 2]], "
                  "Times[Power[x, 2], Power[y, 4]]]]");
}

static void test_eliminate_two_polys(void) {
    /* {x^2+y^2+z^2-1, x y-z+2, z^2-2x+3y} eliminate z -> two polys in {x,y} */
    check_eq("GroebnerBasis[{x^2+y^2+z^2-1, x y-z+2, z^2-2x+3y}, {x, y}, {z}]",
             "List["
             "Plus[28, Times[8, y], Times[41, Power[y, 2]], "
                  "Times[-26, Power[y, 3]], Times[55, Power[y, 4]], "
                  "Times[-8, Power[y, 5]], Times[7, Power[y, 6]], "
                  "Times[-6, Power[y, 7]], Power[y, 8]], "
             "Plus[16, Times[20, x], Times[-10, y], Times[28, Power[y, 2]], "
                  "Times[-57, Power[y, 3]], Times[4, Power[y, 4]], "
                  "Times[-6, Power[y, 5]], Times[6, Power[y, 6]], "
                  "Times[-1, Power[y, 7]]]]");
}

/* ------------------------------------------------------------------ */
/* 8. Equation form                                                    */
/* ------------------------------------------------------------------ */

static void test_equation_form(void) {
    /* GroebnerBasis[{x^2-2y^2==1, x y==3}, {x, y}]
       -> {-9 + y^2 + 2 y^4, 3 x - y - 2 y^3} */
    check_eq("GroebnerBasis[{x^2-2y^2==1, x y==3}, {x, y}]",
             "List["
             "Plus[-9, Power[y, 2], Times[2, Power[y, 4]]], "
             "Plus[Times[3, x], Times[-1, y], Times[-2, Power[y, 3]]]]");
}

/* ------------------------------------------------------------------ */
/* 9. MonomialOrder option                                             */
/* ------------------------------------------------------------------ */

static void test_monomial_order_lex_explicit(void) {
    /* Explicit Lexicographic matches the default. */
    check_eq("GroebnerBasis[{x^2-2y^2, x y-3}, {x, y}, "
             "MonomialOrder -> Lexicographic]",
             "List[Plus[-9, Times[2, Power[y, 4]]], "
                  "Plus[Times[3, x], Times[-2, Power[y, 3]]]]");
}

static void test_monomial_order_degrevlex(void) {
    /* {x^2+y^2+z^2-1, x-z+2, z^2-x y} under DegRevLex
       -> {2+x-z, -2y+y z-z^2, 3+y^2-4z+2z^2, -6+4y+11z-6z^2+3z^3} */
    check_eq("GroebnerBasis[{x^2+y^2+z^2-1, x-z+2, z^2-x y}, {x, y, z}, "
             "MonomialOrder -> DegreeReverseLexicographic]",
             "List["
             "Plus[2, x, Times[-1, z]], "
             "Plus[Times[-2, y], Times[y, z], Times[-1, Power[z, 2]]], "
             "Plus[3, Power[y, 2], Times[-4, z], Times[2, Power[z, 2]]], "
             "Plus[-6, Times[4, y], Times[11, z], Times[-6, Power[z, 2]], "
                  "Times[3, Power[z, 3]]]]");
}

/* Same polys under default Lex -> a totally different basis. */
static void test_monomial_order_lex_differs(void) {
    check_eq("GroebnerBasis[{x^2+y^2+z^2-1, x-z+2, z^2-x y}, {x, y, z}]",
             "List["
             "Plus[12, Times[-28, z], Times[27, Power[z, 2]], "
                  "Times[-12, Power[z, 3]], Times[3, Power[z, 4]]], "
             "Plus[-6, Times[4, y], Times[11, z], "
                  "Times[-6, Power[z, 2]], Times[3, Power[z, 3]]], "
             "Plus[2, x, Times[-1, z]]]");
}

/* ------------------------------------------------------------------ */
/* 10. Edge cases                                                      */
/* ------------------------------------------------------------------ */

static void test_empty_input(void) {
    check_eq("GroebnerBasis[{}, {x, y}]", "List[]");
}

static void test_all_zero_input(void) {
    check_eq("GroebnerBasis[{0, 0}, {x, y}]", "List[]");
}

static void test_constant_in_ideal(void) {
    /* Non-zero constant in the input -> ideal is <1>. */
    check_eq("GroebnerBasis[{5, x + y}, {x, y}]", "List[1]");
}

static void test_single_variable(void) {
    /* Single bare element shorthand on the variable side. */
    check_eq("GroebnerBasis[{x^2-1}, x]",
             "List[Plus[-1, Power[x, 2]]]");
}

/* ------------------------------------------------------------------ */
/* 11. Variable order matters                                          */
/* ------------------------------------------------------------------ */

static void test_var_order_matters(void) {
    /* Same polys, different variable order -> different basis. */
    check_eq("GroebnerBasis[{x^2+y^2-1, y^3-2 x y-3}, {y, x}]",
             "List[Plus[8, Times[4, x], Times[-1, Power[x, 2]], "
                       "Times[-8, Power[x, 3]], Power[x, 4], "
                       "Times[4, Power[x, 5]], Power[x, 6]], "
                  "Plus[-1, Times[2, x], Times[2, Power[x, 2]], "
                       "Times[-2, Power[x, 3]], Times[-1, Power[x, 4]], "
                       "Times[3, y]]]");
}

/* ------------------------------------------------------------------ */
/* 12. Attributes                                                      */
/* ------------------------------------------------------------------ */

static void test_attributes_protected(void) {
    check_true("MemberQ[Attributes[GroebnerBasis], Protected]");
}

/* ------------------------------------------------------------------ */
/* 13. Diagnostic / nimpl fallbacks                                    */
/* ------------------------------------------------------------------ */

static void test_argt_zero(void) {
    mute_stderr_once();
    /* Wrong arity stays unevaluated. */
    check_eq("GroebnerBasis[]", "GroebnerBasis[]");
}

static void test_argt_one(void) {
    mute_stderr_once();
    check_eq("GroebnerBasis[{x^2-1}]",
             "GroebnerBasis[List[Plus[-1, Power[x, 2]]]]");
}

static void test_groebner_walk_matches_buchberger(void) {
    /* Method -> "GroebnerWalk" now computes the basis via the Gröbner walk;
     * by reduced-GB uniqueness the result is identical to Buchberger. */
    check_eq("GroebnerBasis[{x^2-2y^2, x y-3}, {x, y}, "
             "Method -> \"GroebnerWalk\"]",
             "List[Plus[-9, Times[2, Power[y, 4]]], "
                  "Plus[Times[3, x], Times[-2, Power[y, 3]]]]");
}

/* The walk output equals the direct Buchberger output (lex) on a system
 * the engine handles quickly.  Comparing the two bases with == is a
 * self-checking oracle: no hard-coded basis needed. */
static void test_walk_matches_buchberger_3var(void) {
    check_true("GroebnerBasis[{x^2+y^2+z^2-1, x y+y z+z x, x y z-1}, "
               "{x,y,z}] == "
               "GroebnerBasis[{x^2+y^2+z^2-1, x y+y z+z x, x y z-1}, "
               "{x,y,z}, Method -> \"GroebnerWalk\"]");
}

static void test_walk_single_variable(void) {
    /* n_vars == 1: all orders coincide; the walk short-circuits. */
    check_eq("GroebnerBasis[{x^2-1, x^3-x}, {x}, Method -> \"GroebnerWalk\"]",
             "List[Plus[-1, Power[x, 2]]]");
}

static void test_walk_ideal_one(void) {
    /* A unit in the ideal: the basis is {1} in every order. */
    check_eq("GroebnerBasis[{x+1, y+2, 5}, {x, y}, "
             "Method -> \"GroebnerWalk\"]",
             "List[1]");
}

/* The identity weight matrix is exactly Lexicographic. */
static void test_matrix_order_identity_is_lex(void) {
    check_true("GroebnerBasis[{x^2+y^2+z^2-1, x y+y z+z x, x y z-1}, "
               "{x,y,z}, MonomialOrder -> Lexicographic] == "
               "GroebnerBasis[{x^2+y^2+z^2-1, x y+y z+z x, x y z-1}, "
               "{x,y,z}, MonomialOrder -> {{1,0,0},{0,1,0},{0,0,1}}]");
}

/* The grevlex matrix is exactly DegreeReverseLexicographic. */
static void test_matrix_order_grevlex_is_degrevlex(void) {
    check_true("GroebnerBasis[{x^2+y^2+z^2-1, x y+y z+z x, x y z-1}, "
               "{x,y,z}, MonomialOrder -> DegreeReverseLexicographic] == "
               "GroebnerBasis[{x^2+y^2+z^2-1, x y+y z+z x, x y z-1}, "
               "{x,y,z}, MonomialOrder -> {{1,1,1},{0,0,-1},{0,-1,0}}]");
}

/* A genuine weight matrix evaluates without raising nimpl (the issue
 * example).  We only assert it produces a list, not its exact content. */
static void test_matrix_order_weighted_no_error(void) {
    check_eq("Head[GroebnerBasis[{-5x^2+y z-x-1, 2x+3 x y+y^2, "
             "x-3y+x z-2 z^2}, {x,y,z}, "
             "MonomialOrder -> {{1,3,1},{-1,2,0},{4,-3,2}}]]",
             "List");
}

/* The walk reaches the same basis as direct Buchberger for a weight-matrix
 * target order too. */
static void test_matrix_order_walk_matches_direct(void) {
    check_true("GroebnerBasis[{x^2+y^2+z^2-1, x y+y z+z x, x y z-1}, "
               "{x,y,z}, MonomialOrder -> {{1,3,1},{-1,2,0},{4,-3,2}}] == "
               "GroebnerBasis[{x^2+y^2+z^2-1, x y+y z+z x, x y z-1}, "
               "{x,y,z}, MonomialOrder -> {{1,3,1},{-1,2,0},{4,-3,2}}, "
               "Method -> \"GroebnerWalk\"]");
}

/* A rank-deficient matrix is not a term order: nimpl + lex fallback. */
static void test_matrix_order_invalid_falls_back(void) {
    mute_stderr_once();
    check_true("GroebnerBasis[{x^2+y^2+z^2-1, x y+y z+z x, x y z-1}, "
               "{x,y,z}, MonomialOrder -> Lexicographic] == "
               "GroebnerBasis[{x^2+y^2+z^2-1, x y+y z+z x, x y z-1}, "
               "{x,y,z}, MonomialOrder -> {{1,1,1},{2,2,2},{0,0,1}}]");
}

/* CoefficientDomain -> Integers: strong Gröbner basis over Z.  Over Z the
 * leading coefficients are not units, so the strong basis carries extra
 * generators beyond the rational basis (here 4 rather than 2). */
static void test_coeff_domain_integers_strong_basis(void) {
    check_eq("GroebnerBasis[{x^2-2y^2, x y-3}, {x, y}, "
             "CoefficientDomain -> Integers]",
             "List[Plus[-9, Times[2, Power[y, 4]]], "
                  "Plus[Times[3, x], Times[-2, Power[y, 3]]], "
                  "Plus[-3, Times[x, y]], "
                  "Plus[Power[x, 2], Times[-2, Power[y, 2]]]]");
}

/* Strong basis exposes integer GCD combinations: <2x, 3x> = <x> over Z
 * (gcd(2, 3) = 1), and <2x+2, 3x+3> = <x+1>. */
static void test_coeff_domain_integers_gcd_combination(void) {
    check_eq("GroebnerBasis[{2 x, 3 x}, x, CoefficientDomain -> Integers]",
             "List[x]");
    check_eq("GroebnerBasis[{2 x + 2, 3 x + 3}, x, "
             "CoefficientDomain -> Integers]",
             "List[Plus[1, x]]");
}

/* The documented strong-basis example: five generators over Z (the
 * rational basis has three).  `a` is an ordinary variable of Z[x, y, a]. */
static void test_coeff_domain_integers_documented(void) {
    check_eq("GroebnerBasis[{a x^2 + 5 x - 1, 2 x + 3 x y + y^2}, {x, y}, "
             "CoefficientDomain -> Integers]",
             "List[Plus[-4, Times[-12, y], Times[-19, Power[y, 2]], "
                       "Times[-15, Power[y, 3]], Times[a, Power[y, 4]]], "
                  "Plus[18, Times[4, a, x], Times[27, y], "
                       "Times[45, Power[y, 2]], "
                       "Times[2, a, Power[y, 2]], "
                       "Times[-3, a, Power[y, 3]]], "
                  "Plus[Times[2, x], Times[3, x, y], Power[y, 2]], "
                  "Plus[6, Times[2, a, x], Times[9, y], "
                       "Times[a, x, y], Times[15, Power[y, 2]], "
                       "Times[a, Power[y, 2]], "
                       "Times[-1, a, Power[y, 3]]], "
                  "Plus[-1, Times[5, x], Times[a, Power[x, 2]]]]");
}

static void test_deglex_honoured(void) {
    /* DegreeLexicographic is now honoured via the shared named-order weight
     * matrix (it no longer silently falls back to Lexicographic).  The result
     * agrees with the equivalent explicit weight matrix {{1,1},{1,0}} and
     * differs from the Lexicographic basis. */
    check_eq("GroebnerBasis[{x^2-2y^2, x y-3}, {x, y}, "
             "MonomialOrder -> DegreeLexicographic]",
             "List[Plus[-3, Times[x, y]], "
                  "Plus[Power[x, 2], Times[-2, Power[y, 2]]], "
                  "Plus[Times[-3, x], Times[2, Power[y, 3]]]]");
    check_eq("GroebnerBasis[{x^2-2y^2, x y-3}, {x, y}, "
             "MonomialOrder -> DegreeLexicographic] === "
             "GroebnerBasis[{x^2-2y^2, x y-3}, {x, y}, "
             "MonomialOrder -> {{1,1},{1,0}}]",
             "True");
    check_eq("GroebnerBasis[{x^2-2y^2, x y-3}, {x, y}, "
             "MonomialOrder -> DegreeLexicographic] =!= "
             "GroebnerBasis[{x^2-2y^2, x y-3}, {x, y}]",
             "True");
}

/* ------------------------------------------------------------------ */
/* 14. Non-polynomial input -> NULL (unevaluated)                      */
/* ------------------------------------------------------------------ */

static void test_non_polynomial_unevaluated(void) {
    /* Sin[x] is not a polynomial in x: it is an indeterminate of its own
     * (v0.6132), as in the reference: {-1 + Sin[x]^2, x + Sin[x]}. */
    check_eq("GroebnerBasis[{Sin[x] + x, x^2 - 1}, {x}]",
             "List[Plus[-1, Power[Sin[x], 2]], Plus[x, Sin[x]]]");
}

static void test_parametric_coeffs_in_q_of_a(void) {
    /* `a` is auto-discovered as a parameter; the joint variable array
     * is [x, y, a] under lex and the post-Buchberger filter keeps
     * polynomials that mention x or y.  Matches Mathematica's
     * `GroebnerBasis[{a x^2 + 5x - 1, 3 x y + y^2 + 2 x}, {x, y}]`. */
    check_eq("GroebnerBasis[{a x^2+5 x-1, 2 x+3 x y+y^2}, {x, y}]",
             "List[Plus[-4, Times[-12, y], Times[-19, Power[y, 2]], "
                       "Times[-15, Power[y, 3]], "
                       "Times[a, Power[y, 4]]], "
                  "Plus[18, Times[4, a, x], Times[27, y], "
                       "Times[45, Power[y, 2]], "
                       "Times[2, a, Power[y, 2]], "
                       "Times[-3, a, Power[y, 3]]], "
                  "Plus[Times[2, x], Times[3, x, y], "
                       "Power[y, 2]]]");
}

/* CoefficientDomain -> Polynomials[a]: `a` is a coefficient-ring
 * variable, i.e. a parameter.  This is exactly the auto-parameter
 * behaviour, so the basis equals the default. */
static void test_coeff_domain_polynomials_equals_default(void) {
    check_true("GroebnerBasis[{a x^2+5 x-1, 2 x+3 x y+y^2}, {x, y}, "
               "CoefficientDomain -> Polynomials[a]] == "
               "GroebnerBasis[{a x^2+5 x-1, 2 x+3 x y+y^2}, {x, y}]");
}

/* CoefficientDomain -> InexactNumbers[p]: the basis is the exact basis
 * made monic with p-digit real coefficients.  We assert the shape
 * (3 polynomials) and a representative coefficient (the second
 * polynomial's constant term is 18/4 = 4.5) robustly via Round. */
static void test_coeff_domain_inexact_numbers(void) {
    check_true("Length[GroebnerBasis[{a x^2+5 x-1, 2 x+3 x y+y^2}, {x, y}, "
               "CoefficientDomain -> InexactNumbers[20]]] == 3");
    check_true("Round[2 * GroebnerBasis[{a x^2+5 x-1, 2 x+3 x y+y^2}, "
               "{x, y}, CoefficientDomain -> InexactNumbers[20]][[2, 1]]] "
               "== 9");
}

/* Inexact (Real) input coefficients are rationalised before the
 * Buchberger run, so a real-coefficient system still produces a basis
 * (previously gb_from_expr rejected reals and the call was left
 * unevaluated).  1.0 x^2 + ... == x^2 + ... after rationalisation. */
static void test_real_input_rationalised(void) {
    check_true("GroebnerBasis[{1.0 x^2 - 2.0 y^2, x y - 3}, {x, y}] == "
               "GroebnerBasis[{x^2 - 2 y^2, x y - 3}, {x, y}]");
}

/* CoefficientDomain -> RationalFunctions: the parameters live in the
 * coefficient field Q(a), so they are units and a generator whose main
 * leading monomial is a multiple of another's becomes redundant.  Matches
 * Mathematica's documented two-element basis (the default basis with the
 * `2x + 3xy + y^2` generator removed). */
static void test_coeff_domain_rational_functions(void) {
    check_eq("GroebnerBasis[{a x^2+5 x-1, 2 x+3 x y+y^2}, {x, y}, "
             "CoefficientDomain -> RationalFunctions]",
             "List[Plus[-4, Times[-12, y], Times[-19, Power[y, 2]], "
                       "Times[-15, Power[y, 3]], Times[a, Power[y, 4]]], "
                  "Plus[18, Times[4, a, x], Times[27, y], "
                       "Times[45, Power[y, 2]], "
                       "Times[2, a, Power[y, 2]], "
                       "Times[-3, a, Power[y, 3]]]]");
}

/* With no parameters Q(params) == Q, so RationalFunctions must reproduce
 * the ordinary rational basis. */
static void test_coeff_domain_rational_functions_no_params(void) {
    check_true("GroebnerBasis[{x^2 - 2 y^2, x y - 3}, {x, y}, "
               "CoefficientDomain -> RationalFunctions] == "
               "GroebnerBasis[{x^2 - 2 y^2, x y - 3}, {x, y}]");
}

/* Field elimination of x: from {a x + b y, x^2 - y} over Q(a, b) we get
 * b^2 y^2 - a^2 y and a x + b y (substituting x = -b y / a). */
static void test_coeff_domain_rational_functions_eliminates(void) {
    check_eq("GroebnerBasis[{a x + b y, x^2 - y}, {x, y}, "
             "CoefficientDomain -> RationalFunctions]",
             "List[Plus[Times[-1, Power[a, 2], y], "
                       "Times[Power[b, 2], Power[y, 2]]], "
                  "Plus[Times[a, x], Times[b, y]]]");
}

/* ------------------------------------------------------------------ */
/* 15. Ideal-membership sanity                                         */
/* ------------------------------------------------------------------ */

static void test_membership_input_in_ideal(void) {
    /* For the {x^2+y^2-1, y^3-2xy-3} system, every input lies in the
     * ideal generated by the output.  Verify by reducing each input
     * modulo the basis via PolynomialReduce -- but Mathilda doesn't
     * have PolynomialReduce yet, so use the indirect test:
     *   the basis should generate `x^2 + y^2 - 1` and `y^3 - 2 x y - 3`.
     * We check by GroebnerBasis-of-basis being a permutation of the
     * basis (idempotency). */
    check_true(
        "GroebnerBasis[GroebnerBasis[{x^2+y^2-1, y^3-2 x y-3}, {y, x}], {y, x}] "
        " === GroebnerBasis[{x^2+y^2-1, y^3-2 x y-3}, {y, x}]");
}

static void test_idempotency_no_common_roots(void) {
    check_true("GroebnerBasis[GroebnerBasis[{x+y, x^2-1, y^2-2 x}, {x, y}], "
                            "{x, y}] === List[1]");
}

/* ------------------------------------------------------------------ */
/* 16. Memory smoke -- run the biggest spec example several times      */
/* ------------------------------------------------------------------ */

static void test_memory_smoke(void) {
    for (int i = 0; i < 25; i++) {
        Expr* e = parse_expression(
            "GroebnerBasis[{x^2+y^2+z^2-1, x y-z+2, z^2-2 x+3 y}, {x, y, z}]");
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(e);
        expr_free(r);
    }
}

/* ------------------------------------------------------------------ */
/* 17. Two-var Bezout-style example (small ideal)                      */
/* ------------------------------------------------------------------ */

static void test_two_var_simple_basis(void) {
    /* Single binomial generator -> basis is the same single poly. */
    check_eq("GroebnerBasis[{x^2 + y^2}, {x, y}]",
             "List[Plus[Power[x, 2], Power[y, 2]]]");
}

/* ------------------------------------------------------------------ */
/* 18. Issue regressions: implicit parameters, ParameterVariables,     */
/*     subset-of-vars, Sort option, Modulus diagnostic                 */
/* ------------------------------------------------------------------ */

/* Polys mention x but x is omitted from `vars`; x should be auto-
 * promoted to a parameter, and the basis should match the explicit
 * 3-arg-with-parameters form.  Matches Mathematica's
 *   GroebnerBasis[polys, {y, z}]
 * on the example from the issue report. */
static void test_subset_vars_auto_parameter(void) {
    check_eq(
        "GroebnerBasis["
        "  {-5 x^2 + y z - x - 1, 2 x + 3 x y + y^2, x - 3 y + x z - 2 z^2},"
        "  {y, z}]"
        " === "
        "GroebnerBasis["
        "  {-5 x^2 + y z - x - 1, 2 x + 3 x y + y^2, x - 3 y + x z - 2 z^2},"
        "  ParameterVariables -> x]",
        "True");
}

/* `Sort -> True` should be equivalent to reversing the main-variable
 * list before computation.  Hand-built reference. */
static void test_sort_true_reverses_vars(void) {
    check_eq(
        "GroebnerBasis["
        "  {-5 x^2 + y z - x - 1, 2 x + 3 x y + y^2, x - 3 y + x z - 2 z^2},"
        "  {x, y, z}, Sort -> True]"
        " === "
        "GroebnerBasis["
        "  {-5 x^2 + y z - x - 1, 2 x + 3 x y + y^2, x - 3 y + x z - 2 z^2},"
        "  {z, y, x}]",
        "True");
}

/* `Modulus -> p` (p prime) computes a genuine Gröbner basis over GF(p)
 * (gbmod.c).  Small pinned example: the two-circle ideal reduces to a monic
 * lex triangular basis mod 7. */
static void test_modulus_gfp_basis(void) {
    check_eq(
        "GroebnerBasis[{x^2 + y^2 - 1, (x - 1)^2 + y^2 - 1}, {x, y}, Modulus -> 7]"
        " === {1 + y^2, 3 + x}",
        "True");
}

/* The modular basis differs from the rational basis (real GF(p) arithmetic,
 * not the ignored-option stub). */
static void test_modulus_differs_from_rationals(void) {
    check_eq(
        "GroebnerBasis["
        "  {3 x^2 + y z - 5 x - 1, 2 x + 3 x y + y^2, x - 3 y + x z - 2 z^2},"
        "  {x, y, z}, Modulus -> 7]"
        " =!= "
        "GroebnerBasis["
        "  {3 x^2 + y z - 5 x - 1, 2 x + 3 x y + y^2, x - 3 y + x z - 2 z^2},"
        "  {x, y, z}]",
        "True");
}

/* A non-prime modulus is not a field; the option falls back to the rational
 * basis (with a one-shot diagnostic on stderr). */
static void test_modulus_composite_falls_back(void) {
    check_eq(
        "GroebnerBasis[{x^2 + y^2 - 1, x y - 1}, {x, y}, Modulus -> 6]"
        " === "
        "GroebnerBasis[{x^2 + y^2 - 1, x y - 1}, {x, y}]",
        "True");
}

/* `GroebnerBasis[polys, ParameterVariables -> x]` with no positional
 * `vars` argument should auto-derive the main-variable list (= all
 * free symbols in polys minus the parameters). */
static void test_parameter_variables_auto_main(void) {
    check_eq(
        "GroebnerBasis["
        "  {-5 x^2 + y z - x - 1, 2 x + 3 x y + y^2, x - 3 y + x z - 2 z^2},"
        "  ParameterVariables -> x]"
        " === "
        "GroebnerBasis["
        "  {-5 x^2 + y z - x - 1, 2 x + 3 x y + y^2, x - 3 y + x z - 2 z^2},"
        "  {y, z}]",
        "True");
}

/* TimeConstrained must be able to abort an in-flight Buchberger run: cyclic-7
 * in DegreeReverseLexicographic over Q takes far longer than the budget.  The
 * FLINT engine polls the deadline without unwinding, frees its state, and only
 * then unwinds. */
static void test_timeconstrained_aborts_buchberger(void) {
    check_eq(
        "TimeConstrained["
        "  GroebnerBasis[{a + b + c + d + e + f + g,"
        "    a b + b c + c d + d e + e f + f g + g a,"
        "    a b c + b c d + c d e + d e f + e f g + f g a + g a b,"
        "    a b c d + b c d e + c d e f + d e f g + e f g a + f g a b + g a b c,"
        "    a b c d e + b c d e f + c d e f g + d e f g a + e f g a b + f g a b c + g a b c d,"
        "    a b c d e f + b c d e f g + c d e f g a + d e f g a b + e f g a b c + f g a b c d"
        "      + g a b c d e,"
        "    a b c d e f g - 1}, {a, b, c, d, e, f, g},"
        "    MonomialOrder -> DegreeReverseLexicographic],"
        "  1]",
        "$Aborted");
}

/* The issue-2 system, whose Lex basis Buchberger could not reach (the
 * coefficients exploded), now comes from its DRL basis by FGLM: zero-
 * dimensional, a degree-54 eliminant in z, and the inputs reduce to 0. */
static void test_hard_lex_input_by_fglm(void) {
    check_eq(
        "With[{g = GroebnerBasis["
        "    {x y^4 + y z^4 - 2 x^2 y - 3,"
        "     y^4 + x y^2 z + x^2 - 2 x y + y^2 + z^2,"
        "     -x^3 y^2 + x y z^3 + y^4 + x y^2 z - 2 x y},"
        "    {x, y, z}]},"
        "  {Length[g], Exponent[g[[1]], z], Exponent[g[[2]], y],"
        "   Mod[Coefficient[g[[1]], z, 10], 1000003],"
        "   Union[Last /@ (PolynomialReduce[#, g, {x, y, z}] & /@"
        "     {x y^4 + y z^4 - 2 x^2 y - 3, y^4 + x y^2 z + x^2 - 2 x y + y^2 + z^2,"
        "      -x^3 y^2 + x y z^3 + y^4 + x y^2 z - 2 x y})]}]",
        "List[3, 54, 1, 92728, List[0]]");
}

static void test_grevlex_handles_hard_lex_input(void) {
    /* We only assert that the result is a non-empty List; the exact
     * basis is long and not the point. */
    check_true(
        "Length["
        "  GroebnerBasis["
        "    {x y^4 + y z^4 - 2 x^2 y - 3,"
        "     y^4 + x y^2 z + x^2 - 2 x y + y^2 + z^2,"
        "     -x^3 y^2 + x y z^3 + y^4 + x y^2 z - 2 x y},"
        "    {x, y, z}, MonomialOrder -> DegreeReverseLexicographic]] > 0");
}

/* The FLINT engine (DRL, DegLex, Lex direct and by FGLM, GF(p)) gives the
 * classical engines' basis. */
static void test_flint_engine_matches_classical(void) {
    check_eq("GroebnerBasis[{a + b + c + d, a b + b c + c d + d a, a b c + b c d + c d a + d a b,"
             " a b c d - 1}, {a, b, c, d}]",
             "List[Plus[1, Times[-1, Power[c, 2], Power[d, 2]], Times[-1, Power[d, 4]],"
             " Times[Power[c, 2], Power[d, 6]]], Plus[Times[-1, c], Times[-1, d],"
             " Times[Power[c, 3], Power[d, 2]], Times[Power[c, 2], Power[d, 3]]],"
             " Plus[Times[-1, b], Times[-1, d], Times[b, Power[d, 4]], Power[d, 5]],"
             " Plus[Times[b, c], Times[-1, b, d], Times[c, d], Times[-2, Power[d, 2]],"
             " Times[Power[c, 2], Power[d, 4]]], Plus[Power[b, 2], Times[2, b, d], Power[d, 2]],"
             " Plus[a, b, c, d]]");
    check_eq("GroebnerBasis[{x^2 + y + z - 1, x + y^2 + z - 1, x + y + z^2 - 1}, {x, y, z},"
             " MonomialOrder -> DegreeLexicographic]",
             "List[Plus[-1, x, y, Power[z, 2]], Plus[-1, x, Power[y, 2], z], Plus[-1, Power[x, 2], y, z]]");
    check_eq("GroebnerBasis[{x^5 + y^4 + z^3 - 1, x^3 + y^3 + z^2 - 1}, {x, y, z}, Modulus -> 32003] // Length",
             "8");
}


/* ------------------------------------------------------------------ */
/* F4 over GF(p) (gb_f4.c) against the Buchberger engine              */
/* ------------------------------------------------------------------ */

#ifdef USE_FLINT
static void check_f4_matches_buchberger(const char* input) {
    gb_f4_enable(0);
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* ref = expr_to_string_fullform(res);
    expr_free(e); expr_free(res);
    gb_f4_enable(1);
    e = parse_expression(input);
    res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    expr_free(e); expr_free(res);
    if (strcmp(ref, got) != 0)
        fprintf(stdout, "FAIL: %s\n  buchberger: %.300s\n  f4:         %.300s\n", input, ref, got);
    ASSERT_STR_EQ(got, ref);
    free(ref); free(got);
}
#endif

/* The reduced basis is unique: F4 and Buchberger agree on every order, small
 * (< 2^32: one-word accumulators) and 62-bit primes (two-word accumulators),
 * zero- and positive-dimensional ideals, the unit ideal, inputs that are not
 * interreduced (equal and divisible leading monomials). */
static void test_f4_matches_buchberger(void) {
#ifdef USE_FLINT
    static const char* const sys[] = {
        "{a + b + c + d, a b + b c + c d + d a, a b c + b c d + c d a + d a b, a b c d - 1}, {a, b, c, d}",
        "{a + b + c + d + e, a b + b c + c d + d e + e a, a b c + b c d + c d e + d e a + e a b, "
        "a b c d + b c d e + c d e a + d e a b + e a b c, a b c d e - 1}, {a, b, c, d, e}",
        "{x + 2 y + 2 z - 1, x^2 + 2 y^2 + 2 z^2 - x, 2 x y + 2 y z - y}, {x, y, z}",
        "{x^2 + y^2 + z^2 - 1, x - y z}, {x, y, z}",
        "{x^5 + y^4 + z^3 - 1, x^3 + y^3 + z^2 - 1}, {x, y, z}",
        "{x y - 1, x + y, x^2 - y^2 + 3}, {x, y}",
        "{x^2 + x y, x^2 + x y + y, x y^2 - 1, x^2 y}, {x, y}",
        "{3 x^2 y + 2 x y + y + 9 x^2 + 5 x - 3, 2 x^3 y - x y - y + 6 x^3 - 2 x^2 - 3 x + 3, "
        "x^3 y + x^2 y + 3 x^3 + 2 x^2}, {x, y}",
    };
    static const char* const ord[] = {"Lexicographic", "DegreeLexicographic", "DegreeReverseLexicographic"};
    static const char* const mod[] = {"32003", "2147483647", "4611686018427388039"};
    char buf[4096];
    for (size_t i = 0; i < sizeof sys / sizeof sys[0]; i++)
        for (size_t o = 0; o < 3; o++)
            for (size_t m = 0; m < 3; m++) {
                snprintf(buf, sizeof buf, "GroebnerBasis[%s, MonomialOrder -> %s, Modulus -> %s]",
                         sys[i], ord[o], mod[m]);
                check_f4_matches_buchberger(buf);
            }
#endif
}

/* ------------------------------------------------------------------ */
/* The modular shape-position Lex route (gb_flint.c, gf_lex_shape)    */
/* ------------------------------------------------------------------ */

/* FullForm of GroebnerBasis[input] with the shape route on or off. */
static char* gb_fullform(const char* input, int shape_on) {
    gb_flint_shape_route(shape_on);
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    expr_free(e);
    expr_free(res);
    gb_flint_shape_route(1);
    return got;
}

/* The route's basis equals the exact routes' basis (the reduced Lex basis is
 * unique, so equal FullForm), and the route did (or did not) produce it. */
static void check_shape_matches_exact(const char* input, int expect_route) {
    char* exact = gb_fullform(input, 0);
    long before = gb_flint_shape_route(1);
    char* shape = gb_fullform(input, 1);
    long after = gb_flint_shape_route(1);
    if (strcmp(exact, shape) != 0)
        fprintf(stdout, "FAIL: %s\n  exact: %.300s\n  shape: %.300s\n", input, exact, shape);
    ASSERT_STR_EQ(shape, exact);
    if ((after > before) != expect_route) {
        fprintf(stdout, "FAIL: %s: shape route %s\n", input, expect_route ? "did not run" : "ran");
        ASSERT_STR_EQ(expect_route ? "declined" : "ran", expect_route ? "ran" : "declined");
    }
    free(exact); free(shape);
}

/* A deterministic dense polynomial in v[0..nv) of total degree deg, every
 * monomial up to that degree with a nonzero coefficient in [-cmax, cmax]. */
static unsigned long long shape_lcg = 0x9e3779b97f4a7c15ULL;
static long shape_coeff(long cmax) {
    for (;;) {
        shape_lcg = shape_lcg * 6364136223846793005ULL + 1442695040888963407ULL;
        long c = (long)((shape_lcg >> 33) % (unsigned long long)(2 * cmax + 1)) - cmax;
        if (c) return c;
    }
}
static void shape_poly(char* buf, size_t cap, const char* const* v, int nv, int deg, long cmax) {
    int e[6] = {0, 0, 0, 0, 0, 0};
    size_t len = strlen(buf);
    int first = 1;
    for (;;) {
        int tot = 0;
        for (int i = 0; i < nv; i++) tot += e[i];
        if (tot <= deg) {
            len += (size_t)snprintf(buf + len, cap - len, "%s(%ld)", first ? "" : " + ", shape_coeff(cmax));
            for (int i = 0; i < nv; i++)
                if (e[i]) len += (size_t)snprintf(buf + len, cap - len, "*%s^%d", v[i], e[i]);
            first = 0;
        }
        int k = 0;
        while (k < nv && ++e[k] > deg) e[k++] = 0;
        if (k == nv) break;
    }
}

/* Random zero-dimensional systems (generic, so in shape position), 2-4
 * variables, mixed degrees, 1-digit and 6-digit coefficients: the modular
 * route returns the exact routes' basis. */
static void test_shape_route_random_systems(void) {
    static const char* const V[] = {"x", "y", "z", "w"};
    static const int sys[][6] = {           /* nv, degrees..., 6-digit? */
        {2, 2, 3, 0, 0, 0}, {2, 4, 3, 0, 0, 0}, {2, 3, 3, 0, 0, 1}, {3, 2, 2, 2, 0, 0},
        {3, 2, 2, 3, 0, 0}, {3, 2, 2, 2, 0, 1}, {4, 2, 2, 2, 2, 0}, {4, 2, 2, 2, 2, 1},
    };
    for (size_t t = 0; t < sizeof sys / sizeof sys[0]; t++) {
        int nv = sys[t][0];
        long cmax = sys[t][5] ? 999999 : 9;
        char buf[16384];
        snprintf(buf, sizeof buf, "GroebnerBasis[{");
        for (int i = 0; i < nv; i++) {
            if (i) strcat(buf, ", ");
            shape_poly(buf, sizeof buf, V, nv, sys[t][1 + i], cmax);
        }
        strcat(buf, "}, {");
        for (int i = 0; i < nv; i++) { if (i) strcat(buf, ", "); strcat(buf, V[i]); }
        strcat(buf, "}]");
        check_shape_matches_exact(buf, 1);
    }
}

/* The first prime the route uses, NextPrime[2^62], is unlucky for
 * p x^2 + x - y: modulo it the system has 2 solutions, over Q 4.  The next
 * primes agree with each other and replace the reference. */
static void test_shape_route_unlucky_prime(void) {
    check_shape_matches_exact("GroebnerBasis[{NextPrime[2^62] x^2 + x - y, y^2 - 2}, {y, x}]", 1);
    check_shape_matches_exact(
        "GroebnerBasis[{NextPrime[2^62] x^2 + x - y + z, y^2 - 2 + x z, z^2 - 3 x - 1}, {z, y, x}]", 1);
}

/* Where the route must decline, the exact routes answer: not in shape
 * position (cyclic-4, a grid of points), solutions at infinity (the dimension
 * certificate fails), positive-dimensional. */
static void test_shape_route_declines(void) {
    check_shape_matches_exact("GroebnerBasis[{a + b + c + d, a b + b c + c d + d a, "
                              "a b c + b c d + c d a + d a b, a b c d - 1}, {a, b, c, d}]", 0);
    check_shape_matches_exact("GroebnerBasis[{x^2 - 1, y^2 - 1}, {x, y}]", 0);
    check_shape_matches_exact("GroebnerBasis[{x y - 1, y^2 - 2}, {x, y}]", 0);
    check_shape_matches_exact("GroebnerBasis[{x^2 + y^2 + z^2 - 1, x - y z}, {x, y, z}]", 0);
}
/* ------------------------------------------------------------------ */
/* The certified multi-modular route in degree orders (gf_modular_*)  */
/* ------------------------------------------------------------------ */

/* GroebnerBasis[input] with the modular route off (the fraction-free
 * engine) and on must agree exactly, and the route must (or must not) have
 * produced the answer. */
static void check_modular_matches_exact(const char* input, int expect_route) {
    gb_flint_modular_route(0);
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* exact = expr_to_string_fullform(res);
    expr_free(e); expr_free(res);
    long before = gb_flint_modular_route(2);        /* no direct attempt first */
    e = parse_expression(input);
    res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    expr_free(e); expr_free(res);
    long after = gb_flint_modular_route(1);
    if (strcmp(exact, got) != 0)
        fprintf(stdout, "FAIL: %s\n  exact:   %.300s\n  modular: %.300s\n", input, exact, got);
    ASSERT_STR_EQ(got, exact);
    if ((after > before) != expect_route) {
        fprintf(stdout, "FAIL: %s: modular route %s\n", input, expect_route ? "did not run" : "ran");
        ASSERT_STR_EQ(expect_route ? "declined" : "ran", expect_route ? "ran" : "declined");
    }
    free(exact); free(got);
}

/* Random dense systems (zero-dimensional, no solutions at infinity), in DRL
 * and DegLex, small and 6-digit coefficients; a rational input. */
static void test_modular_route_random_systems(void) {
    static const char* const V[] = {"x", "y", "z", "w"};
    static const int sys[][6] = {
        {2, 2, 3, 0, 0, 0}, {3, 2, 2, 2, 0, 0}, {3, 2, 3, 2, 0, 1}, {4, 2, 2, 2, 2, 0},
    };
    static const char* const ord[] = {"DegreeReverseLexicographic", "DegreeLexicographic"};
    for (size_t t = 0; t < sizeof sys / sizeof sys[0]; t++)
        for (size_t o = 0; o < 2; o++) {
            int nv = sys[t][0];
            long cmax = sys[t][5] ? 999999 : 9;
            char buf[16384];
            snprintf(buf, sizeof buf, "GroebnerBasis[{");
            for (int i = 0; i < nv; i++) {
                if (i) strcat(buf, ", ");
                shape_poly(buf, sizeof buf, V, nv, sys[t][1 + i], cmax);
            }
            strcat(buf, "}, {");
            for (int i = 0; i < nv; i++) { if (i) strcat(buf, ", "); strcat(buf, V[i]); }
            strcat(buf, "}, MonomialOrder -> ");
            strcat(buf, ord[o]);
            strcat(buf, "]");
            check_modular_matches_exact(buf, 1);
        }
    check_modular_matches_exact("GroebnerBasis[{x^2/3 + y/7 - 1, y^2 - x/5 + 2/11}, {x, y}, "
                                "MonomialOrder -> DegreeReverseLexicographic]", 1);
}

/* Solutions at infinity (cyclic-5, noon-4): the Hilbert certificate fails
 * for <F^h> and the leading monomials come from the certified basis of the
 * homogenised ideal (DRL); in DegLex there is no such route: the exact
 * engine answers. */
static void test_modular_route_infinity(void) {
    check_modular_matches_exact("GroebnerBasis[{a + b + c + d + e, a b + b c + c d + d e + e a, "
        "a b c + b c d + c d e + d e a + e a b, a b c d + b c d e + c d e a + d e a b + e a b c, "
        "a b c d e - 1}, {a, b, c, d, e}, MonomialOrder -> DegreeReverseLexicographic]", 1);
    check_modular_matches_exact("GroebnerBasis[{10 x (y^2 + z^2) - 11 x + 10, 10 y (x^2 + z^2) - 11 y + 10, "
        "10 z (x^2 + y^2) - 11 z + 10}, {x, y, z}, MonomialOrder -> DegreeReverseLexicographic]", 1);
    check_modular_matches_exact("GroebnerBasis[{x y - 1, y^2 - 2}, {x, y}, "
        "MonomialOrder -> DegreeLexicographic]", 0);
}

/* Positive-dimensional ideals, and a first prime that is unlucky: modulo
 * NextPrime[2^62] the leading term p x^2 vanishes, so its leading monomials
 * differ; the next two agree and replace it. */
static void test_modular_route_unlucky_and_positive_dim(void) {
    check_modular_matches_exact("GroebnerBasis[{x^2 + y^2 + z^2 - 1, x y - z}, {x, y, z}, "
        "MonomialOrder -> DegreeReverseLexicographic]", 1);
    check_modular_matches_exact("GroebnerBasis[{x^3 - y z, y^2 - x z + 1}, {x, y, z}, "
        "MonomialOrder -> DegreeReverseLexicographic]", 1);
    check_modular_matches_exact("GroebnerBasis[{NextPrime[2^62] x^2 + y, y^2 + x - 1}, {x, y}, "
        "MonomialOrder -> DegreeReverseLexicographic]", 1);
    check_modular_matches_exact("GroebnerBasis[{NextPrime[2^62] x^2 + y z, y^2 + x - 1, z^2 - x y + 3}, {x, y, z}, "
        "MonomialOrder -> DegreeReverseLexicographic]", 1);
}

/* Lex over Q through the modular route: non-shape zero-dimensional ideals
 * (FGLM mod p per prime), positive-dimensional ones (F4 in Lex mod p), and an
 * unlucky first prime -- each certified by the DRL basis -- against the exact
 * engines. */
static void test_modular_route_lex(void) {
    check_modular_matches_exact("GroebnerBasis[{a + b + c + d, a b + b c + c d + d a, "
                                "a b c + b c d + c d a + d a b, a b c d - 1}, {a, b, c, d}]", 1);
    check_modular_matches_exact("GroebnerBasis[{10 x y^2 + 10 x z^2 - 11 x + 10, 10 y x^2 + 10 y z^2 - 11 y + 10, "
                                "10 z x^2 + 10 z y^2 - 11 z + 10}, {x, y, z}]", 1);
    check_modular_matches_exact("GroebnerBasis[{x^2 - 1, y^2 - 1, z^2 - x y}, {x, y, z}]", 1);
    check_modular_matches_exact("GroebnerBasis[{x^2 + y^2 + z^2 - 1, x y - z}, {x, y, z}]", 1);
    check_modular_matches_exact("GroebnerBasis[{x^3 - y z, y^2 - x z + 1}, {x, y, z}]", 1);
    check_modular_matches_exact("GroebnerBasis[{NextPrime[2^62] x^2 + y^2 - 1, y^2 - 1, z^2 - x}, {x, y, z}]", 1);
}

/* ------------------------------------------------------------------ */

/* InexactNumbers on a zero-dimensional Lex problem goes through the numeric
 * FGLM: the basis must have the exact basis's exponent vectors, and every
 * coefficient must match the exact basis made monic, to the precision
 * asked for (ISSAC '97 problem 6: four dense quadrics, shape basis of
 * degree 16, exact coefficients of ~900 digits). */
static void test_inexact_numeric_fglm_matches_exact(void) {
    const char* setup =
        "eqs6 = {8 w^2 + 5 w x - 4 w y + 2 w z + 3 w + 5 x^2 + 2 x y - 7 x z - 7 x + 7 y^2 - 8 y z - 7 y + 7 z^2 - 8 z + 8, "
        "3 w^2 - 5 w x - 3 w y - 6 w z + 9 w + 4 x^2 + 2 x y - 2 x z + 7 x + 9 y^2 + 6 y z + 5 y + 7 z^2 + 7 z + 5, "
        "-2 w^2 + 9 w x + 9 w y - 7 w z - 4 w + 8 x^2 + 9 x y - 3 x z + 8 x + 6 y^2 - 7 y z + 4 y - 6 z^2 + 8 z + 2, "
        "7 w^2 + 5 w x + 3 w y - 5 w z - 5 w + 2 x^2 + 9 x y - 7 x z + 4 x - 4 y^2 - 5 y z + 6 y - 4 z^2 - 9 z + 2}; "
        "ex6 = GroebnerBasis[eqs6, {x, y, z, w}]; "
        "mon6 = Map[(#/First[CoefficientRules[#, {x, y, z, w}]][[2]]) &, ex6]; True";
    check_true(setup);
    const char* cmp =
        "Module[{g = GroebnerBasis[eqs6, {x, y, z, w}, CoefficientDomain -> %s], "
        "cr, ce}, Length[g] == 4 && "
        "And @@ Table[cr = CoefficientRules[g[[k]], {x, y, z, w}]; "
        "ce = CoefficientRules[mon6[[k]], {x, y, z, w}]; "
        "(cr[[All, 1]] === ce[[All, 1]]) && "
        "Max[Abs[(cr[[All, 2]] - ce[[All, 2]])/ce[[All, 2]]]] < %s, {k, 4}]]";
    char buf[2048];
    snprintf(buf, sizeof buf, cmp, "InexactNumbers", "10^-14");
    check_true(buf);
    snprintf(buf, sizeof buf, cmp, "InexactNumbers[40]", "10^-38");
    check_true(buf);
    check_true("Precision[GroebnerBasis[eqs6, {x, y, z, w}, CoefficientDomain -> InexactNumbers[40]][[1, 1]]] >= 39");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_bivariate_two_polys);
    TEST(test_no_common_roots_three);
    TEST(test_no_common_roots_four);
    TEST(test_trivariate_finite);
    TEST(test_trivariate_infinite);
    TEST(test_univariate_as_gcd_list);
    TEST(test_univariate_as_gcd_symbol);
    TEST(test_linear_gauss);
    TEST(test_eliminate_xy);
    TEST(test_eliminate_two_polys);
    TEST(test_equation_form);
    TEST(test_monomial_order_lex_explicit);
    TEST(test_monomial_order_degrevlex);
    TEST(test_monomial_order_lex_differs);
    TEST(test_empty_input);
    TEST(test_all_zero_input);
    TEST(test_constant_in_ideal);
    TEST(test_single_variable);
    TEST(test_var_order_matters);
    TEST(test_attributes_protected);
    TEST(test_argt_zero);
    TEST(test_argt_one);
    TEST(test_groebner_walk_matches_buchberger);
    TEST(test_walk_matches_buchberger_3var);
    TEST(test_walk_single_variable);
    TEST(test_walk_ideal_one);
    TEST(test_matrix_order_identity_is_lex);
    TEST(test_matrix_order_grevlex_is_degrevlex);
    TEST(test_matrix_order_weighted_no_error);
    TEST(test_matrix_order_walk_matches_direct);
    TEST(test_matrix_order_invalid_falls_back);
    TEST(test_coeff_domain_integers_strong_basis);
    TEST(test_coeff_domain_integers_gcd_combination);
    TEST(test_coeff_domain_integers_documented);
    TEST(test_deglex_honoured);
    TEST(test_non_polynomial_unevaluated);
    TEST(test_coeff_domain_polynomials_equals_default);
    TEST(test_coeff_domain_inexact_numbers);
    TEST(test_real_input_rationalised);
    TEST(test_coeff_domain_rational_functions);
    TEST(test_coeff_domain_rational_functions_no_params);
    TEST(test_coeff_domain_rational_functions_eliminates);
    TEST(test_parametric_coeffs_in_q_of_a);
    TEST(test_membership_input_in_ideal);
    TEST(test_idempotency_no_common_roots);
    TEST(test_two_var_simple_basis);
    TEST(test_subset_vars_auto_parameter);
    TEST(test_sort_true_reverses_vars);
    TEST(test_modulus_gfp_basis);
    TEST(test_modulus_differs_from_rationals);
    TEST(test_modulus_composite_falls_back);
    TEST(test_parameter_variables_auto_main);
    TEST(test_timeconstrained_aborts_buchberger);
    TEST(test_hard_lex_input_by_fglm);
    TEST(test_grevlex_handles_hard_lex_input);
    TEST(test_flint_engine_matches_classical);
    TEST(test_inexact_numeric_fglm_matches_exact);
    TEST(test_f4_matches_buchberger);
    TEST(test_shape_route_random_systems);
    TEST(test_shape_route_unlucky_prime);
    TEST(test_shape_route_declines);
    TEST(test_modular_route_random_systems);
    TEST(test_modular_route_infinity);
    TEST(test_modular_route_unlucky_and_positive_dim);
    TEST(test_modular_route_lex);
    TEST(test_memory_smoke);

    printf("All GroebnerBasis tests passed!\n");
    return 0;
}
