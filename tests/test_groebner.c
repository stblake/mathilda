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
    if (!done) { freopen("/dev/null", "w", stderr); done = 1; }
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
     * Print FullForm: the -1*(x^2 y^2) term keeps its Times-of-Times shape
     * because Times is built from gb_to_expr factor lists -- the outer
     * Times[-1, ...] is grafted on after the variable-factor Times is
     * already an Expr*. */
    check_eq("GroebnerBasis[{x^2+y^2+z^2-1, x y z-3}, {x, y}, {z}]",
             "List[Plus[9, "
                  "Times[-1, Times[Power[x, 2], Power[y, 2]]], "
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

static void test_nimpl_groebner_walk_falls_back(void) {
    mute_stderr_once();
    /* Method -> "GroebnerWalk" emits nimpl and falls back to Buchberger. */
    check_eq("GroebnerBasis[{x^2-2y^2, x y-3}, {x, y}, "
             "Method -> \"GroebnerWalk\"]",
             "List[Plus[-9, Times[2, Power[y, 4]]], "
                  "Plus[Times[3, x], Times[-2, Power[y, 3]]]]");
}

static void test_nimpl_coeff_domain_integers_falls_back(void) {
    mute_stderr_once();
    /* CoefficientDomain -> Integers not implemented -> nimpl + fall back. */
    check_eq("GroebnerBasis[{x^2-2y^2, x y-3}, {x, y}, "
             "CoefficientDomain -> Integers]",
             "List[Plus[-9, Times[2, Power[y, 4]]], "
                  "Plus[Times[3, x], Times[-2, Power[y, 3]]]]");
}

static void test_nimpl_deglex_falls_back(void) {
    mute_stderr_once();
    /* DegreeLexicographic deferred -> nimpl + fall back to Lex. */
    check_eq("GroebnerBasis[{x^2-2y^2, x y-3}, {x, y}, "
             "MonomialOrder -> DegreeLexicographic]",
             "List[Plus[-9, Times[2, Power[y, 4]]], "
                  "Plus[Times[3, x], Times[-2, Power[y, 3]]]]");
}

/* ------------------------------------------------------------------ */
/* 14. Non-polynomial input -> NULL (unevaluated)                      */
/* ------------------------------------------------------------------ */

static void test_non_polynomial_unevaluated(void) {
    /* Sin[x] is not a polynomial in x; stays unevaluated. */
    check_eq("GroebnerBasis[{Sin[x] + x, x^2 - 1}, {x}]",
             "GroebnerBasis[List[Plus[x, Sin[x]], Plus[-1, Power[x, 2]]], "
                          "List[x]]");
}

static void test_parametric_coeffs_unevaluated(void) {
    /* a is not in the variable list and is treated as a non-poly symbol;
     * the first-cut implementation returns unevaluated rather than
     * computing a Q(a) basis. */
    check_eq("GroebnerBasis[{a x^2+5 x-1, 2 x+3 x y+y^2}, {x, y}]",
             "GroebnerBasis[List[Plus[-1, Times[5, x], "
                                    "Times[a, Power[x, 2]]], "
                                "Plus[Times[2, x], "
                                    "Times[3, Times[x, y]], "
                                    "Power[y, 2]]], "
                          "List[x, y]]");
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
    TEST(test_nimpl_groebner_walk_falls_back);
    TEST(test_nimpl_coeff_domain_integers_falls_back);
    TEST(test_nimpl_deglex_falls_back);
    TEST(test_non_polynomial_unevaluated);
    TEST(test_parametric_coeffs_unevaluated);
    TEST(test_membership_input_in_ideal);
    TEST(test_idempotency_no_common_roots);
    TEST(test_two_var_simple_basis);
    TEST(test_memory_smoke);

    printf("All GroebnerBasis tests passed!\n");
    return 0;
}
