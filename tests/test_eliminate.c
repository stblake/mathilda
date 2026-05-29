/* Unit tests for Eliminate (src/poly/eliminate.c).
 *
 * Eliminate is a thin driver over the lex-order Buchberger engine in
 * src/poly/groebner.c with an elimination block, plus a principal-branch
 * inverse-function pre-pass and a balanced `Equal[lhs, rhs]` reconstruction.
 *
 * The expected outputs are taken from a smoke run of the implementation
 * (cross-checked for mathematical equivalence with Mathematica's reference
 * answers for the documented examples).  Matching FullForm byte-for-byte
 * is the strongest correctness signal we have.
 *
 * Run binary directly: ./eliminate_tests
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

/* "expr evaluates to True" assertion (avoids hard-coded long FullForm). */
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

/* Silence stderr once for diagnostic-emitting tests. */
static void mute_stderr_once(void) {
    static int done = 0;
    if (!done) { freopen("/dev/null", "w", stderr); done = 1; }
}

/* ------------------------------------------------------------------ */
/* 1. Linear pair: headline example from the prompt                    */
/* ------------------------------------------------------------------ */

static void test_linear_pair(void) {
    /* Eliminate[{x == 2+y, y == z}, y] -> 2 + z == x */
    check_eq("Eliminate[{x == 2+y, y == z}, y]",
             "Equal[Plus[2, z], x]");
}

/* ------------------------------------------------------------------ */
/* 2. Linear system, eliminate z                                       */
/* ------------------------------------------------------------------ */

static void test_linear_system_eliminate_one(void) {
    /* Eliminate[{2x+3y+4z == 1, 9x+8y+7z == 2}, z]
     * Mathematica: 1 - 11 y == 22 x.  Sign-equivalent FullForm below. */
    check_eq("Eliminate[{2x+3y+4z == 1, 9x+8y+7z == 2}, z]",
             "Equal[Plus[Times[22, x], Times[11, y]], 1]");
}

/* ------------------------------------------------------------------ */
/* 3. Symmetric-function elimination                                   */
/* ------------------------------------------------------------------ */

static void test_symmetric_functions(void) {
    /* Eliminate[{f == x^5+y^5, a == x+y, b == x*y}, {x, y}]
     * Mathematica: f == a^5 - 5 a^3 b + 5 a b^2 (sign-balanced equiv). */
    check_eq("Eliminate[{f == x^5+y^5, a == x+y, b == x*y}, {x, y}]",
             "Equal[Plus[Power[a, 5], Times[5, Times[a, Power[b, 2]]]], "
                   "Plus[Times[5, Times[Power[a, 3], b]], f]]");
}

/* ------------------------------------------------------------------ */
/* 4. Polynomial system, eliminate z -> two equations remain           */
/* ------------------------------------------------------------------ */

static void test_poly_system_eliminate_z(void) {
    /* Eliminate[{x^2+y^2+z^2 == 1, x-y+z == 2, x^3-y^2 == z+1}, z]
     * One purely-polynomial constraint on x plus a y-eq parametric in x. */
    check_eq("Eliminate[{x^2+y^2+z^2 == 1, x-y+z == 2, x^3-y^2 == z+1}, z]",
             "And[Equal[Plus[27, Times[4, Power[x, 2]], Times[8, Power[x, 4]], "
                            "Times[4, Power[x, 5]], Times[4, Power[x, 6]]], "
                       "Plus[Times[18, x], Times[28, Power[x, 3]]]], "
                  "Equal[Plus[12, Times[2, x], Times[5, Power[x, 2]], y], "
                       "Plus[Times[8, Power[x, 3]], Times[4, Power[x, 4]], "
                            "Times[2, Power[x, 5]]]]]");
}

/* ------------------------------------------------------------------ */
/* 5. Same system, two-variable elimination -> single equation         */
/* ------------------------------------------------------------------ */

static void test_poly_system_eliminate_y_z(void) {
    check_eq("Eliminate[{x^2+y^2+z^2 == 1, x-y+z == 2, x^3-y^2 == z+1}, {y, z}]",
             "Equal[Plus[27, Times[4, Power[x, 2]], Times[8, Power[x, 4]], "
                        "Times[4, Power[x, 5]], Times[4, Power[x, 6]]], "
                   "Plus[Times[18, x], Times[28, Power[x, 3]]]]");
}

/* ------------------------------------------------------------------ */
/* 6. && form input is equivalent to {} form                           */
/* ------------------------------------------------------------------ */

static void test_and_form_equivalent(void) {
    check_true(
        "Eliminate[(2x+3y+4z == 1) && (9x+8y+7z == 2), z] === "
        "Eliminate[{2x+3y+4z == 1, 9x+8y+7z == 2}, z]");
}

/* ------------------------------------------------------------------ */
/* 7. Single-symbol var shorthand: y vs {y}                            */
/* ------------------------------------------------------------------ */

static void test_single_symbol_var_shorthand(void) {
    check_true(
        "Eliminate[{x == 2+y, y == z}, y] === "
        "Eliminate[{x == 2+y, y == z}, {y}]");
}

/* ------------------------------------------------------------------ */
/* 8. Resultant cross-check: incompatible roots -> False               */
/* ------------------------------------------------------------------ */

static void test_resultant_incompatible(void) {
    /* x^2 - 1 = 0 and x^2 - 4 = 0 share no root. */
    check_eq("Eliminate[{x^2 - 1 == 0, x^2 - 4 == 0}, x]",
             "False");
}

/* ------------------------------------------------------------------ */
/* 9. Common-root condition: a == b                                    */
/* ------------------------------------------------------------------ */

static void test_common_root_condition(void) {
    /* x - a = 0 and x - b = 0 share a root iff a == b. */
    check_eq("Eliminate[{x - a == 0, x - b == 0}, x]",
             "Equal[b, a]");
}

/* ------------------------------------------------------------------ */
/* 10. Empty elimination ideal -> True                                 */
/* ------------------------------------------------------------------ */

static void test_empty_ideal_returns_true(void) {
    /* One equation, one unknown to eliminate -> y solvable for any x. */
    check_eq("Eliminate[x + y == 0, y]", "True");
}

/* ------------------------------------------------------------------ */
/* 11. Numeric contradiction -> False                                  */
/* ------------------------------------------------------------------ */

static void test_contradiction_returns_false(void) {
    /* `1 == 2` evaluates to `False` before Eliminate runs; we fold that
     * sentinel rather than emit Eliminate::eqf. */
    check_eq("Eliminate[1 == 2, x]", "False");
}

/* ------------------------------------------------------------------ */
/* 12. Tautology -> True                                               */
/* ------------------------------------------------------------------ */

static void test_tautology_returns_true(void) {
    /* `x == x` evaluates to `True`; same sentinel-folding path. */
    check_eq("Eliminate[x == x, x]", "True");
}

/* ------------------------------------------------------------------ */
/* 13. Transcendental pre-pass: principal-branch inverse rewrite       */
/* ------------------------------------------------------------------ */

static void test_inverse_function_prepass(void) {
    /* `Eliminate[Sin[y] == 1, y]` peels Sin via ArcSin (ifun fires),
     * leaves y == Pi/2; eliminating y from a single solvable equation
     * yields no residual constraint -> True. */
    mute_stderr_once();
    check_eq("Eliminate[Sin[y] == 1, y]", "True");
}

static void test_inverse_function_two_vars(void) {
    /* Sin[x + y] == 1 -> x + y == Pi/2, then eliminating y is trivial. */
    mute_stderr_once();
    check_eq("Eliminate[Sin[x + y] == 1, y]", "True");
}

/* ------------------------------------------------------------------ */
/* 14. Empty variable list -> equations pass through                   */
/* ------------------------------------------------------------------ */

static void test_empty_var_list(void) {
    /* No variables to eliminate -> just return the equation unchanged. */
    check_eq("Eliminate[{x == y + 1}, {}]", "Equal[x, Plus[1, y]]");
}

/* ------------------------------------------------------------------ */
/* 15. Attribute check: Eliminate is Protected                         */
/* ------------------------------------------------------------------ */

static void test_eliminate_protected(void) {
    /* Setting Eliminate[x_] := ... must fail on a Protected head. */
    mute_stderr_once();
    check_true("MemberQ[Attributes[Eliminate], Protected]");
}

/* ------------------------------------------------------------------ */
/* 16. Memory smoke -- run the symmetric-function spec example         */
/* ------------------------------------------------------------------ */

static void test_memory_smoke(void) {
    for (int i = 0; i < 25; i++) {
        Expr* e = parse_expression(
            "Eliminate[{f == x^5+y^5, a == x+y, b == x*y}, {x, y}]");
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(e);
        expr_free(r);
    }
}

/* ------------------------------------------------------------------ */
/* 17. Memory smoke -- larger polynomial system                        */
/* ------------------------------------------------------------------ */

static void test_memory_smoke_larger(void) {
    for (int i = 0; i < 25; i++) {
        Expr* e = parse_expression(
            "Eliminate[{x^2+y^2+z^2 == 1, x-y+z == 2, x^3-y^2 == z+1}, "
            "{y, z}]");
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(e);
        expr_free(r);
    }
}

/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_linear_pair);
    TEST(test_linear_system_eliminate_one);
    TEST(test_symmetric_functions);
    TEST(test_poly_system_eliminate_z);
    TEST(test_poly_system_eliminate_y_z);
    TEST(test_and_form_equivalent);
    TEST(test_single_symbol_var_shorthand);
    TEST(test_resultant_incompatible);
    TEST(test_common_root_condition);
    TEST(test_empty_ideal_returns_true);
    TEST(test_contradiction_returns_false);
    TEST(test_tautology_returns_true);
    TEST(test_inverse_function_prepass);
    TEST(test_inverse_function_two_vars);
    TEST(test_empty_var_list);
    TEST(test_eliminate_protected);
    TEST(test_memory_smoke);
    TEST(test_memory_smoke_larger);

    printf("All Eliminate tests passed!\n");
    return 0;
}
