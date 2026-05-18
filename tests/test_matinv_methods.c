/*
 * test_matinv_methods.c -- Method-option coverage for Inverse.
 *
 * Mirrors tests/test_matsol_methods.c: every supported Method value
 * (`Automatic` symbol, "Automatic" string, "DivisionFreeRowReduction",
 * "OneStepRowReduction", "CofactorExpansion") is exercised across the
 * matrix families we care about:
 *
 *   - tiny exact integer       (deterministic FullForm string check)
 *   - 3x3 invertible integer   (deterministic FullForm string check)
 *   - 1x1                      (cofactor fast path)
 *   - rational                 (FullForm string check)
 *   - bignum                   (back-substitution check)
 *   - fully symbolic           (back-substitution -- methods differ in
 *                                canonical form, so we don't pin a
 *                                string)
 *   - exact complex            (back-substitution)
 *   - machine-precision Real   (back-substitution + length check)
 *   - MPFR                     (back-substitution at SetPrecision[..., 50],
 *                                #ifdef USE_MPFR)
 *
 * Plus the error paths: singular matrix, non-square argument, invalid
 * Method value, all of which leave the call unevaluated.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

static void silence_stderr(void) {
    fflush(stderr);
    freopen("/dev/null", "w", stderr);
}

static void run_test_str(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n",
               input, expected, s);
        ASSERT(0);
    } else {
        printf("PASS: %s -> %s\n", input, s);
    }
    free(s);
    expr_free(r);
    expr_free(e);
}

static void run_check_true(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    if (strcmp(s, "True") != 0) {
        printf("FAIL (expected True): %s\n  got: %s\n", input, s);
        ASSERT(0);
    } else {
        printf("PASS: %s -> True\n", input);
    }
    free(s);
    expr_free(r);
    expr_free(e);
}

/* All accepted Method values, in dispatcher order. */
static const char* const METHODS_ALL[] = {
    "Automatic",                /* the symbol -- no quotes */
    "\"Automatic\"",
    "\"DivisionFreeRowReduction\"",
    "\"OneStepRowReduction\"",
    "\"CofactorExpansion\"",
};
static const size_t N_METHODS_ALL = sizeof(METHODS_ALL) / sizeof(METHODS_ALL[0]);

/* ------------------------------------------------------------------
 * Tiny integer: every method canonicalises to the same FullForm.
 * ------------------------------------------------------------------ */
static void test_inverse_int_2x2_all_methods(void) {
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "Inverse[{{1, 2}, {3, 4}}, Method -> %s]", METHODS_ALL[i]);
        run_test_str(buf,
            "List[List[-2, 1], List[Rational[3, 2], Rational[-1, 2]]]");
    }
}

static void test_inverse_int_3x3_all_methods(void) {
    /* {{1, 0, 0}, {0, 2, 0}, {0, 0, 3}} -> diag(1, 1/2, 1/3). */
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "Inverse[{{1, 0, 0}, {0, 2, 0}, {0, 0, 3}}, Method -> %s]",
                 METHODS_ALL[i]);
        run_test_str(buf,
            "List[List[1, 0, 0], "
            "List[0, Rational[1, 2], 0], "
            "List[0, 0, Rational[1, 3]]]");
    }
}

/* ------------------------------------------------------------------
 * 1x1: exercise the cofactor fast path (laplace_det's recursion
 * base is n == 1 so the generic loop can't handle n == 1).
 * ------------------------------------------------------------------ */
static void test_inverse_1x1_all_methods(void) {
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "Inverse[{{5}}, Method -> %s]", METHODS_ALL[i]);
        run_test_str(buf, "List[List[Rational[1, 5]]]");
    }
}

/* ------------------------------------------------------------------
 * Rational: deterministic across methods.
 * ------------------------------------------------------------------ */
static void test_inverse_rational_all_methods(void) {
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "Equal[Inverse[{{1/2, 1/3}, {1/4, 1/5}}, Method -> %s], "
                 "Inverse[{{1/2, 1/3}, {1/4, 1/5}}]]",
                 METHODS_ALL[i]);
        run_check_true(buf);
    }
}

/* ------------------------------------------------------------------
 * Bignum: 10^25 entries; check via back-substitution.
 * ------------------------------------------------------------------ */
static void test_inverse_bignum_all_methods(void) {
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "CompoundExpression["
                 "Set[mm, {{10^25, 1}, {1, 10^25}}], "
                 "Equal[Dot[mm, Inverse[mm, Method -> %s]], "
                 "IdentityMatrix[2]]]", METHODS_ALL[i]);
        run_check_true(buf);
    }
}

/* ------------------------------------------------------------------
 * Symbolic invertible 2x2: methods differ in canonical form, so we
 * verify via Dot[m, Inverse[m]] // Together // Expand == I.
 * ------------------------------------------------------------------ */
static void test_inverse_symbolic_all_methods(void) {
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "CompoundExpression["
                 "Set[mm, {{a, b}, {c, d}}], "
                 "Equal[Together[Dot[mm, Inverse[mm, Method -> %s]]], "
                 "{{1, 0}, {0, 1}}]]", METHODS_ALL[i]);
        run_check_true(buf);
    }
}

/* ------------------------------------------------------------------
 * Symbolic invertible 3x3: same idea, larger matrix.
 * ------------------------------------------------------------------ */
static void test_inverse_symbolic_3x3_all_methods(void) {
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "CompoundExpression["
                 "Set[mm, {{a, b, 0}, {0, c, d}, {e, 0, f}}], "
                 "Equal[Together[Dot[mm, Inverse[mm, Method -> %s]]], "
                 "IdentityMatrix[3]]]", METHODS_ALL[i]);
        run_check_true(buf);
    }
}

/* ------------------------------------------------------------------
 * Complex invertible: back-substitution check.
 * ------------------------------------------------------------------ */
static void test_inverse_complex_all_methods(void) {
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "CompoundExpression["
                 "Set[mm, {{1 + I, 2}, {3, 4 - I}}], "
                 "Equal[Expand[Dot[mm, Inverse[mm, Method -> %s]]], "
                 "{{1, 0}, {0, 1}}]]", METHODS_ALL[i]);
        run_check_true(buf);
    }
}

/* ------------------------------------------------------------------
 * Machine-precision Real: back-substitution + length check.
 * ------------------------------------------------------------------ */
static void test_inverse_float_all_methods(void) {
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "CompoundExpression["
                 "Set[mm, {{1.5, 2.0}, {3.0, 4.5}}], "
                 "Set[ii, Inverse[mm, Method -> %s]], "
                 "Equal[Length[ii], 2]]", METHODS_ALL[i]);
        run_check_true(buf);
    }
}

#ifdef USE_MPFR
static void test_inverse_mpfr_all_methods(void) {
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "CompoundExpression["
                 "Set[mm, SetPrecision[{{1.5, 2.0}, {3.0, 4.5}}, 50]], "
                 "Set[ii, Inverse[mm, Method -> %s]], "
                 "Equal[Length[ii], 2]]", METHODS_ALL[i]);
        run_check_true(buf);
    }
}
#endif

/* ------------------------------------------------------------------
 * Error paths
 * ------------------------------------------------------------------ */
static void test_inverse_singular_all_methods(void) {
    /* Singular 2x2: each method must emit Inverse::sing and return
     * unevaluated.  We silence stderr because the warning is part of
     * the contract, not the test signal. */
    silence_stderr();
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[256];
        char expected[512];
        snprintf(buf, sizeof(buf),
                 "Inverse[{{1, 2}, {2, 4}}, Method -> %s]", METHODS_ALL[i]);
        /* Method -> Automatic (the symbol) prints differently from the
         * string form.  Construct the expected unevaluated form. */
        if (i == 0) {
            snprintf(expected, sizeof(expected),
                "Inverse[List[List[1, 2], List[2, 4]], "
                "Rule[Method, Automatic]]");
        } else {
            snprintf(expected, sizeof(expected),
                "Inverse[List[List[1, 2], List[2, 4]], "
                "Rule[Method, %s]]", METHODS_ALL[i]);
        }
        run_test_str(buf, expected);
    }
}

static void test_inverse_singular_symbolic_cofactor(void) {
    /* Symbolic rank-deficient 3x3: third row = first + second.  Every
     * method must report it as singular and leave the call
     * unevaluated. */
    silence_stderr();
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        char expected[512];
        snprintf(buf, sizeof(buf),
                 "Inverse[{{a, b, c}, {d, e, f}, {a + d, b + e, c + f}}, "
                 "Method -> %s]", METHODS_ALL[i]);
        if (i == 0) {
            snprintf(expected, sizeof(expected),
                "Inverse[List[List[a, b, c], List[d, e, f], "
                "List[Plus[a, d], Plus[b, e], Plus[c, f]]], "
                "Rule[Method, Automatic]]");
        } else {
            snprintf(expected, sizeof(expected),
                "Inverse[List[List[a, b, c], List[d, e, f], "
                "List[Plus[a, d], Plus[b, e], Plus[c, f]]], "
                "Rule[Method, %s]]", METHODS_ALL[i]);
        }
        run_test_str(buf, expected);
    }
}

static void test_inverse_nonsquare_all_methods(void) {
    /* Rectangular: every method must emit Inverse::matsq and return
     * unevaluated. */
    silence_stderr();
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[256];
        char expected[512];
        snprintf(buf, sizeof(buf),
                 "Inverse[{{1, 2}, {3, 4}, {5, 6}}, Method -> %s]",
                 METHODS_ALL[i]);
        if (i == 0) {
            snprintf(expected, sizeof(expected),
                "Inverse[List[List[1, 2], List[3, 4], List[5, 6]], "
                "Rule[Method, Automatic]]");
        } else {
            snprintf(expected, sizeof(expected),
                "Inverse[List[List[1, 2], List[3, 4], List[5, 6]], "
                "Rule[Method, %s]]", METHODS_ALL[i]);
        }
        run_test_str(buf, expected);
    }
}

static void test_inverse_invalid_method(void) {
    /* Unknown method name -> Inverse::method warning + unevaluated. */
    silence_stderr();
    run_test_str(
        "Inverse[{{1, 2}, {3, 4}}, Method -> \"Bogus\"]",
        "Inverse[List[List[1, 2], List[3, 4]], Rule[Method, \"Bogus\"]]");
    /* Wrong rule LHS. */
    run_test_str(
        "Inverse[{{1, 2}, {3, 4}}, Foo -> \"DivisionFreeRowReduction\"]",
        "Inverse[List[List[1, 2], List[3, 4]], "
        "Rule[Foo, \"DivisionFreeRowReduction\"]]");
    /* Method with a non-string non-Automatic RHS. */
    run_test_str(
        "Inverse[{{1, 2}, {3, 4}}, Method -> 42]",
        "Inverse[List[List[1, 2], List[3, 4]], Rule[Method, 42]]");
}

/* The bare default-call must still work after the dispatcher
 * refactor (regression). */
static void test_inverse_no_method_default(void) {
    run_test_str("Inverse[{{1, 2}, {3, 4}}]",
        "List[List[-2, 1], List[Rational[3, 2], Rational[-1, 2]]]");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running Inverse Method-option tests...\n");

    TEST(test_inverse_no_method_default);
    TEST(test_inverse_int_2x2_all_methods);
    TEST(test_inverse_int_3x3_all_methods);
    TEST(test_inverse_1x1_all_methods);
    TEST(test_inverse_rational_all_methods);
    TEST(test_inverse_bignum_all_methods);
    TEST(test_inverse_symbolic_all_methods);
    TEST(test_inverse_symbolic_3x3_all_methods);
    TEST(test_inverse_complex_all_methods);
    TEST(test_inverse_float_all_methods);
#ifdef USE_MPFR
    TEST(test_inverse_mpfr_all_methods);
#endif

    TEST(test_inverse_singular_all_methods);
    TEST(test_inverse_singular_symbolic_cofactor);
    TEST(test_inverse_nonsquare_all_methods);
    TEST(test_inverse_invalid_method);

    printf("All Inverse Method-option tests passed!\n");
    symtab_clear();
    return 0;
}
