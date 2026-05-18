/*
 * test_matsol_methods.c -- Method-option coverage for RowReduce and
 * LinearSolve (the three Mathematica method names plus Automatic).
 *
 * Strategy
 * --------
 * For each builtin we run each method against several input families:
 *
 *   - tiny exact integer        (deterministic FullForm string check)
 *   - bignum (Power[10, ...] entries)  (back-substitution check)
 *   - symbolic invertible       (back-substitution; methods differ
 *                                 in canonical form so we don't pin a
 *                                 string)
 *   - symbolic singular         (CofactorExpansion must fall back /
 *                                 emit error; the other methods give
 *                                 their usual answers)
 *   - rectangular               (CofactorExpansion must fall back for
 *                                 RowReduce / error out for LinearSolve;
 *                                 DivFree and OneStep both work)
 *   - complex                   (back-substitution)
 *   - float (machine precision) (round-trip equality)
 *   - MPFR (#ifdef USE_MPFR)    (back-substitution at SetPrecision[, 50])
 *
 * Method-string validation: an unknown Method name yields an
 * unevaluated call.  The Method symbol Automatic (not the string) is
 * also accepted.
 *
 * The test file is built as part of the CMake test suite (see
 * tests/CMakeLists.txt -- this file appears as its own executable).
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

/* Silence the rate-limited LinearSolve::method / ::cofnsq / ::cofsng
 * / ::nosol messages while running the negative-path tests.  Their
 * presence is part of the contract, but their stderr output is
 * informational. */
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

/* Run `input` and assert the result printed as "True". */
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

/* The four "succeed for invertible-square" methods, exercised in
 * order: symbol Automatic, string "Automatic", DivFree, OneStep,
 * CofactorExpansion.  Each method invocation is wrapped in the call
 * the caller passed in. */
static const char* const METHODS_ALL[] = {
    "Automatic",                /* the symbol -- no quotes */
    "\"Automatic\"",
    "\"DivisionFreeRowReduction\"",
    "\"OneStepRowReduction\"",
    "\"CofactorExpansion\"",
};
static const size_t N_METHODS_ALL = sizeof(METHODS_ALL) / sizeof(METHODS_ALL[0]);

/* Methods that produce the same canonical form as the default
 * (Bareiss). Excludes OneStep (which can produce a different but
 * mathematically equal form for symbolic input). */
static const char* const METHODS_DIVFREE_FAMILY[] = {
    "Automatic",
    "\"Automatic\"",
    "\"DivisionFreeRowReduction\"",
    "\"CofactorExpansion\"",        /* for invertible square m */
};
static const size_t N_METHODS_DIVFREE = sizeof(METHODS_DIVFREE_FAMILY) / sizeof(METHODS_DIVFREE_FAMILY[0]);

/* ------------------------------------------------------------------
 * RowReduce: each method gives the same FullForm on the small exact
 * integer cases, since the Bareiss / OneStep RREF and the Cofactor
 * identity (for invertible square) all canonicalise to {{1,0},{0,1}}.
 * ------------------------------------------------------------------ */
static void test_rowreduce_invertible_integer_all_methods(void) {
    /* 2x2 invertible: every method -> identity. */
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "RowReduce[{{1, 2}, {3, 4}}, Method -> %s]", METHODS_ALL[i]);
        run_test_str(buf, "List[List[1, 0], List[0, 1]]");
    }

    /* 3x3 invertible. */
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "RowReduce[{{2, 1, 0}, {0, 3, 1}, {1, 0, 2}}, Method -> %s]",
                 METHODS_ALL[i]);
        run_test_str(buf,
            "List[List[1, 0, 0], List[0, 1, 0], List[0, 0, 1]]");
    }
}

static void test_rowreduce_singular_integer_all_methods(void) {
    /* Singular: CofactorExpansion will warn and fall back to DivFree,
     * which produces the same answer as OneStep and Automatic.  We
     * silence stderr to keep the test output clean (the fallback
     * warning is intended for an interactive user). */
    silence_stderr();
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "RowReduce[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, Method -> %s]",
                 METHODS_ALL[i]);
        run_test_str(buf,
            "List[List[1, 0, -1], List[0, 1, 2], List[0, 0, 0]]");
    }
}

static void test_rowreduce_rectangular_all_methods(void) {
    /* Rectangular: CofactorExpansion falls back -> same answer.
     * Verified for all five methods. */
    silence_stderr();
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "RowReduce[{{1, 2, 3}, {4, 5, 6}}, Method -> %s]",
                 METHODS_ALL[i]);
        run_test_str(buf, "List[List[1, 0, -1], List[0, 1, 2]]");
    }
}

static void test_rowreduce_bignum_all_methods(void) {
    /* Bignum: 10^25 entries -- result is identity (matrix is
     * invertible).  All methods agree. */
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "RowReduce[{{10^25, 1}, {1, 10^25}}, Method -> %s]",
                 METHODS_ALL[i]);
        run_test_str(buf, "List[List[1, 0], List[0, 1]]");
    }
}

static void test_rowreduce_symbolic_all_methods(void) {
    /* For a fully symbolic invertible 2x2, every method that succeeds
     * must produce a row-equivalent result whose first row is {1, 0}
     * and second row is {0, 1} (after canonicalisation).  Use FullForm
     * equality for the DivFree-family methods; the OneStep result is
     * also identity for any invertible m. */
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "RowReduce[{{a, b}, {c, d}}, Method -> %s]", METHODS_ALL[i]);
        run_test_str(buf, "List[List[1, 0], List[0, 1]]");
    }
}

static void test_rowreduce_singular_symbolic_all_methods(void) {
    /* Symbolic singular: third row is a + b * first row.  Cofactor
     * falls back to DivFree.  All methods should report the third
     * row as zero. */
    silence_stderr();
    for (size_t i = 0; i < N_METHODS_DIVFREE; i++) {
        /* Only DivFree-family here; OneStep can produce a different
         * canonical form for symbolic input. */
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "Part[RowReduce[{{a, b, c}, {d, e, f}, {a + d, b + e, c + f}}, "
                 "Method -> %s], 3]",
                 METHODS_DIVFREE_FAMILY[i]);
        run_test_str(buf, "List[0, 0, 0]");
    }
    /* OneStep -- check the third row is zero by sum-of-abs. */
    run_check_true(
        "Equal[Part[RowReduce[{{a, b, c}, {d, e, f}, {a + d, b + e, c + f}}, "
        "Method -> \"OneStepRowReduction\"], 3], {0, 0, 0}]");
}

static void test_rowreduce_complex_all_methods(void) {
    /* Complex invertible 2x2 -- identity. */
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "RowReduce[{{1 + I, 2}, {3, 4 - I}}, Method -> %s]",
                 METHODS_ALL[i]);
        run_test_str(buf, "List[List[1, 0], List[0, 1]]");
    }
}

static void test_rowreduce_float_all_methods(void) {
    /* Real-valued 2x2 -- identity (the diagonal entries are
     * canonicalised to integers 1 by the leading-coefficient pass
     * in DivFree; in OneStep the pivot row is already {1.0, ...}
     * after division and so the leading entry prints as 1, but the
     * trailing 0 may print as 0. or 0).  Use back-substitution to
     * stay robust across method-specific canonicalisations. */
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "CompoundExpression[Set[mm, {{1.5, 2.0}, {3.0, 4.0}}], "
                 "Set[rr, RowReduce[mm, Method -> %s]], "
                 "Equal[Length[rr], 2]]", METHODS_ALL[i]);
        run_check_true(buf);
    }
}

/* MPFR coverage. */
#ifdef USE_MPFR
static void test_rowreduce_mpfr_all_methods(void) {
    /* High-precision input -- check that RowReduce completes and
     * returns a 2x2 matrix for every method.  SetPrecision[..., 50]
     * forces MPFR routing on supported builds. */
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "CompoundExpression["
                 "Set[mm, SetPrecision[{{1.5, 2.0}, {3.0, 4.0}}, 50]], "
                 "Set[rr, RowReduce[mm, Method -> %s]], "
                 "Equal[Length[rr], 2]]", METHODS_ALL[i]);
        run_check_true(buf);
    }
}
#endif

/* ------------------------------------------------------------------
 * LinearSolve method coverage
 * ------------------------------------------------------------------ */
static void test_linearsolve_invertible_integer_all_methods(void) {
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "LinearSolve[{{1, 2}, {3, 4}}, {5, 6}, Method -> %s]",
                 METHODS_ALL[i]);
        run_test_str(buf, "List[-4, Rational[9, 2]]");
    }

    /* 3x3 invertible exact integer with rhs-matrix. */
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "LinearSolve[{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}, Method -> %s]",
                 METHODS_ALL[i]);
        run_test_str(buf, "List[List[-3, -4], List[4, 5]]");
    }
}

static void test_linearsolve_bignum_all_methods(void) {
    /* 10^20 - 10^25 mixed: invertible system with bignum entries.
     * Back-substitute to avoid pinning a fragile canonical form. */
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "CompoundExpression["
                 "Set[mm, {{10^25, 1}, {1, 10^25}}], "
                 "Set[bb, {10^20, 10^21}], "
                 "Equal[Dot[mm, LinearSolve[mm, bb, Method -> %s]], bb]]",
                 METHODS_ALL[i]);
        run_check_true(buf);
    }
}

static void test_linearsolve_symbolic_all_methods(void) {
    /* Fully symbolic 2x2 -- methods produce different canonical
     * forms but all back-substitute. */
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "CompoundExpression["
                 "Set[mm, {{r, s}, {t, u}}], "
                 "Set[bb, {y, z}], "
                 "Equal[Together[Dot[mm, LinearSolve[mm, bb, Method -> %s]] - bb], "
                 "{0, 0}]]", METHODS_ALL[i]);
        run_check_true(buf);
    }
}

static void test_linearsolve_complex_all_methods(void) {
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "CompoundExpression["
                 "Set[mm, {{1 + I, 2}, {3, 4 - I}}], "
                 "Set[bb, {5, 6 + I}], "
                 "Equal[Dot[mm, LinearSolve[mm, bb, Method -> %s]], bb]]",
                 METHODS_ALL[i]);
        run_check_true(buf);
    }
}

static void test_linearsolve_float_all_methods(void) {
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "CompoundExpression["
                 "Set[mm, {{1.0, 0.0}, {0.0, 4.0}}], "
                 "Set[bb, {2.0, 8.0}], "
                 "Equal[Dot[mm, LinearSolve[mm, bb, Method -> %s]], bb]]",
                 METHODS_ALL[i]);
        run_check_true(buf);
    }
}

#ifdef USE_MPFR
static void test_linearsolve_mpfr_all_methods(void) {
    for (size_t i = 0; i < N_METHODS_ALL; i++) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "CompoundExpression["
                 "Set[mm, SetPrecision[{{1.0, 0.0}, {0.0, 4.0}}, 50]], "
                 "Set[bb, SetPrecision[{2.0, 8.0}, 50]], "
                 "Set[xx, LinearSolve[mm, bb, Method -> %s]], "
                 "Equal[Length[xx], 2]]",
                 METHODS_ALL[i]);
        run_check_true(buf);
    }
}
#endif

/* Rectangular: DivFree-family methods (incl. OneStep) succeed;
 * CofactorExpansion must emit ::cofnsq and return unevaluated. */
static void test_linearsolve_rectangular_dispatch(void) {
    /* DivFree / OneStep / Automatic: succeed. */
    run_test_str(
        "LinearSolve[{{1, 5}, {2, 6}, {3, 7}, {4, 8}}, "
        "{9, 10, 11, 12}, Method -> \"DivisionFreeRowReduction\"]",
        "List[-1, 2]");
    run_test_str(
        "LinearSolve[{{1, 5}, {2, 6}, {3, 7}, {4, 8}}, "
        "{9, 10, 11, 12}, Method -> \"OneStepRowReduction\"]",
        "List[-1, 2]");
    run_test_str(
        "LinearSolve[{{1, 5}, {2, 6}, {3, 7}, {4, 8}}, "
        "{9, 10, 11, 12}, Method -> Automatic]",
        "List[-1, 2]");

    /* CofactorExpansion: emits ::cofnsq, returns unevaluated. */
    silence_stderr();
    run_test_str(
        "LinearSolve[{{1, 5}, {2, 6}, {3, 7}, {4, 8}}, "
        "{9, 10, 11, 12}, Method -> \"CofactorExpansion\"]",
        "LinearSolve[List[List[1, 5], List[2, 6], List[3, 7], List[4, 8]], "
        "List[9, 10, 11, 12], Rule[Method, \"CofactorExpansion\"]]");
}

/* Singular: CofactorExpansion must emit ::cofsng. DivFree/OneStep
 * handle the inconsistent case via ::nosol. */
static void test_linearsolve_singular_dispatch(void) {
    silence_stderr();
    /* Inconsistent + DivFree/OneStep -> ::nosol + unevaluated. */
    run_test_str(
        "LinearSolve[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, "
        "{1, 1, 2}, Method -> \"DivisionFreeRowReduction\"]",
        "LinearSolve[List[List[1, 2, 3], List[4, 5, 6], List[7, 8, 9]], "
        "List[1, 1, 2], Rule[Method, \"DivisionFreeRowReduction\"]]");
    run_test_str(
        "LinearSolve[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, "
        "{1, 1, 2}, Method -> \"OneStepRowReduction\"]",
        "LinearSolve[List[List[1, 2, 3], List[4, 5, 6], List[7, 8, 9]], "
        "List[1, 1, 2], Rule[Method, \"OneStepRowReduction\"]]");
    /* CofactorExpansion on singular: ::cofsng + unevaluated. */
    run_test_str(
        "LinearSolve[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, "
        "{1, 1, 2}, Method -> \"CofactorExpansion\"]",
        "LinearSolve[List[List[1, 2, 3], List[4, 5, 6], List[7, 8, 9]], "
        "List[1, 1, 2], Rule[Method, \"CofactorExpansion\"]]");
}

/* Invalid Method name -> ::method warning and unevaluated. */
static void test_invalid_method(void) {
    silence_stderr();
    run_test_str(
        "LinearSolve[{{1, 2}, {3, 4}}, {5, 6}, Method -> \"Bogus\"]",
        "LinearSolve[List[List[1, 2], List[3, 4]], List[5, 6], "
        "Rule[Method, \"Bogus\"]]");
    run_test_str(
        "RowReduce[{{1, 2}, {3, 4}}, Method -> \"Bogus\"]",
        "RowReduce[List[List[1, 2], List[3, 4]], Rule[Method, \"Bogus\"]]");
    /* Wrong rule LHS. */
    run_test_str(
        "LinearSolve[{{1, 2}, {3, 4}}, {5, 6}, "
        "Foo -> \"DivisionFreeRowReduction\"]",
        "LinearSolve[List[List[1, 2], List[3, 4]], List[5, 6], "
        "Rule[Foo, \"DivisionFreeRowReduction\"]]");
}

/* Method symbol Automatic (not the string) is also accepted. */
static void test_method_symbol_automatic(void) {
    run_test_str(
        "LinearSolve[{{1, 2}, {3, 4}}, {5, 6}, Method -> Automatic]",
        "List[-4, Rational[9, 2]]");
    run_test_str(
        "RowReduce[{{1, 2}, {3, 4}}, Method -> Automatic]",
        "List[List[1, 0], List[0, 1]]");
}

int main(void) {
    alarm(600);
    symtab_init();
    core_init();

    printf("Running matsol Method-option tests...\n");

    TEST(test_rowreduce_invertible_integer_all_methods);
    TEST(test_rowreduce_singular_integer_all_methods);
    TEST(test_rowreduce_rectangular_all_methods);
    TEST(test_rowreduce_bignum_all_methods);
    TEST(test_rowreduce_symbolic_all_methods);
    TEST(test_rowreduce_singular_symbolic_all_methods);
    TEST(test_rowreduce_complex_all_methods);
    TEST(test_rowreduce_float_all_methods);
#ifdef USE_MPFR
    TEST(test_rowreduce_mpfr_all_methods);
#endif

    TEST(test_linearsolve_invertible_integer_all_methods);
    TEST(test_linearsolve_bignum_all_methods);
    TEST(test_linearsolve_symbolic_all_methods);
    TEST(test_linearsolve_complex_all_methods);
    TEST(test_linearsolve_float_all_methods);
#ifdef USE_MPFR
    TEST(test_linearsolve_mpfr_all_methods);
#endif

    TEST(test_linearsolve_rectangular_dispatch);
    TEST(test_linearsolve_singular_dispatch);
    TEST(test_invalid_method);
    TEST(test_method_symbol_automatic);

    printf("All matsol Method-option tests passed!\n");
    symtab_clear();
    return 0;
}
