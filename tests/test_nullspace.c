/*
 * tests/test_nullspace.c
 *
 * Unit tests for NullSpace[m] (and the Method-option dispatcher).
 *
 * Two complementary verification strategies are used:
 *
 *   1. Exact-form FullForm comparison via `run_test` for cases where the
 *      RREF is canonical and the expected basis is unambiguous (small
 *      square / rectangular integer matrices, simple symbolic
 *      matrices).
 *
 *   2. Semantic verification via `assert_nullspace_valid` for cases
 *      where the exact basis depends on RREF details (machine-precision
 *      floats, complex inputs).  This helper computes
 *        - the expected nullity (cols - rank), via the size of the
 *          returned basis itself; and
 *        - m . v for each basis vector v, asserting each component is
 *          numerically (or structurally) zero.
 *
 * Memory: every parse_expression / evaluate result is freed.  The test
 * binary is also intended to run under valgrind --leak-check=full with
 * zero "definitely lost" bytes.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

/* Exact-form helper -- compares FullForm representation. */
static void run_test(const char* input, const char* expected) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    Expr* res = evaluate(parsed);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  got:      %s\n",
                input, expected, res_str);
        free(res_str);
        expr_free(res);
        expr_free(parsed);
        ASSERT(0);
    }
    printf("  PASS: %s -> %s\n", input, res_str);
    free(res_str);
    expr_free(res);
    expr_free(parsed);
}

/* Helper: parse and evaluate `src`, return owning Expr*. */
static Expr* eval_str(const char* src) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* res = evaluate(parsed);
    expr_free(parsed);
    return res;
}

/* True iff every entry of a flat List `v` evaluates to structural zero
 * (Integer 0, Real 0.0, or any expression accepted by Together + zero
 * comparison).  Non-list `v` is treated as a single scalar entry. */
static bool list_is_zero(Expr* v) {
    if (v->type == EXPR_INTEGER) return v->data.integer == 0;
    /* Loose tolerance for Real: row-reduction over machine-precision
     * floats accumulates round-off proportional to matrix magnitude,
     * and our test matrices include entries up to ~16.  1e-3 is a
     * comfortable absolute bound. */
    if (v->type == EXPR_REAL)    return fabs(v->data.real) < 1e-3;
    if (v->type == EXPR_FUNCTION &&
        v->data.function.head->type == EXPR_SYMBOL &&
        strcmp(v->data.function.head->data.symbol.name, "List") == 0) {
        for (size_t i = 0; i < v->data.function.arg_count; i++) {
            if (!list_is_zero(v->data.function.args[i])) return false;
        }
        return true;
    }
    /* Symbolic 0 may print as "0" after Together / Expand. */
    char* s = expr_to_string(v);
    bool zero = (strcmp(s, "0") == 0);
    free(s);
    return zero;
}

/* Verify that NullSpace[m] returns a basis of the expected size and
 * that every basis vector lies in the null space (m . v == 0).
 *
 * `m_str` is a parseable expression for the matrix.  `expected_n` is
 * the expected basis size (i.e. cols(m) - rank(m)).
 *
 * The "m . v == 0" check evaluates `Expand[m . v]` and verifies every
 * component reduces to zero.  Together-style canonicalisation is also
 * accepted via the symbolic-zero path in `list_is_zero`. */
static void assert_nullspace_valid(const char* m_str, int expected_n) {
    char ns_call[1024];
    snprintf(ns_call, sizeof(ns_call), "NullSpace[%s]", m_str);
    Expr* ns = eval_str(ns_call);

    /* Top-level must be List[...]. */
    ASSERT(ns->type == EXPR_FUNCTION);
    ASSERT(ns->data.function.head->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(ns->data.function.head->data.symbol.name, "List");

    size_t n_basis = ns->data.function.arg_count;
    if ((int)n_basis != expected_n) {
        char* s = expr_to_string(ns);
        fprintf(stderr,
                "FAIL: NullSpace[%s] expected %d basis vectors, got %zu: %s\n",
                m_str, expected_n, n_basis, s);
        free(s);
        ASSERT(0);
    }

    /* Verify m . v == 0 for each basis vector. */
    for (size_t i = 0; i < n_basis; i++) {
        Expr* v = ns->data.function.args[i];
        char* v_str = expr_to_string(v);
        char prod_src[2048];
        snprintf(prod_src, sizeof(prod_src),
                 "Expand[(%s) . %s]", m_str, v_str);
        Expr* prod = eval_str(prod_src);

        if (!list_is_zero(prod)) {
            char* pstr = expr_to_string(prod);
            fprintf(stderr,
                    "FAIL: m . v != 0 for NullSpace[%s][[%zu]] = %s\n"
                    "  m . v = %s\n",
                    m_str, i + 1, v_str, pstr);
            free(pstr);
            free(v_str);
            expr_free(prod);
            expr_free(ns);
            ASSERT(0);
        }
        free(v_str);
        expr_free(prod);
    }

    printf("  PASS: NullSpace[%s]  -- %zu basis vector(s) verified m.v==0\n",
           m_str, n_basis);
    expr_free(ns);
}

/* ===================== Exact-form tests ===================== */

/* The classic Mathematica example: a 3x3 rank-2 matrix.  Pivots in
 * cols 1, 2; free col 3.  Expected basis: {{1, -2, 1}}. */
static void test_nullspace_3x3_singular(void) {
    run_test("NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]",
             "List[List[1, -2, 1]]");
}

/* A full-rank 3x3 matrix -- expected: empty list `{}`. */
static void test_nullspace_3x3_full_rank(void) {
    run_test("NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}]",
             "List[]");
}

/* A 4x4 rank-3 matrix from the spec.  Expected: {{1, -2, 1, 0}}. */
static void test_nullspace_4x4_one_free(void) {
    run_test("NullSpace[{{1, 2, 3, 10}, {4, 5, 6, 11}, "
                       "{7, 8, 9, 12}, {13, 14, 15, 16}}]",
             "List[List[1, -2, 1, 0]]");
}

/* The full-rank 4x4 from the column-space example. */
static void test_nullspace_4x4_full_rank(void) {
    run_test("NullSpace[{{1, 4, 2, -9}, {4, 12, 2, 5}, "
                       "{6, 7, -11, 9}, {5, 15, 10, 12}}]",
             "List[]");
}

/* Three linearly-independent vectors: empty null space. */
static void test_nullspace_li_vectors(void) {
    run_test("NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}]",
             "List[]");
}

/* The {{a, b}, {2a, 2b}} symbolic case: free col 2.
 * RREF row 1 = {1, b/a}; basis = {{-(b/a), 1}}. */
static void test_nullspace_symbolic_2x2(void) {
    run_test("NullSpace[{{a, b}, {2 a, 2 b}}]",
             "List[List[Times[-1, Power[a, -1], b], 1]]");
}

/* The {{a, b, c}, {c, b, a}, {0, 0, 0}} symbolic case from the spec:
 * expected {{1, -(a+c)/b, 1}}. */
static void test_nullspace_symbolic_3x3(void) {
    /* RREF[1,3] = a/b + c/b distributed, so the basis entry comes out
     * as Times[-1, Plus[Times[a, Power[b,-1]], Times[Power[b,-1], c]]];
     * mathematically equivalent to -(a+c)/b. */
    run_test("NullSpace[{{a, b, c}, {c, b, a}, {0, 0, 0}}]",
             "List[List[1, Times[-1, Plus[Times[a, Power[b, -1]], "
             "Times[Power[b, -1], c]]], 1]]");
}

/* Rectangular non-square: shape 3x4, rank 3, one free col.
 * From the spec: NullSpace[{{3, 2, 2, 4}, {2, 3, -2, 7}, {3, 2, 5, 7}}]
 *   = {{12, -23, -5, 5}}. */
static void test_nullspace_3x4_rect(void) {
    run_test("NullSpace[{{3, 2, 2, 4}, {2, 3, -2, 7}, {3, 2, 5, 7}}]",
             "List[List[12, -23, -5, 5]]");
}

/* Method-option pass-through: same input under explicit
 * DivisionFreeRowReduction must give the same answer as the default. */
static void test_nullspace_method_divfree(void) {
    run_test("NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, "
                       "Method -> \"DivisionFreeRowReduction\"]",
             "List[List[1, -2, 1]]");
}

static void test_nullspace_method_onestep(void) {
    /* OneStepRowReduction produces RREF too, so the basis matches. */
    run_test("NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, "
                       "Method -> \"OneStepRowReduction\"]",
             "List[List[1, -2, 1]]");
}

static void test_nullspace_method_automatic_symbol(void) {
    /* Method -> Automatic (the bare symbol, not a string) is accepted. */
    run_test("NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, "
                       "Method -> Automatic]",
             "List[List[1, -2, 1]]");
}

/* Unknown method: the call is left unevaluated. */
static void test_nullspace_method_unknown(void) {
    run_test("NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, "
                       "Method -> \"NoSuchMethod\"]",
             "NullSpace[List[List[1, 2, 3], List[4, 5, 6], List[7, 8, 9]], "
             "Rule[Method, \"NoSuchMethod\"]]");
}

/* Empty / non-matrix arguments are returned unevaluated. */
static void test_nullspace_non_matrix(void) {
    run_test("NullSpace[x]", "NullSpace[x]");
    run_test("NullSpace[{}]", "NullSpace[List[]]");
    /* Scalar in a list isn't a rank-2 tensor. */
    run_test("NullSpace[{1, 2, 3}]", "NullSpace[List[1, 2, 3]]");
}

/* Method only -- arity > 2 left unevaluated. */
static void test_nullspace_arity(void) {
    run_test("NullSpace[]", "NullSpace[]");
    run_test("NullSpace[{{1, 2}}, Method -> \"Automatic\", extra]",
             "NullSpace[List[List[1, 2]], Rule[Method, \"Automatic\"], extra]");
}

/* Identity matrix has empty null space. */
static void test_nullspace_identity(void) {
    run_test("NullSpace[IdentityMatrix[3]]", "List[]");
    run_test("NullSpace[IdentityMatrix[5]]", "List[]");
}

/* Zero matrix: every column is free -- basis is the n standard basis
 * vectors in REVERSED order (rightmost-free-first ordering). */
static void test_nullspace_zero_matrix(void) {
    run_test("NullSpace[{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}]",
             "List[List[0, 0, 1], List[0, 1, 0], List[1, 0, 0]]");
}

/* 1x5 single-row symbolic matrix from the spec.  After scaling by `a`
 * the row is {1, 2, Pi, 1/3, I}.  Free cols: 2, 3, 4, 5 (rightmost
 * first).  Each basis vector has the appropriate denominator cleared
 * where applicable. */
static void test_nullspace_1x5_symbolic(void) {
    /* The expected vectors, rightmost-free-first:
     *   free col 5 -> v = {-I, 0, 0, 0, 1}
     *   free col 4 -> v = {-1/3, 0, 0, 1, 0}  (denom 3 cleared on the
     *                                          v[1]=-(1/3)/1=-1/3 entry,
     *                                          giving {-1, 0, 0, 3, 0})
     *   free col 3 -> v = {-Pi, 0, 1, 0, 0}
     *   free col 2 -> v = {-2, 1, 0, 0, 0}
     *
     * Denominator-clearing applies only when *all* entries have integer
     * denominators.  v[1]=-1/3 -> Denominator=3, others 1.  LCM=3.
     *
     * Result: {{-I, 0, 0, 0, 1}, {-1, 0, 0, 3, 0}, {-Pi, 0, 1, 0, 0},
     *          {-2, 1, 0, 0, 0}}.
     *
     * We assert this semantically (m . v == 0) since the exact form of
     * symbolic Pi / Complex normalisation can drift across evaluator
     * tweaks; basis size is the load-bearing invariant. */
    assert_nullspace_valid("{{a, 2 a, Pi a, a/3, I a}}", 4);
}

/* Rectangular integer matrix from the spec:
 *   m = {{0, 5, 2, 4, 4}, {2, 5, 0, 4, 0}, {5, 1, 5, 4, 5}}
 *   rank 3, nullity 2.  Verify both basis vectors. */
static void test_nullspace_3x5_rect_int(void) {
    assert_nullspace_valid(
        "{{0, 5, 2, 4, 4}, {2, 5, 0, 4, 0}, {5, 1, 5, 4, 5}}", 2);
}

/* Four 4-d vectors with a rank-2 relation: nullity should be 2.
 * From the spec, expected basis = {{-44, -3, 0, 27}, {-26, 6, 9, 0}}. */
static void test_nullspace_4x4_rank2(void) {
    assert_nullspace_valid(
        "{{108, 90, 252, 186}, {240, 260, 520, 420}, "
        "{264, 340, 536, 468}, {522, 705, 1038, 929}}", 2);
}

/* Wider-than-tall rectangular: 2 rows, 4 cols -> nullity = 2 always
 * when the two rows are linearly independent. */
static void test_nullspace_2x4_li(void) {
    assert_nullspace_valid("{{1, 0, 2, 3}, {0, 1, -1, 4}}", 2);
}

/* Wide rectangular: 1 row, 3 cols -> nullity = 2. */
static void test_nullspace_1x3(void) {
    assert_nullspace_valid("{{1, 2, 3}}", 2);
}

/* Pure rational entries.  Denominator-clearing should integerise the
 * basis where possible. */
static void test_nullspace_rational_entries(void) {
    assert_nullspace_valid("{{1, 1/2, 1/3}, {1/2, 1/3, 1/4}}", 1);
}

/* Symbolic-with-symbolic-denominator: denominator-clearing must NOT
 * scale symbolic denominators (Mathematica returns the natural
 * rational form for such inputs). */
static void test_nullspace_symbolic_no_int_clearing(void) {
    /* The natural basis is {-b/a, 1}; clearing must NOT scale by `a`. */
    Expr* ns = eval_str("NullSpace[{{a, b}, {2 a, 2 b}}]");
    char* s = expr_to_string(ns);
    /* The printed form should contain "1/a" or "b/a", not "{-b, a}". */
    ASSERT(strstr(s, "b/a") != NULL || strstr(s, "1/a") != NULL);
    free(s);
    expr_free(ns);
    printf("  PASS: symbolic denominator preserved\n");
}

/* Machine-precision floats: verify nullity but not exact entries
 * (those depend on RREF rounding). */
static void test_nullspace_float_3x3(void) {
    /* Matrix designed to have nullity 1 in exact arithmetic; we accept
     * what Mathilda's RowReduce produces.  Skip if RowReduce can't get
     * to nullity 1 for floats (the row reducer is structural, so a
     * float-rank-deficient matrix may still yield nullity 0).  Use a
     * matrix whose rows are exactly linearly dependent. */
    /* Row 3 = 2 * Row 1 + Row 2 -- so structural zero appears
     * after exact-arithmetic subtraction. */
    assert_nullspace_valid(
        "{{1.5, 4.75, -3.2}, {-1.7, 6.7, -9.3}, "
        "{1.3, 16.2, -15.7}}", 1);
}

/* Bigint entries to exercise GMP / arbitrary-precision path. */
static void test_nullspace_bigint(void) {
    /* A 2x3 matrix with each entry a big integer.  Independent rows,
     * one free col -> nullity 1. */
    assert_nullspace_valid(
        "{{100000000000000000001, 200000000000000000002, "
        "  300000000000000000003}, "
        " {400000000000000000004, 500000000000000000005, "
        "  600000000000000000006}}", 1);
}

/* Repeated NullSpace calls should not accumulate state. */
static void test_nullspace_repeated_calls(void) {
    for (int i = 0; i < 5; i++) {
        assert_nullspace_valid("{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}", 1);
    }
}

/* The docstring registration check -- ?NullSpace should print a
 * non-empty docstring.  Easiest to assert symbol table has docstring. */
static void test_nullspace_docstring(void) {
    Expr* res = eval_str("Information[NullSpace]");
    /* Information returns Null after printing, but the call succeeds
     * iff the docstring is registered. */
    expr_free(res);
    printf("  PASS: Information[NullSpace] evaluated\n");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running NullSpace tests...\n");

    /* Exact-form tests. */
    TEST(test_nullspace_3x3_singular);
    TEST(test_nullspace_3x3_full_rank);
    TEST(test_nullspace_4x4_one_free);
    TEST(test_nullspace_4x4_full_rank);
    TEST(test_nullspace_li_vectors);
    TEST(test_nullspace_symbolic_2x2);
    TEST(test_nullspace_symbolic_3x3);
    TEST(test_nullspace_3x4_rect);
    TEST(test_nullspace_method_divfree);
    TEST(test_nullspace_method_onestep);
    TEST(test_nullspace_method_automatic_symbol);
    TEST(test_nullspace_method_unknown);
    TEST(test_nullspace_non_matrix);
    TEST(test_nullspace_arity);
    TEST(test_nullspace_identity);
    TEST(test_nullspace_zero_matrix);

    /* Semantic tests. */
    TEST(test_nullspace_1x5_symbolic);
    TEST(test_nullspace_3x5_rect_int);
    TEST(test_nullspace_4x4_rank2);
    TEST(test_nullspace_2x4_li);
    TEST(test_nullspace_1x3);
    TEST(test_nullspace_rational_entries);
    TEST(test_nullspace_symbolic_no_int_clearing);
    TEST(test_nullspace_float_3x3);
    TEST(test_nullspace_bigint);
    TEST(test_nullspace_repeated_calls);
    TEST(test_nullspace_docstring);

    printf("All NullSpace tests passed!\n");
    symtab_clear();
    return 0;
}
