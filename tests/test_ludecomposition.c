/*
 * tests/test_ludecomposition.c
 *
 * Unit tests for LUDecomposition[m].
 *
 * Verification strategies (mirroring tests/test_qrdecomposition.c):
 *
 *   1. Public identity:
 *          m[[p]] == l . u
 *      where l = LowerTriangularize[lu, -1] + IdentityMatrix[n] and
 *      u = UpperTriangularize[lu].  Checked by evaluating the
 *      difference through Together + Chop and asserting every leaf is
 *      a literal-zero form (Integer 0, Rational[0, _], Real 0.0,
 *      Complex[0, 0]).
 *
 *   2. Shape: the returned list has three elements; lu is n x n; p is
 *      a 1-indexed length-n permutation; c is either an Integer 0 or
 *      a Real.
 *
 *   3. p is a valid permutation of 1..n.
 *
 *   4. For the exact / symbolic cases we additionally print-compare
 *      against known canonical strings (regression).
 *
 * Each test frees every parsed / evaluated Expr.  The binary is meant
 * to run cleanly under `valgrind --leak-check=full`.
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

/* Parse + evaluate.  Caller owns the returned Expr. */
static Expr* run(const char* src) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* res = evaluate(parsed);
    expr_free(parsed);
    return res;
}

/* Print-compare regression check. */
static void run_test(const char* input, const char* expected) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    Expr* res = evaluate(parsed);
    char* res_str = expr_to_string(res);
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

/* Treat e as "approximately zero" -- mirrors test_qrdecomposition.c. */
static int is_zero_entry(Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL)    return e->data.real == 0.0
                                        || (e->data.real >  -1e-9
                                         && e->data.real <   1e-9);
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (strcmp(h, "Complex") == 0 && e->data.function.arg_count == 2) {
            return is_zero_entry(e->data.function.args[0])
                && is_zero_entry(e->data.function.args[1]);
        }
        if (strcmp(h, "Rational") == 0 && e->data.function.arg_count == 2) {
            return is_zero_entry(e->data.function.args[0]);
        }
    }
    return 0;
}

static int all_zero_tensor(Expr* t) {
    if (!t) return 0;
    if (t->type == EXPR_FUNCTION
        && t->data.function.head->type == EXPR_SYMBOL
        && strcmp(t->data.function.head->data.symbol, "List") == 0) {
        for (size_t i = 0; i < t->data.function.arg_count; i++) {
            if (!all_zero_tensor(t->data.function.args[i])) return 0;
        }
        return 1;
    }
    return is_zero_entry(t);
}

/* Assert m[[p]] == l . u where {lu, p, _} = LUDecomposition[m_src].
 *
 * `LowerTriangularize` / `UpperTriangularize` aren't builtin in
 * Mathilda yet, so we construct L (unit lower) and U (upper) via
 * Table / If on the indices. */
static void assert_lu_identity(const char* m_src) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
             "Block[{m, lu, p, l, u, n, ii, jj}, "
             "  m = %s; "
             "  {lu, p} = LUDecomposition[m][[1;;2]]; "
             "  n = Length[m]; "
             "  l = Table[If[ii > jj, lu[[ii, jj]], "
             "                  If[ii == jj, 1, 0]], "
             "             {ii, n}, {jj, n}]; "
             "  u = Table[If[ii <= jj, lu[[ii, jj]], 0], "
             "             {ii, n}, {jj, n}]; "
             "  Together[m[[p]] - l . u]"
             "]",
             m_src);
    Expr* res = run(buf);
    if (!all_zero_tensor(res)) {
        char* s = expr_to_string(res);
        fprintf(stderr,
                "FAIL: LU identity failed for %s\n  diff: %s\n",
                m_src, s);
        free(s);
        expr_free(res);
        ASSERT(0);
    }
    printf("  PASS: m[[p]] == l.u  for %s\n", m_src);
    expr_free(res);
}

/* Assert p is a permutation of {1..n}. */
static void assert_perm_valid(const char* m_src, int n) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
             "Sort[LUDecomposition[%s][[2]]] == Range[%d]",
             m_src, n);
    Expr* res = run(buf);
    ASSERT(res->type == EXPR_SYMBOL && strcmp(res->data.symbol, "True") == 0);
    expr_free(res);
    printf("  PASS: p is a permutation of 1..%d  for %s\n", n, m_src);
}

/* Assert the shape: result is {lu, p, c} with lu n x n and p length n. */
static void assert_shape(const char* m_src, int n) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
             "{Length[LUDecomposition[%s]], "
             " Dimensions[LUDecomposition[%s][[1]]], "
             " Length[LUDecomposition[%s][[2]]]}",
             m_src, m_src, m_src);
    Expr* res = run(buf);
    char* s = expr_to_string(res);
    char expect[128];
    snprintf(expect, sizeof(expect), "{3, {%d, %d}, %d}", n, n, n);
    if (strcmp(s, expect) != 0) {
        fprintf(stderr,
                "FAIL: bad shape for %s\n  expected: %s\n  got:      %s\n",
                m_src, expect, s);
        free(s);
        expr_free(res);
        ASSERT(0);
    }
    printf("  PASS: shape  %s -> %s\n", m_src, s);
    free(s);
    expr_free(res);
}

/* The full triple of checks. */
static void assert_lu_valid(const char* m_src, int n) {
    assert_shape(m_src, n);
    assert_perm_valid(m_src, n);
    assert_lu_identity(m_src);
}

/* ============================================================
 *   Exact tests
 * ============================================================ */

static void test_lu_3x3_vandermonde(void) {
    /* The Mathematica documentation example. */
    assert_lu_valid("{{1, 1, 1}, {2, 4, 8}, {3, 9, 27}}", 3);
    run_test("LUDecomposition[{{1, 1, 1}, {2, 4, 8}, {3, 9, 27}}]",
             "{{{1, 1, 1}, {2, 2, 6}, {3, 3, 6}}, {1, 2, 3}, 0}");
}

static void test_lu_2x2_symbolic(void) {
    assert_lu_valid("{{a, b}, {c, d}}", 2);
    /* The classic textbook closed form.  Mathilda's Together canonicalises
     * `-(b c / a) + d` to `(-b c + a d) / a`, equivalent to the
     * Mathematica documentation example. */
    run_test("LUDecomposition[{{a, b}, {c, d}}]",
             "{{{a, b}, {c/a, (-b c + a d)/a}}, {1, 2}, 0}");
}

static void test_lu_2x2_integer(void) {
    assert_lu_valid("{{4, 3}, {6, 3}}", 2);
}

static void test_lu_3x3_rational(void) {
    assert_lu_valid("{{1, 1/2, 1/3}, {1/2, 1/3, 1/4}, {1/3, 1/4, 1/5}}", 3);
}

static void test_lu_3x3_general_symbolic(void) {
    assert_lu_valid("{{a, b, c}, {d, e, f}, {g, h, i}}", 3);
}

static void test_lu_4x4_integer(void) {
    assert_lu_valid("{{2, -1, 0, 1}, {3, 0, 4, -2}, "
                    " {1, 5, -3, 2}, {4, 2, 1, 0}}", 4);
}

static void test_lu_identity(void) {
    assert_lu_valid("IdentityMatrix[5]", 5);
    run_test("LUDecomposition[IdentityMatrix[3]][[1;;2]]",
             "{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {1, 2, 3}}");
}

static void test_lu_complex(void) {
    /* Complex 2x2.  Verify identity through the round-trip; exact
     * print regression is too brittle on complex output. */
    assert_lu_valid("{{2 + 4 I, 9 + 9 I}, {2 + 9 I, 1 + 3 I}}", 2);
}

static void test_lu_singular_zero_pivot(void) {
    /* First column has a zero on the diagonal but a non-zero below;
     * the symbolic kernel should pivot rather than emit ::sing. */
    assert_lu_valid("{{0, 1, 2}, {1, 0, 3}, {4, 5, 0}}", 3);
}

static void test_lu_truly_singular(void) {
    /* Third row is the sum of rows 1 and 2; matrix is rank 2.  We
     * just check the call returns a well-shaped result -- the
     * singular warning is informational. */
    Expr* res = run("LUDecomposition[{{1, 2, 3}, {4, 5, 6}, {5, 7, 9}}]");
    ASSERT(res->type == EXPR_FUNCTION
        && res->data.function.head->type == EXPR_SYMBOL
        && strcmp(res->data.function.head->data.symbol, "List") == 0
        && res->data.function.arg_count == 3);
    expr_free(res);
    printf("  PASS: singular matrix returns a {lu, p, c} triple\n");
}

static void test_lu_p_indexing_works(void) {
    /* Cross-check m[[p]] indexing produces the expected permuted rows.
     * We rebuild L (unit lower) and U (upper) explicitly via Table so
     * the test doesn't depend on builtins we haven't added yet. */
    run_test("Block[{m, lu, p, l, u, n, ii, jj}, "
             "  m = {{2, -1}, {3, 4}}; "
             "  {lu, p} = LUDecomposition[m][[1;;2]]; "
             "  n = Length[m]; "
             "  l = Table[If[ii > jj, lu[[ii, jj]], "
             "                  If[ii == jj, 1, 0]], "
             "             {ii, n}, {jj, n}]; "
             "  u = Table[If[ii <= jj, lu[[ii, jj]], 0], "
             "             {ii, n}, {jj, n}]; "
             "  m[[p]] == l . u"
             "]",
             "True");
}

static void test_lu_bad_inputs(void) {
    /* Non-square - call should be left unevaluated (returns the
     * literal LUDecomposition[...] expression). */
    Expr* res = run("LUDecomposition[{{1, 2, 3}, {4, 5, 6}}]");
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.head->data.symbol,
                  "LUDecomposition") == 0);
    expr_free(res);

    /* Empty matrix. */
    res = run("LUDecomposition[{}]");
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(strcmp(res->data.function.head->data.symbol,
                  "LUDecomposition") == 0);
    expr_free(res);
    printf("  PASS: bad inputs left unevaluated\n");
}

static void test_lu_attributes(void) {
    /* Protected attribute must be set. */
    Expr* res = run("MemberQ[Attributes[LUDecomposition], Protected]");
    ASSERT(res->type == EXPR_SYMBOL
        && strcmp(res->data.symbol, "True") == 0);
    expr_free(res);
    printf("  PASS: LUDecomposition has Protected attribute\n");
}

int main(void) {
    symtab_init();
    core_init();

    /* ---- shape / identity ---- */
    TEST(test_lu_3x3_vandermonde);
    TEST(test_lu_2x2_symbolic);
    TEST(test_lu_2x2_integer);
    TEST(test_lu_3x3_rational);
    TEST(test_lu_3x3_general_symbolic);
    TEST(test_lu_4x4_integer);
    TEST(test_lu_identity);

    /* ---- complex ---- */
    TEST(test_lu_complex);

    /* ---- singular handling ---- */
    TEST(test_lu_singular_zero_pivot);
    TEST(test_lu_truly_singular);

    /* ---- regression ---- */
    TEST(test_lu_p_indexing_works);

    /* ---- error paths ---- */
    TEST(test_lu_bad_inputs);
    TEST(test_lu_attributes);

    printf("All LUDecomposition tests passed.\n");
    return 0;
}
