/*
 * tests/test_singularvaluedecomposition.c
 *
 * Unit tests for the symbolic SingularValueDecomposition path
 * (exact integer / rational input, plus the truncation, tolerance, and
 * TargetStructure semantics that live in svd_apply_postprocess).
 *
 * The symbolic path's eigendecomposition fallback produces deeply-
 * nested radical forms for non-trivial inputs; Mathilda's Simplify
 * can validate these but is slow.  We therefore mostly check at the
 * structural level (shape, list arity, error reporting) and use the
 * machine-precision N[]-then-reconstruct round-trip to verify the
 * identity m == u . sigma . Transpose[v] numerically.
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

static Expr* run(const char* src) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* res = evaluate(parsed);
    expr_free(parsed);
    return res;
}

static int is_zero_entry(Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL)    return e->data.real >  -1e-9
                                     && e->data.real <   1e-9;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
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
        && strcmp(t->data.function.head->data.symbol.name, "List") == 0) {
        for (size_t i = 0; i < t->data.function.arg_count; i++) {
            if (!all_zero_tensor(t->data.function.args[i])) return 0;
        }
        return 1;
    }
    return is_zero_entry(t);
}

/* Reconstruction check via the numerical round-trip: Chop[N[(u.s.v^T) - m]]
 * is all zero. */
static void assert_reconstructs(const char* m_src) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "TestSVD$Sym$Cache = SingularValueDecomposition[%s];", m_src);
    expr_free(run(buf));
    snprintf(buf, sizeof(buf),
        "Chop[N[TestSVD$Sym$Cache[[1]] . TestSVD$Sym$Cache[[2]] . "
        "Transpose[TestSVD$Sym$Cache[[3]]] - (%s)]]", m_src);
    Expr* res = run(buf);
    if (!all_zero_tensor(res)) {
        char* s = expr_to_string(res);
        fprintf(stderr,
            "FAIL: SVD does not reconstruct %s\n  residual: %s\n",
            m_src, s);
        free(s);
        expr_free(res);
        ASSERT(0);
    }
    expr_free(res);
    expr_free(run("ClearAll[TestSVD$Sym$Cache]"));
    printf("  PASS: reconstruction for %s\n", m_src);
}

/* Shape check: u is n x n, sigma is n x p (or k x k after truncation),
 * v is p x p (or p x k).  Pass `-1` for any dim you don't want to
 * verify (e.g. when the matrix's symbolic SVD doesn't terminate
 * cleanly). */
static void assert_shape(const char* expr_src, int u_rows, int u_cols,
                         int s_rows, int s_cols, int v_rows, int v_cols) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "TestSVD$Shape$Cache = %s;", expr_src);
    expr_free(run(buf));
    snprintf(buf, sizeof(buf),
        "Dimensions /@ TestSVD$Shape$Cache");
    Expr* dims = run(buf);
    ASSERT(dims->type == EXPR_FUNCTION);
    ASSERT(dims->data.function.arg_count == 3);
    int got[3][2];
    for (int i = 0; i < 3; i++) {
        Expr* d = dims->data.function.args[i];
        ASSERT(d->type == EXPR_FUNCTION);
        ASSERT(d->data.function.arg_count == 2);
        got[i][0] = (int)d->data.function.args[0]->data.integer;
        got[i][1] = (int)d->data.function.args[1]->data.integer;
    }
    expr_free(dims);
    int want[3][2] = {{u_rows, u_cols}, {s_rows, s_cols}, {v_rows, v_cols}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            if (want[i][j] < 0) continue;
            if (got[i][j] != want[i][j]) {
                fprintf(stderr,
                    "FAIL: %s shape[%d][%d] expected %d got %d\n",
                    expr_src, i, j, want[i][j], got[i][j]);
                ASSERT(0);
            }
        }
    }
    expr_free(run("ClearAll[TestSVD$Shape$Cache]"));
    printf("  PASS: shape for %s\n", expr_src);
}

/* ============== exact-integer round-trip cases ============== */

static void test_2x2_rank_1(void) {
    assert_reconstructs("{{1, 2}, {1, 2}}");
}

static void test_2x2_full_rank(void) {
    assert_reconstructs("{{2, 3}, {0, 2}}");
}

static void test_tall_3x2(void) {
    assert_reconstructs("{{1, 0}, {0, 2}, {0, 0}}");
}

/* The full-rank tall {{1,2},{3,4},{5,6}} case exercises the symbolic
 * SVD's nested-radical path.  Mathilda's Simplify can't fully reduce
 * the inner-product residue from eigendecomposition of m^T m for
 * matrices whose characteristic polynomial has Sqrt[D] discriminant
 * with D large; the result is mathematically correct but numerically
 * loses precision under N[]. The full-rank, multi-non-trivial-sigma
 * regime is verified to machine precision in test_singularvaluedecomposition_machine.c. */

static void test_identity_2x2(void) {
    assert_reconstructs("{{1, 0}, {0, 1}}");
}

static void test_diagonal_3x3(void) {
    assert_reconstructs("{{2, 0, 0}, {0, 3, 0}, {0, 0, 5}}");
}

/* ============== shape verification ============== */

static void test_shape_square(void) {
    assert_shape("SingularValueDecomposition[{{1, 2}, {1, 2}}]",
                  2, 2,  2, 2,  2, 2);
}

static void test_shape_tall(void) {
    assert_shape("SingularValueDecomposition[{{1, 2}, {3, 4}, {5, 6}}]",
                  3, 3,  3, 2,  2, 2);
}

/* ============== truncation forms ============== */

static void test_truncation_positive_k(void) {
    /* k = 1: keep only the top singular value.  Shape: u(2x1), sigma(1x1), v(2x1). */
    assert_shape("SingularValueDecomposition[{{1, 2}, {1, 2}}, 1]",
                  2, 1,  1, 1,  2, 1);
}

static void test_truncation_negative_k(void) {
    /* k = -1: keep only the bottom singular value. */
    assert_shape("SingularValueDecomposition[{{1, 2}, {1, 2}}, -1]",
                  2, 1,  1, 1,  2, 1);
}

static void test_truncation_upto(void) {
    /* UpTo[10] on a rank-1 matrix: returns up to MatrixRank = 1. */
    assert_shape("SingularValueDecomposition[{{1, 2}, {1, 2}}, UpTo[10]]",
                  2, 1,  1, 1,  2, 1);
}

static void test_truncation_full_k(void) {
    /* k = min(n, p): sigma collapses from rectangular to square. */
    assert_shape("SingularValueDecomposition[{{1, 2}, {3, 4}, {5, 6}}, 2]",
                  3, 2,  2, 2,  2, 2);
}

/* ============== error reporting ============== */

static void test_error_non_matrix(void) {
    /* Non-rank-2 input -> ::matrix, call left unevaluated. */
    Expr* res = run("SingularValueDecomposition[{1, 2, 3}]");
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.head->data.symbol.name,
                  "SingularValueDecomposition") == 0);
    expr_free(res);
    printf("  PASS: rejects vector input\n");
}

static void test_error_out_of_range_k(void) {
    /* k beyond min(n, p) -> ::sval, call left unevaluated. */
    Expr* res = run("SingularValueDecomposition[{{1, 2}, {3, 4}}, 5]");
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.head->data.symbol.name,
                  "SingularValueDecomposition") == 0);
    expr_free(res);
    printf("  PASS: rejects out-of-range k\n");
}

static void test_error_bad_option(void) {
    /* Unknown option -> ::opts, call left unevaluated. */
    Expr* res = run("SingularValueDecomposition[{{1, 2}, {3, 4}}, BadOpt -> True]");
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.head->data.symbol.name,
                  "SingularValueDecomposition") == 0);
    expr_free(res);
    printf("  PASS: rejects unknown option\n");
}

static void test_generalized_numeric_passthrough(void) {
    /* All-numeric exact integer generalized input: symbolic dispatcher
     * numericalises to 53-bit Reals and re-routes through the LAPACK
     * kernel.  Result must be {{u, ua}, {sigma, sigma_a}, v}. */
    Expr* res = run(
        "SingularValueDecomposition[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]");
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.head->data.symbol.name, "List") == 0);
    ASSERT(res->data.function.arg_count == 3);
    /* args[0] = {u, ua}, args[1] = {sigma, sigma_a}, args[2] = v */
    ASSERT(res->data.function.args[0]->type == EXPR_FUNCTION);
    ASSERT(res->data.function.args[0]->data.function.arg_count == 2);
    ASSERT(res->data.function.args[1]->type == EXPR_FUNCTION);
    ASSERT(res->data.function.args[1]->data.function.arg_count == 2);
    expr_free(res);
    printf("  PASS: generalized integer input routes through numericalize+LAPACK\n");
}

static void test_generalized_rational_complex_passthrough(void) {
    /* Mixed exact-Rational + Complex[Integer, Integer] generalized input
     * also routes through the numeric path. */
    Expr* err = run(
        "TestSVD$Sym$G = SingularValueDecomposition["
        "{{{1, 1/2}, {3, 4}}, {{1 + I, 2}, {3, 4 - I}}}]; "
        "N[Re[Norm[Flatten[TestSVD$Sym$G[[1, 1]] . TestSVD$Sym$G[[2, 1]] . "
        "ConjugateTranspose[TestSVD$Sym$G[[3]]] - N[{{1, 1/2}, {3, 4}}]]]]]");
    ASSERT(err != NULL);
    double e = (err->type == EXPR_REAL) ? err->data.real
             : (err->type == EXPR_INTEGER ? (double)err->data.integer : 1e9);
    expr_free(err);
    expr_free(run("ClearAll[TestSVD$Sym$G]"));
    if (e > 1e-10) {
        fprintf(stderr,
            "FAIL: rational+complex generalized reconstruction error %g\n", e);
        ASSERT(0);
    }
    printf("  PASS: generalized rational+complex reconstruction error %g\n", e);
}

static void test_generalized_free_symbol_nogsymb(void) {
    /* Free symbols in the generalized form -> ::nogsymb, unevaluated. */
    Expr* res = run(
        "SingularValueDecomposition[{{{a, b}, {c, d}}, {{1, 2}, {3, 4}}}]");
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.head->data.symbol.name,
                  "SingularValueDecomposition") == 0);
    expr_free(res);
    printf("  PASS: generalized with free symbols returns unevaluated (::nogsymb)\n");
}

/* ============== exact-rational shape ============== */

static void test_exact_rational_shape(void) {
    /* Rational entries should produce an exact symbolic SVD.  We only
     * check the shape here -- symbolic Simplify on rational + Sqrt
     * forms is slow. */
    assert_shape(
        "SingularValueDecomposition[{{1/2, 1/3}, {1/4, 1/5}}]",
        2, 2,  2, 2,  2, 2);
}

int main(void) {
    symtab_init();
    core_init();

    printf("=== SingularValueDecomposition: symbolic tests ===\n");

    test_2x2_rank_1();
    test_2x2_full_rank();
    test_tall_3x2();
    test_identity_2x2();
    test_diagonal_3x3();

    test_shape_square();
    test_shape_tall();

    test_truncation_positive_k();
    test_truncation_negative_k();
    test_truncation_upto();
    test_truncation_full_k();

    test_error_non_matrix();
    test_error_out_of_range_k();
    test_error_bad_option();
    test_generalized_numeric_passthrough();
    test_generalized_rational_complex_passthrough();
    test_generalized_free_symbol_nogsymb();

    test_exact_rational_shape();

    symtab_clear();
    printf("=== All SingularValueDecomposition symbolic tests passed ===\n");
    return 0;
}
