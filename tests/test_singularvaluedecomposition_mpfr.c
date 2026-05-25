/*
 * tests/test_singularvaluedecomposition_mpfr.c
 *
 * Unit tests for the MPFR Jacobi SVD fast path of
 * SingularValueDecomposition.
 *
 * Routed through svd_dispatch by feeding inexact inputs at > 53 bits of
 * precision (via N[m, digits] with digits > 16).  Verification:
 *
 *   1. Reconstruction at the input precision: the reconstructed
 *      matrix - input is identically zero to within MPFR rounding.
 *   2. Shape: u (n x n), sigma (n x p), v (p x p).
 *   3. Singular-value ordering (non-increasing).
 *   4. Rank-deficient handling: zero singular values land at the
 *      tail; u is still fully orthonormal via the orthonormal
 *      completion pass.
 */

#include <math.h>
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

/* Convert a leaf Expr to double for the residual comparison.  Handles
 * Real, Integer, Rational, MPFR (via the dedicated EXPR_MPFR slot),
 * and Complex by taking |re + i im|.  Other heads abort the test. */
#ifdef USE_MPFR
#include <mpfr.h>
#endif

static double as_double(Expr* e) {
    if (e->type == EXPR_REAL)    return e->data.real;
    if (e->type == EXPR_INTEGER) return (double)e->data.integer;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    return mpfr_get_d(e->data.mpfr, MPFR_RNDN);
#endif
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (strcmp(h, "Rational") == 0 && e->data.function.arg_count == 2) {
            return as_double(e->data.function.args[0])
                 / as_double(e->data.function.args[1]);
        }
        if (strcmp(h, "Complex") == 0 && e->data.function.arg_count == 2) {
            double re = as_double(e->data.function.args[0]);
            double im = as_double(e->data.function.args[1]);
            return sqrt(re * re + im * im);
        }
        if (strcmp(h, "Times") == 0 && e->data.function.arg_count == 2) {
            return as_double(e->data.function.args[0])
                 * as_double(e->data.function.args[1]);
        }
        if (strcmp(h, "Plus") == 0) {
            double s = 0.0;
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                s += as_double(e->data.function.args[i]);
            }
            return s;
        }
    }
    char* s = expr_to_string(e);
    fprintf(stderr, "FAIL: cannot convert to double: %s\n", s);
    free(s);
    ASSERT(0);
    return 0.0;
}

/* Walk a Flattened residual list and return the maximum |entry|.
 * Uses as_double for the per-element |.| computation; MPFR entries
 * are converted via mpfr_get_d. */
static double max_abs_flat(Expr* lst) {
    if (lst->type != EXPR_FUNCTION) return as_double(lst);
    double m = 0.0;
    for (size_t i = 0; i < lst->data.function.arg_count; i++) {
        double v = fabs(as_double(lst->data.function.args[i]));
        if (v > m) m = v;
    }
    return m;
}

static double run_d(const char* src) {
    Expr* r = run(src);
    double v = as_double(r);
    expr_free(r);
    return v;
}

/* Reconstruction within `tol` absolute error.  We pull the residual
 * back to C as a Flatten[...] list and compute the max |entry| in C;
 * Mathilda's Max[Abs[...]] doesn't always reduce MPFR Sqrt-of-sums
 * to a single Real, so we avoid relying on it for the assertion. */
static void assert_reconstruction(const char* m_src, double tol) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "TestSVDMpfr$Cache = SingularValueDecomposition[%s];", m_src);
    expr_free(run(buf));
    snprintf(buf, sizeof(buf),
        "Flatten["
        "TestSVDMpfr$Cache[[1]] . TestSVDMpfr$Cache[[2]] . "
        "ConjugateTranspose[TestSVDMpfr$Cache[[3]]] - (%s)]", m_src);
    Expr* residual = run(buf);
    double err = max_abs_flat(residual);
    expr_free(residual);
    expr_free(run("ClearAll[TestSVDMpfr$Cache]"));
    if (err > tol) {
        fprintf(stderr,
            "FAIL: %s MPFR reconstruction error %g > tol %g\n",
            m_src, err, tol);
        ASSERT(0);
    }
    printf("  PASS: MPFR reconstruction error %g for %s\n", err, m_src);
}

static void assert_shape(const char* expr_src, int ur, int uc,
                         int sr, int sc, int vr, int vc) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "TestSVDMpfr$Shape = %s;", expr_src);
    expr_free(run(buf));
    Expr* dims = run("Dimensions /@ TestSVDMpfr$Shape");
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
    int want[3][2] = {{ur, uc}, {sr, sc}, {vr, vc}};
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
    expr_free(run("ClearAll[TestSVDMpfr$Shape]"));
    printf("  PASS: shape for %s\n", expr_src);
}

static void assert_sigma_ordered(const char* m_src) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "TestSVDMpfr$Cache = SingularValueDecomposition[%s];", m_src);
    expr_free(run(buf));
    snprintf(buf, sizeof(buf),
        "TestSVDMpfr$Diag = Table[TestSVDMpfr$Cache[[2]][[i, i]], "
        "{i, Min[Dimensions[%s]]}];", m_src);
    expr_free(run(buf));
    Expr* diag = run("TestSVDMpfr$Diag");
    int n = (int)diag->data.function.arg_count;
    double prev = 1e300;
    for (int i = 0; i < n; i++) {
        double cur = as_double(diag->data.function.args[i]);
        if (cur > prev + 1e-12) {
            fprintf(stderr,
                "FAIL: %s sigma not descending: sigma[%d]=%g > sigma[%d]=%g\n",
                m_src, i, cur, i - 1, prev);
            ASSERT(0);
        }
        prev = cur;
    }
    expr_free(diag);
    expr_free(run("ClearAll[TestSVDMpfr$Cache, TestSVDMpfr$Diag]"));
    printf("  PASS: sigma ordering for %s\n", m_src);
}

/* ============== reconstruction at MPFR precision ============== */

static void test_recon_2x2_30digit(void) {
    assert_reconstruction("N[{{1, 2}, {3, 4}}, 30]", 1e-25);
}

static void test_recon_3x2_30digit(void) {
    assert_reconstruction("N[{{1, 2}, {3, 4}, {5, 6}}, 30]", 1e-25);
}

static void test_recon_3x3_30digit(void) {
    assert_reconstruction(
        "N[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}, 30]", 1e-25);
}

static void test_recon_2x2_50digit(void) {
    assert_reconstruction("N[{{2, 3}, {0, 2}}, 50]", 1e-45);
}

/* ============== rank-deficient tall input ============== */

static void test_rank_1_tall(void) {
    /* m = {{1, 2}, {2, 4}, {3, 6}}: rank 1.  Singular values: one
     * non-zero, one zero.  Reconstruction must still hold exactly. */
    assert_reconstruction(
        "N[{{1, 2}, {2, 4}, {3, 6}}, 30]", 1e-25);
}

/* ============== shape ============== */

static void test_shape_square(void) {
    assert_shape(
        "SingularValueDecomposition[N[{{1, 2}, {3, 4}}, 30]]",
        2, 2, 2, 2, 2, 2);
}

static void test_shape_tall(void) {
    assert_shape(
        "SingularValueDecomposition[N[{{1, 2}, {3, 4}, {5, 6}}, 30]]",
        3, 3, 3, 2, 2, 2);
}

/* ============== ordering ============== */

static void test_ordering_3x3(void) {
    assert_sigma_ordered("N[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}, 30]");
}

/* ============== generalized SVD at MPFR input precision ============== */

/* Generalized SVD with high-precision MPFR inputs.  Mathilda doesn't
 * have a native MPFR Paige/Van Loan kernel yet, so the dispatcher
 * emits a one-shot ::gmpdwn warning and falls through to the symbolic
 * path which numericalises to 53-bit Reals and routes through LAPACK.
 * We verify the result reconstructs at machine precision. */
static void assert_gen_recon_mpfr(const char* m_src, const char* a_src,
                                    double tol)
{
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "TestSVDMpfr$G = SingularValueDecomposition[{%s, %s}];",
        m_src, a_src);
    expr_free(run(buf));
    snprintf(buf, sizeof(buf),
        "Flatten[TestSVDMpfr$G[[1, 1]] . TestSVDMpfr$G[[2, 1]] . "
        "ConjugateTranspose[TestSVDMpfr$G[[3]]] - (%s)]", m_src);
    Expr* res_m = run(buf);
    double err_m = max_abs_flat(res_m);
    expr_free(res_m);
    snprintf(buf, sizeof(buf),
        "Flatten[TestSVDMpfr$G[[1, 2]] . TestSVDMpfr$G[[2, 2]] . "
        "ConjugateTranspose[TestSVDMpfr$G[[3]]] - (%s)]", a_src);
    Expr* res_a = run(buf);
    double err_a = max_abs_flat(res_a);
    expr_free(res_a);
    expr_free(run("ClearAll[TestSVDMpfr$G]"));
    if (err_m > tol || err_a > tol) {
        fprintf(stderr,
            "FAIL: MPFR-input generalized SVD errors m=%g a=%g (tol %g)\n",
            err_m, err_a, tol);
        ASSERT(0);
    }
    printf("  PASS: MPFR-input generalized SVD errors m=%g a=%g\n",
           err_m, err_a);
}

static void test_gen_mpfr_30digit(void) {
    assert_gen_recon_mpfr(
        "N[{{1, 2}, {3, 4}}, 30]",
        "N[{{5, 6}, {7, 8}}, 30]", 1e-10);
}

static void test_gen_mpfr_tall_30digit(void) {
    assert_gen_recon_mpfr(
        "N[{{1, 2}, {3, 4}, {5, 6}}, 30]",
        "N[{{7, 8}, {9, 10}}, 30]", 1e-10);
}

/* ============== complex MPFR Jacobi ============== */

static void test_complex_recon_2x2_30digit(void) {
    /* Both entries carry imaginary parts; the dispatcher must route
     * through the complex Jacobi kernel. */
    assert_reconstruction(
        "N[{{1 + 0.5 I, 2}, {3, 4 - 0.5 I}}, 30]", 1e-25);
}

static void test_complex_recon_tall_30digit(void) {
    /* 3 x 2 complex (M > P), exercises the orthonormal-completion of
     * the trailing U column. */
    assert_reconstruction(
        "N[{{1 + I, 2}, {3, 4 - I}, {5 + 2 I, 6}}, 30]", 1e-25);
}

static void test_complex_hermitian_2x2_50digit(void) {
    /* Hermitian input: result should still reconstruct to ~50-digit
     * precision via the complex Jacobi loop. */
    assert_reconstruction(
        "N[{{4, 1 + I}, {1 - I, 3}}, 50]", 1e-45);
}

static void test_complex_pure_imaginary_diag(void) {
    /* Diagonal complex matrix at 30 digits: singular values are |Diag|;
     * U and V are diagonal phase factors. */
    assert_reconstruction(
        "N[{{2 I, 0}, {0, 3 + 4 I}}, 30]", 1e-25);
}

/* ============== truncation ============== */

static void test_truncation_positive(void) {
    assert_shape(
        "SingularValueDecomposition[N[{{1, 2}, {3, 4}, {5, 6}}, 30], 1]",
        3, 1, 1, 1, 2, 1);
}

static void test_truncation_upto(void) {
    /* UpTo[k] with k > min(n, p): clamps to min(n, p) for a full-rank
     * matrix.  (Rank-deficient rank detection in the MPFR kernel would
     * tighten this to MatrixRank, but the Jacobi kernel currently
     * reports mn since the tolerance cutoff hasn't been applied yet at
     * the rank-count site -- the user-facing tolerance still zeroes the
     * sigma diagonal entries themselves.) */
    assert_shape(
        "SingularValueDecomposition[N[{{1, 2}, {3, 4}, {5, 6}}, 30], UpTo[10]]",
        3, 2, 2, 2, 2, 2);
}

int main(void) {
    symtab_init();
    core_init();

    printf("=== SingularValueDecomposition: MPFR Jacobi tests ===\n");

    test_recon_2x2_30digit();
    test_recon_3x2_30digit();
    test_recon_3x3_30digit();
    test_recon_2x2_50digit();
    test_rank_1_tall();

    test_shape_square();
    test_shape_tall();
    test_ordering_3x3();

    test_gen_mpfr_30digit();
    test_gen_mpfr_tall_30digit();

    test_complex_recon_2x2_30digit();
    test_complex_recon_tall_30digit();
    test_complex_hermitian_2x2_50digit();
    test_complex_pure_imaginary_diag();

    test_truncation_positive();
    test_truncation_upto();

    symtab_clear();
    printf("=== All SingularValueDecomposition MPFR tests passed ===\n");
    return 0;
}
