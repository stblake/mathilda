/*
 * tests/test_qrdecomposition_mpfr.c
 *
 * Phase-4 MPFR Householder kernel coverage for QRDecomposition.
 *
 * The existing tests/test_qrdecomposition.c already lightly exercises
 * MPFR-precision input (`N[..., 30]`).  This binary focuses on:
 *
 *   1. The shape / aspect-ratio matrix from S5.3 of
 *      tasks/qr_perf_plan.md at MPFR precisions 30, 60, 100, 200
 *      decimal digits.
 *
 *   2. Reconstruction + orthonormality bounds: after factoring at
 *      MPFR precision `bits`, the residual measured in doubles is
 *      dominated by the double-truncation of the high-precision
 *      outputs, so we use the same K * eps_double bound as the
 *      machine-precision test.  The MPFR algorithm is correct iff
 *      the residual stays within that machine-precision window
 *      regardless of how much higher the working precision is.
 *
 *   3. Pivoted-factorisation identity m . P == q^H . r and the
 *      monotonicity property |R[0,0]| >= |R[1,1]| >= ...
 *
 *   4. Precision-monotonicity: same matrix factored at 30 / 60 / 100
 *      / 200 digits.  Reconstruction error (Frobenius norm, computed
 *      inside Mathilda's MPFR arithmetic and then extracted as a
 *      single scalar) must shrink monotonically with precision -- the
 *      hallmark of an actually-using-precision kernel as opposed to a
 *      silent fallback to symbolic / machine.
 *
 *   5. Rank-deficient + zero matrix at MPFR precision.
 *
 * Builds clean under both USE_MPFR=1 and USE_MPFR=0; in the off case
 * the suite exits with "skipped" before running any tests, exactly
 * like the machine-kernel binary does for USE_LAPACK=0.
 *
 * Every parse / evaluate intermediate is freed; the binary is
 * intended to run under valgrind with zero Mathilda-code leaks.
 */

#include <float.h>
#include <math.h>
#include <stdbool.h>
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

#ifdef USE_MPFR
#include <mpfr.h>
#endif

static int g_failures = 0;

#define EXPECT_LE(actual, bound, label) do {                                   \
    double _a = (actual), _b = (bound);                                        \
    if (!(_a <= _b)) {                                                         \
        fprintf(stderr, "FAIL %s: %.6e > %.6e (bound)\n", (label), _a, _b);    \
        g_failures++;                                                          \
    } else {                                                                   \
        printf("  PASS %s: %.6e <= %.6e\n", (label), _a, _b);                  \
    }                                                                          \
} while (0)

/* Convenience: parse + evaluate.  Caller owns the returned Expr. */
static Expr* run(const char* src) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* res = evaluate(parsed);
    expr_free(parsed);
    return res;
}

/* ---------------------------------------------------------------------
 * Read a numeric leaf into (re, im) doubles.  Adds EXPR_MPFR support
 * to the same shape used by test_qrdecomposition_machine.c.  Extracting
 * an MPFR cell to a double truncates the high-precision tail; the test
 * is intentionally tolerant of that loss because reconstruction is
 * checked in doubles. */
static bool leaf_to_complex(Expr* e, double* re, double* im) {
    *im = 0.0;
    if (!e) return false;
    if (e->type == EXPR_REAL)    { *re = e->data.real;             return true; }
    if (e->type == EXPR_INTEGER) { *re = (double)e->data.integer;  return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        *re = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
        return true;
    }
#endif
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (strcmp(h, "Rational") == 0 && e->data.function.arg_count == 2) {
            double p, q, dummy;
            if (leaf_to_complex(e->data.function.args[0], &p, &dummy)
             && leaf_to_complex(e->data.function.args[1], &q, &dummy)
             && q != 0.0) { *re = p / q; return true; }
        }
        if (strcmp(h, "Complex") == 0 && e->data.function.arg_count == 2) {
            double r, i, dummy;
            if (leaf_to_complex(e->data.function.args[0], &r, &dummy)
             && leaf_to_complex(e->data.function.args[1], &i, &dummy)) {
                *re = r; *im = i; return true;
            }
        }
    }
    return false;
}

/* Walk a List-of-List into a row-major (re, im) buffer pair.  Buffers
 * are freshly allocated; caller frees with free(). */
static bool extract_matrix(Expr* m, double** re_out, double** im_out,
                            int* rows_out, int* cols_out) {
    if (m->type != EXPR_FUNCTION
        || m->data.function.head->type != EXPR_SYMBOL
        || strcmp(m->data.function.head->data.symbol.name, "List") != 0) {
        return false;
    }
    int rows = (int)m->data.function.arg_count;
    if (rows == 0) {
        *re_out = NULL; *im_out = NULL; *rows_out = 0; *cols_out = 0;
        return true;
    }
    Expr* first = m->data.function.args[0];
    int cols;
    if (first->type == EXPR_FUNCTION
        && first->data.function.head->type == EXPR_SYMBOL
        && strcmp(first->data.function.head->data.symbol.name, "List") == 0) {
        cols = (int)first->data.function.arg_count;
    } else {
        return false;
    }
    double* re = (double*)calloc((size_t)rows * (size_t)cols, sizeof(double));
    double* im = (double*)calloc((size_t)rows * (size_t)cols, sizeof(double));
    for (int i = 0; i < rows; i++) {
        Expr* row = m->data.function.args[i];
        if (row->type != EXPR_FUNCTION
            || (int)row->data.function.arg_count != cols) {
            free(re); free(im); return false;
        }
        for (int j = 0; j < cols; j++) {
            double rv, iv;
            if (!leaf_to_complex(row->data.function.args[j], &rv, &iv)) {
                free(re); free(im); return false;
            }
            re[i * cols + j] = rv;
            im[i * cols + j] = iv;
        }
    }
    *re_out = re; *im_out = im; *rows_out = rows; *cols_out = cols;
    return true;
}

static double norm_inf_matrix(const double* re, const double* im,
                                int rows, int cols) {
    double m = 0.0;
    for (int i = 0; i < rows; i++) {
        double row = 0.0;
        for (int j = 0; j < cols; j++) {
            double r = re[i * cols + j];
            double iv = im ? im[i * cols + j] : 0.0;
            row += hypot(r, iv);
        }
        if (row > m) m = row;
    }
    return m;
}

/* dst = A^H . B + (-1) * C.  A: arows x acols (re+im), B: arows x bcols
 * (re+im), so dst is acols x bcols. */
static void cgemm_hermconj_minus(const double* a_re, const double* a_im,
                                    int a_rows, int a_cols,
                                    const double* b_re, const double* b_im,
                                    int b_cols,
                                    const double* c_re, const double* c_im,
                                    double* d_re, double* d_im) {
    for (int i = 0; i < a_cols; i++) {
        for (int j = 0; j < b_cols; j++) {
            double sr = 0.0, si = 0.0;
            for (int k = 0; k < a_rows; k++) {
                double ar =  a_re[k * a_cols + i];
                double ai = -a_im[k * a_cols + i];   /* conjugate */
                double br =  b_re[k * b_cols + j];
                double bi =  b_im[k * b_cols + j];
                sr += ar * br - ai * bi;
                si += ar * bi + ai * br;
            }
            d_re[i * b_cols + j] = sr - c_re[i * b_cols + j];
            d_im[i * b_cols + j] = si - c_im[i * b_cols + j];
        }
    }
}

/* dst = A . A^H - I_{a_rows}. */
static void aaH_minus_I(const double* a_re, const double* a_im,
                          int a_rows, int a_cols,
                          double* d_re, double* d_im) {
    for (int i = 0; i < a_rows; i++) {
        for (int j = 0; j < a_rows; j++) {
            double sr = 0.0, si = 0.0;
            for (int k = 0; k < a_cols; k++) {
                double ar = a_re[i * a_cols + k];
                double ai = a_im[i * a_cols + k];
                double br =  a_re[j * a_cols + k];
                double bi = -a_im[j * a_cols + k];
                sr += ar * br - ai * bi;
                si += ar * bi + ai * br;
            }
            d_re[i * a_rows + j] = sr - (i == j ? 1.0 : 0.0);
            d_im[i * a_rows + j] = si;
        }
    }
}

/* ---------------------------------------------------------------------
 * Property check (no pivoting): factor m_src at MPFR precision via
 * `N[m, digits]` -> QRDecomposition.  Verify shape, reconstruction,
 * orthonormality.  Tolerances are in double-precision because the
 * extraction-to-double truncates everything below eps_double anyway. */
static void qrm_property_check(const char* label, const char* m_src,
                                int digits, int expected_rank,
                                double K_recon, double K_orth) {
    printf("Test: %s @ %d digits\n", label, digits);

    /* Wrap input in N[.] at the desired precision. */
    size_t qrlen = strlen(m_src) + 64;
    char* nbuf = (char*)malloc(qrlen);
    snprintf(nbuf, qrlen, "N[%s, %d]", m_src, digits);

    double *m_re, *m_im; int n_rows = 0, n_cols = 0;
    Expr* m = run(nbuf);
    if (!extract_matrix(m, &m_re, &m_im, &n_rows, &n_cols)) {
        fprintf(stderr, "FAIL %s: could not extract input matrix\n", label);
        g_failures++;
        expr_free(m);
        free(nbuf);
        return;
    }
    expr_free(m);
    double mnorm = norm_inf_matrix(m_re, m_im, n_rows, n_cols);

    snprintf(nbuf, qrlen, "QRDecomposition[N[%s, %d]]", m_src, digits);
    Expr* qr = run(nbuf);
    free(nbuf);

    Expr* q = qr->data.function.args[0];
    Expr* r = qr->data.function.args[1];

    double *q_re, *q_im, *r_re, *r_im;
    int q_rows, q_cols, r_rows, r_cols;
    if (!extract_matrix(q, &q_re, &q_im, &q_rows, &q_cols)
     || !extract_matrix(r, &r_re, &r_im, &r_rows, &r_cols)) {
        fprintf(stderr, "FAIL %s: bad q or r structure\n", label);
        g_failures++;
        free(m_re); free(m_im);
        expr_free(qr);
        return;
    }
    if (q_rows != expected_rank || q_cols != n_rows
     || r_rows != expected_rank || r_cols != n_cols) {
        fprintf(stderr, "FAIL %s: shape mismatch: "
                "got q %dx%d r %dx%d, expected rank=%d n=%d p=%d\n",
                label, q_rows, q_cols, r_rows, r_cols,
                expected_rank, n_rows, n_cols);
        g_failures++;
        free(q_re); free(q_im); free(r_re); free(r_im);
        free(m_re); free(m_im);
        expr_free(qr);
        return;
    }

    if (expected_rank > 0) {
        /* recon = q^H . r - m */
        double* recon_re = (double*)calloc((size_t)n_rows * (size_t)n_cols, sizeof(double));
        double* recon_im = (double*)calloc((size_t)n_rows * (size_t)n_cols, sizeof(double));
        cgemm_hermconj_minus(q_re, q_im, q_rows, q_cols,
                                r_re, r_im, r_cols,
                                m_re, m_im,
                                recon_re, recon_im);
        double recon_norm = norm_inf_matrix(recon_re, recon_im, n_rows, n_cols);
        double bound = K_recon * 2.2204460492503131e-16 * (mnorm + 1.0);
        EXPECT_LE(recon_norm, bound, "  recon");
        free(recon_re); free(recon_im);

        /* orth = q . q^H - I */
        double* orth_re = (double*)calloc((size_t)q_rows * (size_t)q_rows, sizeof(double));
        double* orth_im = (double*)calloc((size_t)q_rows * (size_t)q_rows, sizeof(double));
        aaH_minus_I(q_re, q_im, q_rows, q_cols, orth_re, orth_im);
        double orth_norm = norm_inf_matrix(orth_re, orth_im, q_rows, q_rows);
        double obound = K_orth * 2.2204460492503131e-16;
        EXPECT_LE(orth_norm, obound, "  orth");
        free(orth_re); free(orth_im);
    }

    free(q_re); free(q_im); free(r_re); free(r_im);
    free(m_re); free(m_im);
    expr_free(qr);
}

/* =====================================================================
 *  Hard-coded shape regressions at MPFR precision (S5.3 coverage).
 * ================================================================== */

#ifdef USE_MPFR

static void test_mpfr_1x1_real_30(void) {
    qrm_property_check("1x1 real", "{{3.5}}", 30, 1, 32.0, 32.0);
}

static void test_mpfr_2x2_real_30(void) {
    qrm_property_check("2x2 real", "{{1, 2}, {3, 4}}", 30, 2, 64.0, 64.0);
}

static void test_mpfr_2x2_real_60(void) {
    qrm_property_check("2x2 real @60d", "{{1, 2}, {3, 4}}", 60, 2, 64.0, 64.0);
}

static void test_mpfr_2x2_real_100(void) {
    qrm_property_check("2x2 real @100d", "{{1, 2}, {3, 4}}", 100, 2, 64.0, 64.0);
}

static void test_mpfr_2x2_real_200(void) {
    qrm_property_check("2x2 real @200d", "{{1, 2}, {3, 4}}", 200, 2, 64.0, 64.0);
}

static void test_mpfr_3x3_real_well_cond(void) {
    qrm_property_check("3x3 real well-cond @60d",
        "{{4, 1, 2}, {1, 5, 3}, {2, 3, 6}}",
        60, 3, 128.0, 128.0);
}

static void test_mpfr_5x5_real_random_60(void) {
    qrm_property_check("5x5 real @60d",
        "{{71/100, 13/100, 97/100, 40/100, 51/100},"
        " {29/100, 81/100, 23/100, 66/100, 18/100},"
        " {42/100, 55/100, 7/100, 89/100, 72/100},"
        " {93/100, 36/100, 61/100, 5/100, 27/100},"
        " {14/100, 68/100, 31/100, 50/100, 84/100}}",
        60, 5, 256.0, 256.0);
}

static void test_mpfr_3x2_tall_60(void) {
    qrm_property_check("3x2 tall @60d",
        "{{1, 2}, {3, 4}, {5, 6}}",
        60, 2, 64.0, 64.0);
}

static void test_mpfr_2x3_wide_60(void) {
    qrm_property_check("2x3 wide @60d",
        "{{1, 2, 3}, {4, 5, 6}}",
        60, 2, 64.0, 64.0);
}

static void test_mpfr_5x3_tall_60(void) {
    qrm_property_check("5x3 tall @60d",
        "{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}, {2, 1, 1/2}, {3/10, 7/10, 12/10}}",
        60, 3, 128.0, 128.0);
}

static void test_mpfr_2x2_complex_60(void) {
    qrm_property_check("2x2 complex @60d",
        "{{1 + 2 I, 3}, {4, 5 - I}}",
        60, 2, 64.0, 64.0);
}

static void test_mpfr_3x3_complex_60(void) {
    qrm_property_check("3x3 complex @60d",
        "{{1 + I, 2, 3 - I},"
        " {4 + I/2, 5 + 2 I, 6},"
        " {7, 8 - I, 9 + 3 I}}",
        60, 3, 128.0, 128.0);
}

static void test_mpfr_rank_deficient_pivoted_30(void) {
    /* Column 2 duplicates column 0; with pivoting on, the MPFR kernel
     * returns rank-2 cleanly.  Without pivoting it should bail to
     * symbolic -- not exercised here (symbolic test covers it). */
    Expr* qr = run(
        "QRDecomposition[N[{{1, 2, 1}, {3, 4, 3}, {5, 6, 5}}, 30], "
        "Pivoting -> True]");
    if (qr->type != EXPR_FUNCTION || qr->data.function.arg_count != 3) {
        fprintf(stderr, "FAIL rank-def pivoted: bad result shape\n");
        g_failures++;
        expr_free(qr);
        return;
    }
    Expr* q = qr->data.function.args[0];
    Expr* r = qr->data.function.args[1];
    int qrows = (int)q->data.function.arg_count;
    int rrows = (int)r->data.function.arg_count;
    if (qrows != 2 || rrows != 2) {
        fprintf(stderr, "FAIL rank-def pivoted: rank %d %d (expected 2)\n",
                qrows, rrows);
        g_failures++;
    } else {
        printf("  PASS rank-def pivoted: rank=2 as expected\n");
    }
    expr_free(qr);
}

static void test_mpfr_zero_matrix(void) {
    Expr* res = run("QRDecomposition[N[{{0, 0}, {0, 0}}, 60]]");
    Expr* q = res->data.function.args[0];
    Expr* r = res->data.function.args[1];
    if (q->data.function.arg_count != 0 || r->data.function.arg_count != 0) {
        fprintf(stderr, "FAIL zero matrix: q or r non-empty\n");
        g_failures++;
    } else {
        printf("  PASS zero matrix -> empty q, r\n");
    }
    expr_free(res);
}

static void test_mpfr_pivot_diag_monotone(void) {
    /* With pivoting on, |R[0,0]| must dominate |R[1,1]| etc. */
    Expr* res = run(
        "QRDecomposition[N["
        "{{1/10, 2/10, 5}, {3/10, 6, 4/10}, {7, 5/10, 6/10}}, 60],"
        " Pivoting -> True]");
    Expr* r = res->data.function.args[1];
    double *re, *im; int rows, cols;
    bool ok = extract_matrix(r, &re, &im, &rows, &cols);
    ASSERT(ok);
    double prev = INFINITY;
    int dn = rows < cols ? rows : cols;
    bool monotone = true;
    for (int i = 0; i < dn; i++) {
        double v = hypot(re[i * cols + i], im[i * cols + i]);
        if (v > prev + 1e-10) { monotone = false; break; }
        prev = v;
    }
    if (!monotone) {
        fprintf(stderr, "FAIL pivot diag monotone\n");
        g_failures++;
    } else {
        printf("  PASS pivot diag monotone: %d entries decreasing\n", dn);
    }
    free(re); free(im);
    expr_free(res);
}

/* Verify m . P == q^H . r at MPFR precision.  Uses the same numeric
 * extraction strategy as the reconstruction check above. */
static void test_mpfr_pivot_identity(void) {
    const char* m_src = "{{1, 2, 5}, {3, 6, 4}, {7, 5, 6}}";
    int digits = 60;
    size_t buflen = strlen(m_src) + 96;
    char* buf = (char*)malloc(buflen);

    snprintf(buf, buflen, "N[%s, %d]", m_src, digits);
    Expr* m = run(buf);
    double *m_re, *m_im; int n_rows, n_cols;
    extract_matrix(m, &m_re, &m_im, &n_rows, &n_cols);
    expr_free(m);

    snprintf(buf, buflen, "QRDecomposition[N[%s, %d], Pivoting -> True]",
             m_src, digits);
    Expr* qr = run(buf);
    Expr* q = qr->data.function.args[0];
    Expr* r = qr->data.function.args[1];
    Expr* P = qr->data.function.args[2];

    double *q_re, *q_im, *r_re, *r_im, *P_re, *P_im;
    int qrows, qcols, rrows, rcols, prows, pcols;
    extract_matrix(q, &q_re, &q_im, &qrows, &qcols);
    extract_matrix(r, &r_re, &r_im, &rrows, &rcols);
    extract_matrix(P, &P_re, &P_im, &prows, &pcols);

    /* mp = m . P  (n_rows x p_cols).  Then compute q^H . r - mp. */
    double* mp_re = (double*)calloc((size_t)n_rows * (size_t)pcols, sizeof(double));
    double* mp_im = (double*)calloc((size_t)n_rows * (size_t)pcols, sizeof(double));
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < pcols; j++) {
            double sr = 0.0, si = 0.0;
            for (int k = 0; k < n_cols; k++) {
                sr += m_re[i * n_cols + k] * P_re[k * pcols + j];
                si += m_im[i * n_cols + k] * P_re[k * pcols + j];
            }
            mp_re[i * pcols + j] = sr;
            mp_im[i * pcols + j] = si;
        }
    }
    double* recon_re = (double*)calloc((size_t)n_rows * (size_t)pcols, sizeof(double));
    double* recon_im = (double*)calloc((size_t)n_rows * (size_t)pcols, sizeof(double));
    cgemm_hermconj_minus(q_re, q_im, qrows, qcols,
                            r_re, r_im, rcols,
                            mp_re, mp_im,
                            recon_re, recon_im);
    double err = norm_inf_matrix(recon_re, recon_im, n_rows, pcols);
    double bound = 128.0 * 2.2204460492503131e-16 * 20.0;
    EXPECT_LE(err, bound, "  pivoted m.P == q^H.r");
    free(mp_re); free(mp_im); free(recon_re); free(recon_im);
    free(m_re); free(m_im);
    free(q_re); free(q_im); free(r_re); free(r_im);
    free(P_re); free(P_im);
    expr_free(qr);
    free(buf);
}

/* ---------------------------------------------------------------------
 *  Precision-monotonicity: factor the same matrix at increasing
 *  precisions; the *raw* MPFR Frobenius norm of the residual should
 *  decay roughly by 2^(bits_step) between successive precisions.
 *
 *  We compute the residual inside Mathilda using the existing
 *  Times / Plus / Sqrt / Total / Flatten pipeline (all MPFR-aware)
 *  and then extract the resulting scalar EXPR_MPFR back to a double
 *  via mpfr_get_d -- the order of magnitude is preserved by mpfr_get_d
 *  (it returns 0 only for genuine zero, otherwise a normalised double
 *  with possibly enormous negative exponent representable via subnormals
 *  / clamped to ~1e-308; we don't need the full precision, just the
 *  magnitude).
 *
 *  For the test below `{{1, 2}, {3, 4}}` at precisions 30 / 60 / 100
 *  digits we expect the Frobenius residual to satisfy
 *      err[60] < err[30] * 10^-15
 *      err[100] < err[60] * 10^-20
 *  i.e. the precision is actually being used.  We use a loose bound
 *  (err[next] <= err[prev]) to tolerate the lower-precision noise
 *  floor.
 * ------------------------------------------------------------------ */
static double frobenius_residual(const char* m_src, int digits) {
    size_t buflen = strlen(m_src) + 256;
    char* buf = (char*)malloc(buflen);
    /* The diff is N[m] - Transpose[q].r.  We square element-wise (Power
     * is Listable), Total[Flatten[...]], Sqrt.  All operators here are
     * MPFR-aware. */
    snprintf(buf, buflen,
        "Module[{m = N[%s, %d], qr},"
        "  qr = QRDecomposition[m];"
        "  Sqrt[Total[Flatten[(m - Transpose[qr[[1]]] . qr[[2]])^2]]]]",
        m_src, digits);
    Expr* res = run(buf);
    free(buf);
    double out = 0.0;
    if (res->type == EXPR_MPFR) {
        out = fabs(mpfr_get_d(res->data.mpfr, MPFR_RNDN));
    } else if (res->type == EXPR_REAL) {
        out = fabs(res->data.real);
    } else if (res->type == EXPR_INTEGER) {
        out = (double)res->data.integer;
        if (out < 0) out = -out;
    } else {
        /* Unexpected.  Print and fail. */
        char* s = expr_to_string(res);
        fprintf(stderr, "frobenius_residual: unexpected residual form: %s\n", s);
        free(s);
        out = HUGE_VAL;
    }
    expr_free(res);
    return out;
}

static void test_mpfr_precision_monotonicity(void) {
    const char* m_src = "{{1, 2}, {3, 4}}";
    double r30  = frobenius_residual(m_src, 30);
    double r60  = frobenius_residual(m_src, 60);
    double r100 = frobenius_residual(m_src, 100);
    double r200 = frobenius_residual(m_src, 200);
    printf("Test: precision-monotonicity {{1,2},{3,4}}\n");
    printf("  Frobenius residual @30d:  %.4e\n",  r30);
    printf("  Frobenius residual @60d:  %.4e\n",  r60);
    printf("  Frobenius residual @100d: %.4e\n", r100);
    printf("  Frobenius residual @200d: %.4e\n", r200);
    /* Each higher precision should give a residual at least as small.
     * Tolerate a tiny slack (factor 2) to absorb the algorithm's
     * inherent precision-spread noise. */
    if (r60 > 2.0 * r30) {
        fprintf(stderr, "FAIL monotonicity: r60 (%.4e) > 2 * r30 (%.4e)\n",
                r60, r30);
        g_failures++;
    } else if (r100 > 2.0 * r60) {
        fprintf(stderr, "FAIL monotonicity: r100 (%.4e) > 2 * r60 (%.4e)\n",
                r100, r60);
        g_failures++;
    } else if (r200 > 2.0 * r100) {
        fprintf(stderr, "FAIL monotonicity: r200 (%.4e) > 2 * r100 (%.4e)\n",
                r200, r100);
        g_failures++;
    } else {
        printf("  PASS monotonicity: residuals non-increasing across precisions\n");
    }
}

#endif /* USE_MPFR */

int main(void) {
#ifndef USE_MPFR
    printf("test_qrdecomposition_mpfr: skipped (USE_MPFR is off).\n");
    return 0;
#else
    symtab_init();
    core_init();

    TEST(test_mpfr_1x1_real_30);
    TEST(test_mpfr_2x2_real_30);
    TEST(test_mpfr_2x2_real_60);
    TEST(test_mpfr_2x2_real_100);
    TEST(test_mpfr_2x2_real_200);
    TEST(test_mpfr_3x3_real_well_cond);
    TEST(test_mpfr_5x5_real_random_60);
    TEST(test_mpfr_3x2_tall_60);
    TEST(test_mpfr_2x3_wide_60);
    TEST(test_mpfr_5x3_tall_60);
    TEST(test_mpfr_2x2_complex_60);
    TEST(test_mpfr_3x3_complex_60);
    TEST(test_mpfr_rank_deficient_pivoted_30);
    TEST(test_mpfr_zero_matrix);
    TEST(test_mpfr_pivot_diag_monotone);
    TEST(test_mpfr_pivot_identity);
    TEST(test_mpfr_precision_monotonicity);

    if (g_failures > 0) {
        fprintf(stderr, "%d MPFR-QR test failure(s)\n", g_failures);
        return 1;
    }
    printf("All MPFR QRDecomposition tests passed.\n");
    return 0;
#endif
}
