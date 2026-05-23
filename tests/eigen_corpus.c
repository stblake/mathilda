#include "eigen_corpus.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "print.h"
#include "test_utils.h"

/* --- helpers --------------------------------------------------------- */

/* Read a real / integer leaf into a double.  Used for parsing the
 * scalar entries of an Eigenvalues result List. */
static int leaf_to_double(Expr* e, double* out) {
    if (!e) return 0;
    if (e->type == EXPR_REAL)    { *out = e->data.real;             return 1; }
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer;  return 1; }
    /* Negative-int wrapped under Times[-1, n]. */
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, "Times") == 0
        && e->data.function.arg_count == 2) {
        double a, b;
        if (leaf_to_double(e->data.function.args[0], &a)
            && leaf_to_double(e->data.function.args[1], &b)) {
            *out = a * b;
            return 1;
        }
    }
    return 0;
}

/* --- public API ------------------------------------------------------ */

double corpus_norm_inf_real(const double* A, size_t n) {
    double m = 0.0;
    for (size_t i = 0; i < n; i++) {
        double s = 0.0;
        for (size_t j = 0; j < n; j++) s += fabs(A[i * n + j]);
        if (s > m) m = s;
    }
    return m;
}

static char* matrix_to_input(const char* head, const double* A, size_t n) {
    /* Each entry expands to SetPrecision[X.XXX...e+XX, MachinePrecision]
     * (≈ 60 chars).  %.17e is the round-trip-safe format for a double,
     * but Mathilda's parser promotes 17-digit literals to MPFR; the
     * SetPrecision wrap forces the value back to machine precision so
     * the LAPACK fast path is exercised as the test name implies. */
    size_t cap = strlen(head) + 16 + n * n * 72;
    char* buf = (char*)malloc(cap);
    size_t pos = 0;
    pos += snprintf(buf + pos, cap - pos, "%s[{", head);
    for (size_t i = 0; i < n; i++) {
        if (i) pos += snprintf(buf + pos, cap - pos, ", ");
        pos += snprintf(buf + pos, cap - pos, "{");
        for (size_t j = 0; j < n; j++) {
            if (j) pos += snprintf(buf + pos, cap - pos, ", ");
            pos += snprintf(buf + pos, cap - pos,
                            "SetPrecision[%.17e, MachinePrecision]",
                            A[i * n + j]);
        }
        pos += snprintf(buf + pos, cap - pos, "}");
    }
    snprintf(buf + pos, cap - pos, "}]");
    return buf;
}

char* corpus_matrix_to_eigenvalues_input(const double* A, size_t n) {
    return matrix_to_input("Eigenvalues", A, n);
}

char* corpus_matrix_to_eigenvectors_input(const double* A, size_t n) {
    return matrix_to_input("Eigenvectors", A, n);
}

size_t corpus_eval_eigenvalues_real(const double* A, size_t n,
                                     double** lambdas) {
    char* input = corpus_matrix_to_eigenvalues_input(A, n);
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    free(input);

    ASSERT(r != NULL);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(r->data.function.head->data.symbol, "List") == 0);

    size_t out_count = r->data.function.arg_count;
    *lambdas = (double*)malloc(sizeof(double) * (out_count ? out_count : 1));
    for (size_t i = 0; i < out_count; i++) {
        ASSERT(leaf_to_double(r->data.function.args[i], &(*lambdas)[i]));
    }
    expr_free(r);
    return out_count;
}

size_t corpus_eval_eigenvectors_real(const double* A, size_t n,
                                      double** V) {
    char* input = corpus_matrix_to_eigenvectors_input(A, n);
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    free(input);

    ASSERT(r != NULL);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(r->data.function.head->data.symbol, "List") == 0);

    size_t out_count = r->data.function.arg_count;
    *V = (double*)malloc(sizeof(double) * (out_count ? out_count * n : 1));
    for (size_t i = 0; i < out_count; i++) {
        Expr* row = r->data.function.args[i];
        ASSERT(row->type == EXPR_FUNCTION);
        ASSERT(row->data.function.head->type == EXPR_SYMBOL);
        ASSERT(strcmp(row->data.function.head->data.symbol, "List") == 0);
        ASSERT(row->data.function.arg_count == n);
        for (size_t j = 0; j < n; j++) {
            ASSERT(leaf_to_double(row->data.function.args[j],
                                  &(*V)[i * n + j]));
        }
    }
    expr_free(r);
    return out_count;
}

void corpus_check_real_eigenvalues(const char* label,
                                    const char* input,
                                    const double* expected, size_t n,
                                    double tol) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    ASSERT(r != NULL);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(r->data.function.head->data.symbol, "List") == 0);
    ASSERT(r->data.function.arg_count == n);
    int ok = 1;
    for (size_t i = 0; i < n; i++) {
        double got;
        if (!leaf_to_double(r->data.function.args[i], &got)) {
            ok = 0; break;
        }
        if (fabs(got - expected[i]) > tol) {
            printf("FAIL: %s\n  index %zu: expected %.17g, got %.17g (tol %g)\n",
                   label, i, expected[i], got, tol);
            ok = 0;
        }
    }
    if (ok) printf("PASS: %s\n", label);
    expr_free(r); expr_free(e);
    ASSERT(ok);
}

double corpus_assert_residual_real(const double* A, size_t n,
                                    double lambda, const double* v,
                                    double tol) {
    double res = 0.0;
    for (size_t i = 0; i < n; i++) {
        double Avi = 0.0;
        for (size_t j = 0; j < n; j++) Avi += A[i * n + j] * v[j];
        double e = fabs(Avi - lambda * v[i]);
        if (e > res) res = e;
    }
    if (res > tol) {
        printf("FAIL: residual %g exceeds tol %g (lambda=%g)\n",
               res, tol, lambda);
    }
    ASSERT(res <= tol);
    return res;
}

void corpus_assert_orthonormal_real(const double* V, size_t n, double tol) {
    /* V is row-major n x n, row i = i-th eigenvector.  V V^T should
     * be the identity for an orthonormal basis. */
    double worst = 0.0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            double dot = 0.0;
            for (size_t k = 0; k < n; k++) dot += V[i * n + k] * V[j * n + k];
            double target = (i == j) ? 1.0 : 0.0;
            double e = fabs(dot - target);
            if (e > worst) worst = e;
        }
    }
    if (worst > tol) {
        printf("FAIL: orthonormality residual %g exceeds tol %g\n",
               worst, tol);
    }
    ASSERT(worst <= tol);
}

/* === Complex corpus helpers ===========================================
 *
 * Complex matrices are stored as parallel re / im arrays, mirroring the
 * MatD layout in src/mateigen.c so the corpus helpers double as a
 * loose contract test against that internal representation. */

double corpus_norm_inf_complex(const double* A_re, const double* A_im,
                                size_t n) {
    double m = 0.0;
    for (size_t i = 0; i < n; i++) {
        double s = 0.0;
        for (size_t j = 0; j < n; j++) {
            s += hypot(A_re[i * n + j], A_im[i * n + j]);
        }
        if (s > m) m = s;
    }
    return m;
}

/* Render a complex scalar.  Imaginary == 0 -> Real; otherwise Complex[r,i].
 * `%.17e` is round-trip-safe for a double but trips Mathilda's MPFR
 * auto-promotion; the SetPrecision wrap forces the value back to
 * machine precision so the LAPACK numerical path is exercised. */
static int append_complex_entry(char* buf, size_t cap, size_t pos,
                                  double r, double im) {
    if (im == 0.0) {
        return pos + snprintf(buf + pos, cap - pos,
                              "SetPrecision[%.17e, MachinePrecision]", r);
    }
    return pos + snprintf(buf + pos, cap - pos,
                          "Complex[SetPrecision[%.17e, MachinePrecision], "
                          "SetPrecision[%.17e, MachinePrecision]]",
                          r, im);
}

static char* matrix_to_input_complex(const char* head,
                                       const double* A_re,
                                       const double* A_im, size_t n) {
    size_t cap = strlen(head) + 16 + n * n * 144;
    char* buf = (char*)malloc(cap);
    size_t pos = 0;
    pos += snprintf(buf + pos, cap - pos, "%s[{", head);
    for (size_t i = 0; i < n; i++) {
        if (i) pos += snprintf(buf + pos, cap - pos, ", ");
        pos += snprintf(buf + pos, cap - pos, "{");
        for (size_t j = 0; j < n; j++) {
            if (j) pos += snprintf(buf + pos, cap - pos, ", ");
            pos = append_complex_entry(buf, cap, pos,
                                         A_re[i * n + j], A_im[i * n + j]);
        }
        pos += snprintf(buf + pos, cap - pos, "}");
    }
    snprintf(buf + pos, cap - pos, "}]");
    return buf;
}

char* corpus_matrix_to_eigenvalues_input_complex(const double* A_re,
                                                  const double* A_im,
                                                  size_t n) {
    return matrix_to_input_complex("Eigenvalues", A_re, A_im, n);
}

char* corpus_matrix_to_eigenvectors_input_complex(const double* A_re,
                                                   const double* A_im,
                                                   size_t n) {
    return matrix_to_input_complex("Eigenvectors", A_re, A_im, n);
}

/* Extract real + imag parts from a leaf Expr that may be Integer, Real,
 * or Complex[re, im].  Returns 1 on success, 0 if unrecognised. */
static int leaf_to_complex(Expr* e, double* out_re, double* out_im) {
    *out_im = 0.0;
    double v;
    if (leaf_to_double(e, &v)) { *out_re = v; return 1; }
    if (e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, "Complex") == 0
        && e->data.function.arg_count == 2) {
        double r, m;
        if (leaf_to_double(e->data.function.args[0], &r)
            && leaf_to_double(e->data.function.args[1], &m)) {
            *out_re = r; *out_im = m;
            return 1;
        }
    }
    return 0;
}

size_t corpus_eval_hermitian_eigenvalues(const double* A_re,
                                          const double* A_im,
                                          size_t n, double** lambdas) {
    char* input = corpus_matrix_to_eigenvalues_input_complex(A_re, A_im, n);
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    free(input);

    ASSERT(r != NULL);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(r->data.function.head->data.symbol, "List") == 0);

    size_t out_count = r->data.function.arg_count;
    *lambdas = (double*)malloc(sizeof(double) * (out_count ? out_count : 1));
    for (size_t i = 0; i < out_count; i++) {
        double rr, ii;
        ASSERT(leaf_to_complex(r->data.function.args[i], &rr, &ii));
        if (fabs(ii) > 1e-10) {
            printf("FAIL: Hermitian eigenvalue %zu has non-zero imag part %g\n",
                   i, ii);
            ASSERT(0);
        }
        (*lambdas)[i] = rr;
    }
    expr_free(r);
    return out_count;
}

size_t corpus_eval_eigenvectors_complex(const double* A_re,
                                         const double* A_im,
                                         size_t n,
                                         double** V_re, double** V_im) {
    char* input = corpus_matrix_to_eigenvectors_input_complex(A_re, A_im, n);
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    free(input);

    ASSERT(r != NULL);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(r->data.function.head->data.symbol, "List") == 0);

    size_t out_count = r->data.function.arg_count;
    size_t cap = out_count ? out_count * n : 1;
    *V_re = (double*)malloc(sizeof(double) * cap);
    *V_im = (double*)malloc(sizeof(double) * cap);
    for (size_t i = 0; i < out_count; i++) {
        Expr* row = r->data.function.args[i];
        ASSERT(row->type == EXPR_FUNCTION);
        ASSERT(row->data.function.head->type == EXPR_SYMBOL);
        ASSERT(strcmp(row->data.function.head->data.symbol, "List") == 0);
        ASSERT(row->data.function.arg_count == n);
        for (size_t j = 0; j < n; j++) {
            double rr, ii;
            ASSERT(leaf_to_complex(row->data.function.args[j], &rr, &ii));
            (*V_re)[i * n + j] = rr;
            (*V_im)[i * n + j] = ii;
        }
    }
    expr_free(r);
    return out_count;
}

double corpus_assert_residual_complex_real_lambda(const double* A_re,
                                                    const double* A_im,
                                                    size_t n,
                                                    double lambda,
                                                    const double* v_re,
                                                    const double* v_im,
                                                    double tol) {
    double res = 0.0;
    for (size_t i = 0; i < n; i++) {
        double Avi_re = 0.0, Avi_im = 0.0;
        for (size_t j = 0; j < n; j++) {
            double ar = A_re[i * n + j];
            double ai = A_im[i * n + j];
            double vr = v_re[j];
            double vi = v_im[j];
            Avi_re += ar * vr - ai * vi;
            Avi_im += ar * vi + ai * vr;
        }
        double dr = Avi_re - lambda * v_re[i];
        double di = Avi_im - lambda * v_im[i];
        double e = hypot(dr, di);
        if (e > res) res = e;
    }
    if (res > tol) {
        printf("FAIL: complex residual %g exceeds tol %g (lambda=%g)\n",
               res, tol, lambda);
    }
    ASSERT(res <= tol);
    return res;
}

double corpus_assert_residual_complex(const double* A_re,
                                       const double* A_im,
                                       size_t n,
                                       double lambda_re, double lambda_im,
                                       const double* v_re,
                                       const double* v_im,
                                       double tol) {
    double res = 0.0;
    for (size_t i = 0; i < n; i++) {
        double Avi_re = 0.0, Avi_im = 0.0;
        for (size_t j = 0; j < n; j++) {
            double ar = A_re[i * n + j];
            double ai = A_im[i * n + j];
            double vr = v_re[j];
            double vi = v_im[j];
            Avi_re += ar * vr - ai * vi;
            Avi_im += ar * vi + ai * vr;
        }
        /* lambda * v[i] */
        double lvr = lambda_re * v_re[i] - lambda_im * v_im[i];
        double lvi = lambda_re * v_im[i] + lambda_im * v_re[i];
        double dr = Avi_re - lvr;
        double di = Avi_im - lvi;
        double e = hypot(dr, di);
        if (e > res) res = e;
    }
    if (res > tol) {
        printf("FAIL: complex residual %g exceeds tol %g (lambda=%g+%gi)\n",
               res, tol, lambda_re, lambda_im);
    }
    ASSERT(res <= tol);
    return res;
}

size_t corpus_eval_eigenvalues_complex(const double* A_re,
                                        const double* A_im,
                                        size_t n,
                                        double** eval_re,
                                        double** eval_im) {
    char* input = corpus_matrix_to_eigenvalues_input_complex(A_re, A_im, n);
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    free(input);

    ASSERT(r != NULL);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(r->data.function.head->data.symbol, "List") == 0);

    size_t out_count = r->data.function.arg_count;
    *eval_re = (double*)malloc(sizeof(double) * (out_count ? out_count : 1));
    *eval_im = (double*)malloc(sizeof(double) * (out_count ? out_count : 1));
    for (size_t i = 0; i < out_count; i++) {
        double rr, ii;
        ASSERT(leaf_to_complex(r->data.function.args[i], &rr, &ii));
        (*eval_re)[i] = rr;
        (*eval_im)[i] = ii;
    }
    expr_free(r);
    return out_count;
}

void corpus_assert_unitary(const double* V_re, const double* V_im,
                            size_t n, double tol) {
    double worst = 0.0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            /* <V_i, V_j> = sum_k conj(V_ik) V_jk */
            double dr = 0.0, di = 0.0;
            for (size_t k = 0; k < n; k++) {
                double vi_r = V_re[i * n + k], vi_i = V_im[i * n + k];
                double vj_r = V_re[j * n + k], vj_i = V_im[j * n + k];
                /* conj(V_ik) * V_jk = (vi_r - i vi_i)(vj_r + i vj_i) */
                dr += vi_r * vj_r + vi_i * vj_i;
                di += vi_r * vj_i - vi_i * vj_r;
            }
            double target_re = (i == j) ? 1.0 : 0.0;
            double e = hypot(dr - target_re, di);
            if (e > worst) worst = e;
        }
    }
    if (worst > tol) {
        printf("FAIL: unitarity residual %g exceeds tol %g\n", worst, tol);
    }
    ASSERT(worst <= tol);
}
