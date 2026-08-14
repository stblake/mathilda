/* mlutil.c -- shared plumbing for the src/ml builtins. See mlutil.h. */
#include <stdlib.h>
#include "expr.h"
#include "sym_names.h"
#include "numarray.h"
#include <math.h>
#include "mlutil.h"

Expr* ml_list_of_reals(const double* v, size_t n) {
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    if (!a) return NULL;
    for (size_t i = 0; i < n; i++) a[i] = expr_new_real(v[i]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), a, n);
    free(a);
    return out;
}

Expr* ml_list_matrix(const double* x, size_t n, size_t dim) {
    Expr** rows = malloc(sizeof(Expr*) * (n ? n : 1));
    if (!rows) return NULL;
    for (size_t i = 0; i < n; i++) rows[i] = ml_list_of_reals(x + i * dim, dim);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, n);
    free(rows);
    return out;
}

bool ml_read_data(Expr* e, size_t* n, size_t* dim, double** buf, bool* was_vector) {
    int r = 0, c = 0; double* b = NULL;
    if (na_load_matrix(e, 0, 0, &r, &c, &b) && r > 0 && c > 0) {
        *n = (size_t)r; *dim = (size_t)c; *buf = b; *was_vector = false;
        return true;
    }
    int len = 0; double* v = NULL;
    if (na_load_vector(e, 0, &len, &v) && len > 0) {
        *n = (size_t)len; *dim = 1; *buf = v; *was_vector = true;
        return true;
    }
    return false;
}

double ml_sqdist(const double* a, const double* b, size_t dim) {
    double s = 0.0;
    for (size_t i = 0; i < dim; i++) { double d = a[i] - b[i]; s += d * d; }
    return s;
}

/* Lower Cholesky factor, in place into `l` (dim x dim, row-major). Returns false if
 * the matrix is not positive definite.
 *
 * Written here rather than routed through mat_lapack_dpotrf deliberately. LAPACK is
 * the right tool for a large factorisation, but `dim` is a feature count -- single
 * digits in every realistic use -- and this runs once per component per EM
 * iteration, where the cost is entirely in the O(n * k * dim^2) E-step below. Twenty
 * lines with no conditional dependency on USE_LAPACK is the better trade at this
 * size; if dim ever grows enough to matter, the call site is one function. */
bool ml_chol(const double* a, double* l, size_t dim) {
    for (size_t i = 0; i < dim * dim; i++) l[i] = 0.0;
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j <= i; j++) {
            double s = a[i * dim + j];
            for (size_t p = 0; p < j; p++) s -= l[i * dim + p] * l[j * dim + p];
            if (i == j) {
                if (!(s > 0.0)) return false;      /* also catches NaN */
                l[i * dim + i] = sqrt(s);
            } else {
                l[i * dim + j] = s / l[j * dim + j];
            }
        }
    }
    return true;
}

/* Squared Mahalanobis distance (x - mu)' S^-1 (x - mu) via forward substitution
 * against the Cholesky factor: solve L y = (x - mu), then return y'y. */
double ml_mahalanobis(const double* l, const double* mu, const double* x,
                             size_t dim, double* y) {
    for (size_t i = 0; i < dim; i++) {
        double s = x[i] - mu[i];
        for (size_t p = 0; p < i; p++) s -= l[i * dim + p] * y[p];
        y[i] = s / l[i * dim + i];
    }
    double q = 0.0;
    for (size_t i = 0; i < dim; i++) q += y[i] * y[i];
    return q;
}

