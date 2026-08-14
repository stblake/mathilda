/* mlutil.c -- shared plumbing for the src/ml builtins. See mlutil.h. */
#include <stdlib.h>
#include "expr.h"
#include "sym_names.h"
#include "numarray.h"
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
