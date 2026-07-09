/*
 * src/linalg/numarray.c
 *
 * Implementation of the dtype-aware Expr <-> double-buffer marshalling used by
 * the BLAS and LAPACK bridges. See numarray.h for the contract.
 *
 * The numeric-leaf recognition mirrors mach_leaf_to_double in
 * src/linalg/qrdecomp_machine.c (Integer / BigInt / Real / MPFR / Rational /
 * Complex); it is reproduced here rather than shared so the bridges do not
 * take a dependency on the QR kernel's internals.
 */

#include "numarray.h"
#include "ndarray.h"
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ------------------------------------------------------------------ */
/*  Leaf conversion                                                    */
/* ------------------------------------------------------------------ */

bool na_read_scalar(const Expr* e, double* re, double* im)
{
    if (!e) return false;
    *im = 0.0;
    switch (e->type) {
        case EXPR_REAL:    *re = e->data.real;             return true;
        case EXPR_INTEGER: *re = (double)e->data.integer;  return true;
        case EXPR_BIGINT:  *re = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *re = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        case EXPR_FUNCTION: {
            const Expr* h = e->data.function.head;
            if (h && h->type == EXPR_SYMBOL && e->data.function.arg_count == 2) {
                const char* name = h->data.symbol;
                if (strcmp(name, "Rational") == 0) {
                    double p, q, d;
                    if (na_read_scalar(e->data.function.args[0], &p, &d)
                        && na_read_scalar(e->data.function.args[1], &q, &d)
                        && q != 0.0) {
                        *re = p / q;
                        return true;
                    }
                    return false;
                }
                if (strcmp(name, "Complex") == 0) {
                    double r, i, d;
                    if (na_read_scalar(e->data.function.args[0], &r, &d)
                        && na_read_scalar(e->data.function.args[1], &i, &d)) {
                        *re = r;
                        *im = i;
                        return true;
                    }
                    return false;
                }
            }
            return false;
        }
        default:
            return false;
    }
}

/* ------------------------------------------------------------------ */
/*  Shape probing                                                      */
/* ------------------------------------------------------------------ */

static bool na_is_list(const Expr* e)
{
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, "List") == 0;
}

/* Store (re, im) at logical element index `lin` into a flat buffer, honouring
 * the complex stride. */
static void na_put(double* buf, size_t lin, double re, double im, bool cplx)
{
    if (cplx) { buf[2 * lin] = re; buf[2 * lin + 1] = im; }
    else      { buf[lin] = re; }
}

/* ------------------------------------------------------------------ */
/*  Vector load                                                        */
/* ------------------------------------------------------------------ */

bool na_load_vector(const Expr* e, bool want_complex, int* n, double** buf)
{
    if (!e) return false;
    size_t stride = want_complex ? 2 : 1;

    /* NDArray rank-1 fast path (real data only). */
    if (is_ndarray(e)) {
        if (e->data.ndarray.rank != 1) return false;
        int len = (int)e->data.ndarray.dims[0];
        if (len <= 0) return false;
        double* out = (double*)malloc(stride * (size_t)len * sizeof(double));
        if (!out) return false;
        for (int i = 0; i < len; i++)
            na_put(out, (size_t)i, e->data.ndarray.data[i], 0.0, want_complex);
        *n = len; *buf = out;
        return true;
    }

    if (!na_is_list(e)) return false;
    int len = (int)e->data.function.arg_count;
    if (len <= 0) return false;

    double* out = (double*)malloc(stride * (size_t)len * sizeof(double));
    if (!out) return false;
    for (int i = 0; i < len; i++) {
        double re, im;
        if (!na_read_scalar(e->data.function.args[i], &re, &im)
            || (!want_complex && im != 0.0)) {
            free(out);
            return false;
        }
        na_put(out, (size_t)i, re, im, want_complex);
    }
    *n = len; *buf = out;
    return true;
}

/* ------------------------------------------------------------------ */
/*  Matrix load                                                        */
/* ------------------------------------------------------------------ */

/* Column-major index i + j*rows, or row-major index i*cols + j. */
static size_t na_index(int i, int j, int rows, int cols, bool colmajor)
{
    return colmajor ? ((size_t)i + (size_t)j * (size_t)rows)
                    : ((size_t)i * (size_t)cols + (size_t)j);
}

bool na_load_matrix(const Expr* e, bool want_complex, bool colmajor,
                    int* rows, int* cols, double** buf)
{
    if (!e) return false;
    size_t stride = want_complex ? 2 : 1;

    /* NDArray rank-2 fast path (real data, source is row-major). */
    if (is_ndarray(e)) {
        if (e->data.ndarray.rank != 2) return false;
        int r = (int)e->data.ndarray.dims[0];
        int c = (int)e->data.ndarray.dims[1];
        if (r <= 0 || c <= 0) return false;
        double* out = (double*)malloc(stride * (size_t)r * (size_t)c * sizeof(double));
        if (!out) return false;
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++) {
                double v = e->data.ndarray.data[(size_t)i * (size_t)c + (size_t)j];
                na_put(out, na_index(i, j, r, c, colmajor), v, 0.0, want_complex);
            }
        *rows = r; *cols = c; *buf = out;
        return true;
    }

    if (!na_is_list(e)) return false;
    int r = (int)e->data.function.arg_count;
    if (r <= 0) return false;
    if (!na_is_list(e->data.function.args[0])) return false;
    int c = (int)e->data.function.args[0]->data.function.arg_count;
    if (c <= 0) return false;

    double* out = (double*)malloc(stride * (size_t)r * (size_t)c * sizeof(double));
    if (!out) return false;
    for (int i = 0; i < r; i++) {
        const Expr* row = e->data.function.args[i];
        if (!na_is_list(row) || (int)row->data.function.arg_count != c) {
            free(out);                     /* ragged / non-list row */
            return false;
        }
        for (int j = 0; j < c; j++) {
            double re, im;
            if (!na_read_scalar(row->data.function.args[j], &re, &im)
                || (!want_complex && im != 0.0)) {
                free(out);
                return false;
            }
            na_put(out, na_index(i, j, r, c, colmajor), re, im, want_complex);
        }
    }
    *rows = r; *cols = c; *buf = out;
    return true;
}

/* ------------------------------------------------------------------ */
/*  Result construction                                                */
/* ------------------------------------------------------------------ */

Expr* na_scalar(double re, double im)
{
    if (im == 0.0) return expr_new_real(re);
    Expr* args[2] = { expr_new_real(re), expr_new_real(im) };
    return expr_new_function(expr_new_symbol("Complex"), args, 2);
}

Expr* na_build_vector(const double* buf, int n, bool is_complex)
{
    if (n <= 0) return NULL;
    if (!is_complex) {
        double* data = (double*)malloc((size_t)n * sizeof(double));
        if (!data) return NULL;
        memcpy(data, buf, (size_t)n * sizeof(double));
        int64_t dims[1] = { n };
        return expr_new_ndarray(1, dims, data);   /* moves `data` */
    }
    Expr** el = (Expr**)malloc((size_t)n * sizeof(Expr*));
    for (int i = 0; i < n; i++)
        el[i] = na_scalar(buf[2 * i], buf[2 * i + 1]);
    Expr* out = expr_new_function(expr_new_symbol("List"), el, (size_t)n);
    free(el);
    return out;
}

Expr* na_build_matrix(const double* buf, int rows, int cols, bool is_complex,
                      bool colmajor)
{
    if (rows <= 0 || cols <= 0) return NULL;

    if (!is_complex) {
        double* data = (double*)malloc((size_t)rows * (size_t)cols * sizeof(double));
        if (!data) return NULL;
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                data[(size_t)i * (size_t)cols + (size_t)j] =
                    buf[na_index(i, j, rows, cols, colmajor)];
        int64_t dims[2] = { rows, cols };
        return expr_new_ndarray(2, dims, data);    /* moves `data` */
    }

    Expr** row_exprs = (Expr**)malloc((size_t)rows * sizeof(Expr*));
    for (int i = 0; i < rows; i++) {
        Expr** el = (Expr**)malloc((size_t)cols * sizeof(Expr*));
        for (int j = 0; j < cols; j++) {
            size_t off = 2 * na_index(i, j, rows, cols, colmajor);
            el[j] = na_scalar(buf[off], buf[off + 1]);
        }
        row_exprs[i] = expr_new_function(expr_new_symbol("List"), el, (size_t)cols);
        free(el);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), row_exprs, (size_t)rows);
    free(row_exprs);
    return out;
}
