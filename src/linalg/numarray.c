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
                const char* name = h->data.symbol.name;
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
        && strcmp(e->data.function.head->data.symbol.name, "List") == 0;
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

    /* NDArray rank-1 fast path (all dtypes; float32 widens to double, complex
     * marshals into the interleaved buffer). */
    if (is_ndarray(e)) {
        if (e->data.ndarray.rank != 1) return false;
        int len = (int)e->data.ndarray.dims[0];
        if (len <= 0) return false;
        NDType dt = e->data.ndarray.dtype;
        double* out = (double*)malloc(stride * (size_t)len * sizeof(double));
        if (!out) return false;

        /* THE COPY IS ALREADY THE RIGHT BYTES. A float64 buffer asked for as
         * real doubles is bit-identical to what the loop below would write, so
         * the loop is 10^6 indirect ndt_get calls producing a memcpy. This is
         * the same boundary cost performance.md §2 removed from na_load_matrix
         * ("na_load_matrix converted row-major to column-major one ndt_get at a
         * time; it is now a cache-blocked transpose, or a memcpy where no
         * transpose is needed") -- the vector loader was simply never given the
         * same treatment, and every BLAS bridge entry point, Norm and Normalize
         * pay for it. Measured on Norm of a 10^6 vector: 4.76 ms against
         * np.linalg.norm's 242 us, of which the marshalling was most.
         *
         * ndt_get sets im = 0.0 unconditionally for a real dtype, so the
         * `im != 0.0` rejection below cannot fire here and nothing is skipped. */
        if (dt == NDT_FLOAT64 && !want_complex) {
            memcpy(out, e->data.ndarray.data, (size_t)len * sizeof(double));
            *n = len; *buf = out;
            return true;
        }

        for (int i = 0; i < len; i++) {
            double re, im;
            ndt_get(e->data.ndarray.data, (size_t)i, dt, &re, &im);
            if (!want_complex && im != 0.0) { free(out); return false; }
            na_put(out, (size_t)i, re, im, want_complex);
        }
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
/* Cache-blocked transpose of an r x c row-major float64 matrix into a
 * column-major buffer (equivalently: out[i + j*r] = in[i*c + j]).
 *
 * Tiled for the same reason ndstruct_transpose is: a naive i/j double loop
 * strides the full matrix on one side and thrashes cache, while a 32x32 f64
 * tile is 8 KB per window -- comfortably inside L1 -- so each tile is read and
 * written with unit stride on one side and short hops on the other. */
#define NA_TRANSPOSE_TILE 32

static void na_transpose_f64(const double* in, double* out, int r, int c)
{
    for (int i0 = 0; i0 < r; i0 += NA_TRANSPOSE_TILE) {
        int imax = (i0 + NA_TRANSPOSE_TILE < r) ? i0 + NA_TRANSPOSE_TILE : r;
        for (int j0 = 0; j0 < c; j0 += NA_TRANSPOSE_TILE) {
            int jmax = (j0 + NA_TRANSPOSE_TILE < c) ? j0 + NA_TRANSPOSE_TILE : c;
            for (int i = i0; i < imax; i++)
                for (int j = j0; j < jmax; j++)
                    out[(size_t)i + (size_t)j * (size_t)r] =
                        in[(size_t)i * (size_t)c + (size_t)j];
        }
    }
}

/* The inverse of na_transpose_f64: a column-major rows x cols source into a
 * row-major destination (out[i*cols + j] = in[i + j*rows]). Same tiling. */
static void na_untranspose_f64(const double* in, double* out, int rows, int cols)
{
    for (int i0 = 0; i0 < rows; i0 += NA_TRANSPOSE_TILE) {
        int imax = (i0 + NA_TRANSPOSE_TILE < rows) ? i0 + NA_TRANSPOSE_TILE : rows;
        for (int j0 = 0; j0 < cols; j0 += NA_TRANSPOSE_TILE) {
            int jmax = (j0 + NA_TRANSPOSE_TILE < cols) ? j0 + NA_TRANSPOSE_TILE : cols;
            for (int i = i0; i < imax; i++)
                for (int j = j0; j < jmax; j++)
                    out[(size_t)i * (size_t)cols + (size_t)j] =
                        in[(size_t)i + (size_t)j * (size_t)rows];
        }
    }
}

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

    /* NDArray rank-2 fast path (all dtypes; source is row-major). */
    if (is_ndarray(e)) {
        if (e->data.ndarray.rank != 2) return false;
        int r = (int)e->data.ndarray.dims[0];
        int c = (int)e->data.ndarray.dims[1];
        if (r <= 0 || c <= 0) return false;
        NDType dt = e->data.ndarray.dtype;
        double* out = (double*)malloc(stride * (size_t)r * (size_t)c * sizeof(double));
        if (!out) return false;

        /* float64 -> real double is the overwhelmingly common case and needs no
         * marshalling at all, only a layout change: row-major to column-major IS
         * a transpose, and same-layout is a memcpy.
         *
         * Worth specialising because this loop runs once per element on the way
         * into every LAPACK call. Profiling LinearSolve at 1000x1000 put
         * na_load_matrix + ndt_get at ~38% of the non-idle samples, against the
         * factorisation itself -- an out-of-line ndt_get per element, plus a
         * cache-hostile scattered write when transposing. The values are
         * identical either way; only the traversal changes. */
        if (!want_complex && dt == NDT_FLOAT64) {
            const double* src = (const double*)e->data.ndarray.data;
            if (!colmajor) {
                memcpy(out, src, (size_t)r * (size_t)c * sizeof(double));
            } else {
                na_transpose_f64(src, out, r, c);
            }
            *rows = r; *cols = c; *buf = out;
            return true;
        }

        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++) {
                double re, im;
                ndt_get(e->data.ndarray.data,
                        (size_t)i * (size_t)c + (size_t)j, dt, &re, &im);
                if (!want_complex && im != 0.0) { free(out); return false; }
                na_put(out, na_index(i, j, r, c, colmajor), re, im, want_complex);
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

/* Normalise negative zero in a machine linalg result.
 *
 * LAPACK's factorisations produce -0.0 wherever a zero entry came out of a
 * subtraction or a sign flip: `Inverse[{{2., 0.}, {0., 2.}}]` came back as
 * {{0.5, -0.0}, {0.0, 0.5}}. The sign of a zero in a matrix result carries no
 * information -- the entry is exactly zero, and which way it was reached is an
 * artefact of the pivoting -- but it PRINTS, and it printed differently from
 * the exact path, which answers the same expression with 0.0.
 *
 * That divergence was invisible while PACK_MIN_ELEMENTS was 250, because a
 * matrix small enough to read was never packed and so never reached LAPACK.
 * Lowering it to 4 on 2026-08-02 put a 2x2 on this path and made the two
 * surfaces print different answers for the same input, which is exactly what
 * pack.h's "a representation may never change a value" exists to prevent.
 *
 * Only ever turns -0.0 into 0.0, so no non-zero value can be touched.
 *
 * The assignment is NOT a no-op and should not be "simplified" away: -0.0
 * compares equal to 0.0 and has a different bit pattern, so the store is the
 * whole point. Written as a compare-and-store rather than `x + 0.0` because
 * that identity holds only under round-to-nearest. */
static void na_fix_negative_zero(double* data, size_t n)
{
    for (size_t i = 0; i < n; i++) if (data[i] == 0.0) data[i] = 0.0;
}

Expr* na_build_vector(const double* buf, int n, bool is_complex)
{
    if (n <= 0) return NULL;
    if (!is_complex) {
        double* data = (double*)malloc((size_t)n * sizeof(double));
        if (!data) return NULL;
        memcpy(data, buf, (size_t)n * sizeof(double));
        na_fix_negative_zero(data, (size_t)n);
        int64_t dims[1] = { n };
        return expr_new_ndarray_raw(1, dims, data, NDT_FLOAT64);   /* moves `data` */
    }
    /* Complex result → complex64 NDArray (interleaved re,im matches `buf`),
     * keeping NDArray a closed system under the linalg bridges. */
    double* data = (double*)malloc((size_t)n * 2 * sizeof(double));
    if (!data) return NULL;
    memcpy(data, buf, (size_t)n * 2 * sizeof(double));
    na_fix_negative_zero(data, (size_t)n * 2);
    int64_t dims[1] = { n };
    return expr_new_ndarray_raw(1, dims, data, NDT_COMPLEX64);     /* moves `data` */
}

Expr* na_build_matrix(const double* buf, int rows, int cols, bool is_complex,
                      bool colmajor)
{
    if (rows <= 0 || cols <= 0) return NULL;

    if (!is_complex) {
        double* data = (double*)malloc((size_t)rows * (size_t)cols * sizeof(double));
        if (!data) return NULL;
        /* Same layout change as na_load_matrix, in the other direction: a
         * column-major source is a transpose, a row-major one a memcpy. */
        if (colmajor) {
            na_untranspose_f64(buf, data, rows, cols);
        } else {
            memcpy(data, buf, (size_t)rows * (size_t)cols * sizeof(double));
        }
        na_fix_negative_zero(data, (size_t)rows * (size_t)cols);
        int64_t dims[2] = { rows, cols };
        return expr_new_ndarray_raw(2, dims, data, NDT_FLOAT64);    /* moves `data` */
    }

    /* Complex result → complex64 NDArray (row-major, interleaved re,im),
     * keeping NDArray a closed system under the linalg bridges. */
    double* data = (double*)malloc((size_t)rows * (size_t)cols * 2 * sizeof(double));
    if (!data) return NULL;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) {
            size_t off = 2 * na_index(i, j, rows, cols, colmajor);
            size_t dst = 2 * ((size_t)i * (size_t)cols + (size_t)j);
            data[dst]     = buf[off];
            data[dst + 1] = buf[off + 1];
        }
    na_fix_negative_zero(data, (size_t)rows * (size_t)cols * 2);
    int64_t dims[2] = { rows, cols };
    return expr_new_ndarray_raw(2, dims, data, NDT_COMPLEX64);     /* moves `data` */
}
