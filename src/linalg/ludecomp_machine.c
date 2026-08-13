/* src/linalg/ludecomp_machine.c
 *
 * Machine-precision LU decomposition kernel.
 *
 * Invoked by lu_dispatch when common_scan_inexact reports an inexact
 * input at <= 53 bits of precision (i.e. the matrix is dominated by
 * IEEE doubles).  On any failure -- USE_LAPACK off, a leaf the loader
 * can't reduce to a double, a fatal LAPACK info -- it returns NULL and
 * lu_dispatch falls back to the symbolic Doolittle path.
 *
 * High-level flow:
 *
 *   1. lu_mach_load_matrix   -- walk the n x n list-of-list expression
 *                               and copy every leaf into a column-major
 *                               double buffer A_cm of length n*n (real)
 *                               or 2*n*n (complex; interleaved re,im).
 *
 *   2. dlange / zlange       -- one-shot L-infinity norm of the
 *                               original matrix (needed by *gecon).
 *
 *   3. dgetrf / zgetrf       -- in-place factorisation; on success the
 *                               combined Doolittle L (strict-lower,
 *                               unit diag) and U (upper) live in A_cm,
 *                               and ipiv carries the 1-indexed pivot
 *                               vector.
 *
 *   4. dgecon / zgecon       -- reciprocal L-infinity condition
 *                               estimate, anorm = the norm from step 2.
 *
 *   5. Build outputs         -- LU re-emitted as a row-major
 *                               List-of-List of Real / Complex entries,
 *                               perm derived from ipiv (LAPACK's ipiv is
 *                               a series of pairwise swaps; Mathematica's
 *                               p is the resulting permutation), c =
 *                               1 / rcond.
 *
 * Memory contract.  Standard builtin contract: this file does NOT call
 * expr_free on the input `m` -- the evaluator owns it.  Every malloc
 * not handed to the caller is matched by a free along every exit path.
 */

#include "ludecomp_internal.h"
#include "linalg.h"
#include "lapack.h"
#include "ndarray.h"
#include "numarray.h"   /* na_load_matrix / na_build_matrix -- packed I/O */
#include "ndlinalg.h"   /* linalg_delist_and_reeval -- packed fast-path net */
#include "pack.h"       /* ndbuild_open -- packed buffer LU result */
#include "expr.h"
#include "print.h"
#include "sym_names.h"
#include "common.h"

#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_LAPACK
static void lu_mach_warn_once(uint64_t* counter, const char* msg)
{
    if (*counter) return;
    *counter = 1;
    fprintf(stderr, "%s", msg);
}

/* ---------------------------------------------------------------------
 * Numeric-leaf -> double conversion.  Recognises Integer, BigInt,
 * Real, MPFR, exact Rational[p, q], and Complex[re, im].  Anything
 * else triggers a fall-back to the symbolic kernel.
 *
 * Mirrors mach_leaf_to_double in qrdecomp_machine.c; kept local so the
 * two modules stay independent translation units. */
static bool lu_mach_leaf_to_double(Expr* e, double* out_re, double* out_im)
{
    *out_im = 0.0;
    if (!e) return false;
    switch (e->type) {
        case EXPR_REAL:    *out_re = e->data.real;                    return true;
        case EXPR_INTEGER: *out_re = (double)e->data.integer;         return true;
        case EXPR_BIGINT:  *out_re = mpz_get_d(e->data.bigint);       return true;
#ifdef USE_MPFR
        case EXPR_MPFR:
            *out_re = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
            return true;
#endif
        case EXPR_FUNCTION:
            if (e->data.function.head->type == EXPR_SYMBOL) {
                const char* h = e->data.function.head->data.symbol.name;
                if (h == SYM_Rational && e->data.function.arg_count == 2) {
                    double p, q, dummy;
                    if (lu_mach_leaf_to_double(e->data.function.args[0], &p, &dummy)
                        && lu_mach_leaf_to_double(e->data.function.args[1], &q, &dummy)
                        && q != 0.0) {
                        *out_re = p / q;
                        return true;
                    }
                    return false;
                }
                if (h == SYM_Complex && e->data.function.arg_count == 2) {
                    double r, i, dummy;
                    if (lu_mach_leaf_to_double(e->data.function.args[0], &r, &dummy)
                        && lu_mach_leaf_to_double(e->data.function.args[1], &i, &dummy)) {
                        *out_re = r;
                        *out_im = i;
                        return true;
                    }
                    return false;
                }
            }
            return false;
        default:
            return false;
    }
}

static bool lu_mach_leaf_is_complex(Expr* e)
{
    double r, i;
    return lu_mach_leaf_to_double(e, &r, &i) && i != 0.0;
}

/* Load the rows x cols matrix into a column-major double buffer.
 * Column-major layout: A[i + j * rows] is m[i, j].  Same shape as
 * mach_load_matrix in qrdecomp_machine.c -- duplicated locally to
 * keep the modules independent. */
static bool lu_mach_load_matrix(Expr* m, int rows, int cols,
                                 double** out_A, bool* out_is_complex)
{
    if (m->type != EXPR_FUNCTION) return false;

    bool is_complex = false;
    for (int i = 0; i < rows && !is_complex; i++) {
        Expr* row = m->data.function.args[i];
        if (row->type != EXPR_FUNCTION) return false;
        for (int j = 0; j < cols && !is_complex; j++) {
            if (lu_mach_leaf_is_complex(row->data.function.args[j])) {
                is_complex = true;
            }
        }
    }

    size_t stride = is_complex ? 2 : 1;
    size_t total  = stride * (size_t)rows * (size_t)cols;
    double* A = (double*)malloc(total * sizeof(double));
    if (!A) return false;

    for (int i = 0; i < rows; i++) {
        Expr* row = m->data.function.args[i];
        for (int j = 0; j < cols; j++) {
            double re, im;
            if (!lu_mach_leaf_to_double(row->data.function.args[j], &re, &im)) {
                free(A);
                return false;
            }
            size_t off = stride * ((size_t)i + (size_t)j * (size_t)rows);
            A[off] = re;
            if (is_complex) A[off + 1] = im;
        }
    }

    *out_A = A;
    *out_is_complex = is_complex;
    return true;
}

/* Build a Real / Complex Mathilda scalar from (re, im).  Following the
 * same convention as qrdecomp_machine.c -- a Complex with literal-zero
 * imaginary part collapses to a bare Real. */
static Expr* lu_mach_make_scalar(double re, double im, bool is_complex)
{
    if (!is_complex || im == 0.0) return expr_new_real(re);
    Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
    args[0] = expr_new_real(re);
    args[1] = expr_new_real(im);
    Expr* z = expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
    free(args);
    return z;
}

/* Build the row-major LU output (rows x cols) from the column-major
 * LAPACK buffer.
 *
 * Real data goes straight into a packed float64 buffer via ndbuild_open
 * -- one contiguous copy instead of rows*cols fresh Real Exprs, which for
 * a 400x400 factor is the difference between ~5 ms and ~30 ms.  This is a
 * PACKED LIST (Head List, fully transparent under Part/Tr/Apply), NOT a
 * visible NDArray: a visible NDArray here is not transparent under Apply
 * (`Times @@ Diagonal[lu]` would thread instead of reduce) and would be
 * inconsistent with Inverse / LinearSolve, which return packed Lists.
 *
 * Under the packing threshold (small matrices) ndbuild_open returns NULL
 * and we build the nested List as before.  Complex data has no packed
 * form, so it always takes the nested-List path. */
static Expr* lu_mach_build_lu(const double* A_cm, int rows, int cols,
                               bool is_complex)
{
    if (!is_complex) {
        int64_t dims[2] = { rows, cols };
        void* raw = NULL;
        Expr* packed = ndbuild_open(2, dims, NDT_FLOAT64, &raw);
        if (packed) {
            double* out = (double*)raw;
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    out[(size_t)i * (size_t)cols + (size_t)j] =
                        A_cm[(size_t)i + (size_t)j * (size_t)rows];
            return packed;   /* col-major LAPACK -> row-major packed buffer */
        }
        /* else: under threshold / packing off -- fall through to List. */
    }

    size_t stride = is_complex ? 2 : 1;
    Expr** row_exprs = (Expr**)malloc(sizeof(Expr*) * (size_t)rows);
    for (int i = 0; i < rows; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)cols);
        for (int j = 0; j < cols; j++) {
            size_t off = stride * ((size_t)i + (size_t)j * (size_t)rows);
            double re = A_cm[off];
            double im = is_complex ? A_cm[off + 1] : 0.0;
            elems[j] = lu_mach_make_scalar(re, im, is_complex);
        }
        row_exprs[i] = expr_new_function(
            expr_new_symbol(SYM_List), elems, (size_t)cols);
        free(elems);
    }
    Expr* lu = expr_new_function(
        expr_new_symbol(SYM_List), row_exprs, (size_t)rows);
    free(row_exprs);
    return lu;
}

/* Convert LAPACK's `ipiv` (a series of pairwise swaps, length
 * min(rows, cols)) into the 1-indexed permutation vector p of length
 * `rows` satisfying m[[p]] == l . u.
 *
 * The LAPACK convention is: at step k (1-indexed), row k of the
 * working matrix was swapped with row ipiv[k - 1] (also 1-indexed,
 * referring to a position in [k, rows]).  Applying those swaps in
 * order to the identity vector [1..rows] yields the permutation
 * Mathematica calls `p`. */
static Expr* lu_mach_build_perm(const int* ipiv, int rows, int ipiv_len)
{
    int* perm = (int*)malloc((size_t)rows * sizeof(int));
    for (int i = 0; i < rows; i++) perm[i] = i + 1;
    for (int k = 0; k < ipiv_len; k++) {
        int target = ipiv[k] - 1;
        if (target != k && target >= 0 && target < rows) {
            int tmp = perm[k]; perm[k] = perm[target]; perm[target] = tmp;
        }
    }
    Expr** elems = (Expr**)malloc((size_t)rows * sizeof(Expr*));
    for (int i = 0; i < rows; i++) elems[i] = expr_new_integer(perm[i]);
    Expr* p = expr_new_function(
        expr_new_symbol(SYM_List), elems, (size_t)rows);
    free(elems);
    free(perm);
    return p;
}
#endif /* USE_LAPACK */

/* ---------------------------------------------------------------------
 * Kernel.  Returns NULL on any failure path; the caller (lu_dispatch)
 * treats NULL as "fall back to symbolic" and never frees the input.
 *
 * Supports rectangular input: dgetrf / zgetrf accept separate m, n.
 * The condition estimate (dgecon / zgecon) requires a square matrix
 * and is only computed when rows == cols; for non-square input we
 * return c = exact Integer 0 (the estimate is mathematically
 * undefined for non-square A).
 * ------------------------------------------------------------------ */
Expr* lu_machine_dispatch(Expr* m, int rows, int cols)
{
#ifndef USE_LAPACK
    (void)m; (void)rows; (void)cols;
    return NULL;
#else
    static uint64_t lapack_warn_counter = 0;
    static uint64_t sing_warn_counter   = 0;

    int steps = (rows < cols) ? rows : cols;
    bool square = (rows == cols);

    double* A_cm = NULL;
    bool is_complex = false;
    if (!lu_mach_load_matrix(m, rows, cols, &A_cm, &is_complex)) return NULL;

    /* Snapshot the original-matrix L-infinity norm before dgetrf
     * overwrites A_cm with the factorisation.  Only needed for square
     * input (it feeds dgecon). */
    double anorm = 0.0;
    if (square) {
        anorm = is_complex
            ? mat_lapack_zlange('I', rows, cols, A_cm, rows)
            : mat_lapack_dlange('I', rows, cols, A_cm, rows);
        if (anorm < 0.0) {
            /* LAPACK stub fired -- USE_LAPACK was off at compile time
             * for this translation unit. */
            free(A_cm);
            return NULL;
        }
    }

    int* ipiv = (int*)malloc((size_t)steps * sizeof(int));
    if (!ipiv) { free(A_cm); return NULL; }

    int info = is_complex
        ? mat_lapack_zgetrf(rows, cols, A_cm, rows, ipiv)
        : mat_lapack_dgetrf(rows, cols, A_cm, rows, ipiv);
    if (info < 0) {
        lu_mach_warn_once(&lapack_warn_counter,
            "LUDecomposition: LAPACK fast path unavailable; "
            "falling back to symbolic kernel.\n");
        free(A_cm); free(ipiv);
        return NULL;
    }
    /* info > 0 means a zero pivot was found at U[info, info] -- the
     * matrix is singular but the factorisation still completed.  We
     * surface the standard singular warning and continue, matching
     * Mathematica's LUDecomposition::sing behaviour. */
    if (info > 0) {
        if (!sing_warn_counter) {
            sing_warn_counter = 1;
            char* s = expr_to_string(m);
            fprintf(stderr,
                "LUDecomposition::sing: Matrix %s is singular.\n", s);
            free(s);
        }
    }

    /* Condition-number estimate (square only).  When the matrix is
     * singular dgecon returns rcond = 0 (or very small); we still emit
     * a Real so the output shape is uniform.  For non-square input
     * the condition number is mathematically undefined, so we return
     * exact Integer 0 in the c slot.
     *
     * Badly-conditioned warning (::luc):  if cond_est > 1 / DBL_EPSILON
     * the matrix is so ill-conditioned that LAPACK's factorisation
     * may have lost most of its significant digits.  Mathematica
     * emits an LUDecomposition::luc warning in this case; we mirror
     * it.  The warning is suppressed for matrices that already
     * tripped ::sing (info > 0) to avoid double-warning. */
    static uint64_t luc_warn_counter = 0;
    Expr* c;
    if (!square) {
        c = expr_new_integer(0);
    } else {
        double rcond = 0.0;
        double cond_est;
        if (info == 0) {
            int cinfo = is_complex
                ? mat_lapack_zgecon('I', rows, A_cm, rows, anorm, &rcond)
                : mat_lapack_dgecon('I', rows, A_cm, rows, anorm, &rcond);
            if (cinfo != 0 || rcond <= 0.0) {
                cond_est = HUGE_VAL;
            } else {
                cond_est = 1.0 / rcond;
            }
        } else {
            cond_est = HUGE_VAL;
        }
        if (info == 0 && cond_est > 1.0 / DBL_EPSILON
                      && !luc_warn_counter) {
            luc_warn_counter = 1;
            char* s = expr_to_string(m);
            fprintf(stderr,
                "LUDecomposition::luc: Result for LUDecomposition of "
                "badly conditioned matrix %s may contain significant "
                "numerical errors.\n", s);
            free(s);
        }
        c = expr_new_real(cond_est);
    }

    Expr* lu = lu_mach_build_lu(A_cm, rows, cols, is_complex);
    Expr* p  = lu_mach_build_perm(ipiv, rows, steps);

    free(A_cm);
    free(ipiv);

    Expr** items = (Expr**)malloc(sizeof(Expr*) * 3);
    items[0] = lu; items[1] = p; items[2] = c;
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, 3);
    free(items);
    return result;
#endif /* USE_LAPACK */
}

/* Packed fast path: LUDecomposition of an already-buffered matrix.
 *
 * The ndla_* pattern (see ndlinalg.h): read the buffer in place instead
 * of materialising the packed input to a nested List and re-evaluating.
 * The delist alone cost a 400x400 factorisation ~30 ms in rows*cols fresh
 * Real Exprs -- more than everything else in the call combined -- so this
 * is where LUDecomposition's gap against scipy.linalg.lu_factor lived.
 *
 * Handles the square real float64 case; rectangular, complex, non-float,
 * or no-LAPACK falls back to linalg_delist_and_reeval (the universal
 * safety net), which lands on the List machine path above. The result's
 * lu factor inherits the input's presentation (packed List in -> packed
 * List out; visible NDArray in -> NDArray out), matching the rest of the
 * numeric linear-algebra surface. */
Expr* ndla_ludecomposition(Expr* res)
{
#ifndef USE_LAPACK
    return linalg_delist_and_reeval(res);
#else
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return linalg_delist_and_reeval(res);
    Expr* arg = res->data.function.args[0];
    if (!is_ndarray(arg) || arg->data.ndarray.rank != 2
        || arg->data.ndarray.dims[0] != arg->data.ndarray.dims[1]
        || arg->data.ndarray.dims[0] == 0
        || arg->data.ndarray.dtype != NDT_FLOAT64
        || !mathilda_lapack_probe())
        return linalg_delist_and_reeval(res);

    int n, cc; double* A = NULL;
    if (!na_load_matrix(arg, false, /*colmajor=*/true, &n, &cc, &A))
        return linalg_delist_and_reeval(res);

    /* L-infinity norm of the ORIGINAL matrix, before dgetrf overwrites A. */
    double anorm = mat_lapack_dlange('I', n, n, A, n);
    int* ipiv = (int*)malloc(sizeof(int) * (size_t)n);
    if (anorm < 0.0 || !ipiv) { free(A); free(ipiv); return linalg_delist_and_reeval(res); }

    int info = mat_lapack_dgetrf(n, n, A, n, ipiv);
    if (info < 0) { free(A); free(ipiv); return linalg_delist_and_reeval(res); }

    static uint64_t sing_warn = 0, luc_warn = 0;
    if (info > 0 && !sing_warn) {
        sing_warn = 1;
        char* s = expr_to_string(arg);
        fprintf(stderr, "LUDecomposition::sing: Matrix %s is singular.\n", s);
        free(s);
    }

    double rcond = 0.0, cond_est;
    if (info == 0) {
        int cinfo = mat_lapack_dgecon('I', n, A, n, anorm, &rcond);
        cond_est = (cinfo != 0 || rcond <= 0.0) ? HUGE_VAL : 1.0 / rcond;
        if (cond_est > 1.0 / DBL_EPSILON && !luc_warn) {
            luc_warn = 1;
            char* s = expr_to_string(arg);
            fprintf(stderr,
                "LUDecomposition::luc: Result for LUDecomposition of badly "
                "conditioned matrix %s may contain significant numerical "
                "errors.\n", s);
            free(s);
        }
    } else {
        cond_est = HUGE_VAL;
    }

    Expr* lu = na_build_matrix(A, n, n, false, /*colmajor=*/true);
    if (lu && is_ndarray(lu))                        /* inherit presentation */
        lu->data.ndarray.present_as = arg->data.ndarray.present_as;
    Expr* p = lu_mach_build_perm(ipiv, n, n);
    Expr* c = expr_new_real(cond_est);
    free(A); free(ipiv);

    Expr** items = (Expr**)malloc(sizeof(Expr*) * 3);
    items[0] = lu; items[1] = p; items[2] = c;
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, 3);
    free(items);
    return result;
#endif /* USE_LAPACK */
}
