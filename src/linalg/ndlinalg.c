/*
 * src/linalg/ndlinalg.c
 *
 * NDArray fast-path implementations for the linear-algebra builtins. See
 * ndlinalg.h for the contract. Every fast path:
 *   - loads the flat double buffer with numarray.c's na_load_* (which has its
 *     own EXPR_NDARRAY fast path: ndt_get straight into a col-/row-major buffer,
 *     widening float32 and marshalling complex),
 *   - computes with an in-house double kernel or a mat_lapack_* routine,
 *   - rebuilds an NDArray result via na_build_* / na_scalar (closed system),
 *   - or defers to linalg_delist_and_reeval(res) for shapes/dtypes it does not
 *     specialise (the routine's ordinary numeric-List path then runs).
 *
 * The real Det/Inverse/LinearSolve/MatrixRank paths use an in-house partial-
 * pivot double LU (nd_lu_real) so they stay O(n^3) in every build, with or
 * without LAPACK. Complex heavy paths use mat_lapack_z* when the runtime probe
 * reports LAPACK, else they defer.
 */

#include "ndlinalg.h"
#include "numarray.h"
#include "ndarray.h"
#include "pack.h"
#include "lapack.h"
#include "eval.h"
#include "print.h"
#include "sym_names.h"
#include <math.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  Generic guard helpers                                              */
/* ------------------------------------------------------------------ */

bool linalg_call_has_ndarray(const Expr* res)
{
    if (!res || res->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        if (is_ndarray(res->data.function.args[i])) return true;
    }
    return false;
}

/* Pack a plain List operand of machine reals, so that a call assembled from
 * SMALL matrices still reaches the NDArray path.
 *
 * linalg_call_has_ndarray asks "is this already a buffer", and a matrix below
 * the packing threshold never is: a 6x6 is 36 elements, well under it, so it
 * stays a nested List and the generic symbolic path runs. That is the whole of
 * the gap in benchmarks/ experiment 22 -- Det on a 6x6 costs 1.28 s over 5000
 * repeats unpacked and 0.0012 s packed, for the identical value. The question
 * worth asking is "is this a matrix of machine numbers", which is this.
 *
 * Deliberately REAL machine buffers only: pack_force sniffs the dtype, and an
 * all-Integer or all-Rational matrix would either pack to an int64 buffer or
 * decline, in both cases putting an exact matrix on a floating-point path whose
 * answer is not the exact one Det and Inverse owe it. Anything that does not
 * sniff to NDT_FLOAT64 is handed back unpacked.
 *
 * Returns a freshly owned copy of `arg` as a packed array, or NULL to decline.
 * Never mutates or frees `arg`. */
Expr* linalg_pack_machine_operand(const Expr* arg)
{
    if (!arg || is_ndarray(arg)) return NULL;
    Expr* packed = pack_force(expr_copy((Expr*)arg), false, NDT_FLOAT64);
    if (!packed) return NULL;
    if (!is_ndarray(packed) || packed->data.ndarray.dtype != NDT_FLOAT64) {
        expr_free(packed);
        return NULL;
    }
    return packed;
}

/* The same, for a whole one-argument call: returns a rebuilt `head[packed]`
 * ready to hand to an ndla_* fast path, or NULL to decline. The caller owns the
 * result and frees it; `res` is untouched. */
Expr* linalg_pack_machine_call1(const Expr* res)
{
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;
    Expr* packed = linalg_pack_machine_operand(res->data.function.args[0]);
    if (!packed) return NULL;
    Expr* args[1] = { packed };
    return expr_new_function(expr_copy((Expr*)res->data.function.head), args, 1);
}

/* Two-argument form, for LinearSolve[m, v] and Dot[a, b]. BOTH operands must be
 * machine-real -- a mixed call (one packed, one exact) would hand the fast path
 * an operand pair whose natural result is not the exact one the generic path
 * would produce, so it is declined rather than half-converted. */
Expr* linalg_pack_machine_call2(const Expr* res)
{
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* a = linalg_pack_machine_operand(res->data.function.args[0]);
    if (!a) return NULL;
    Expr* b = linalg_pack_machine_operand(res->data.function.args[1]);
    if (!b) { expr_free(a); return NULL; }
    Expr* args[2] = { a, b };
    return expr_new_function(expr_copy((Expr*)res->data.function.head), args, 2);
}

/*
 * Give a freshly built array result the PRESENTATION the call's operands imply.
 *
 * na_build_vector / na_build_matrix have no source array to inherit from, so
 * they always produce the visible NDArray[...] form. That was right while the
 * only way into these paths was for the user to type NDArray[...] -- naming the
 * head is a request for machine-buffer semantics, so a machine buffer comes
 * back. Once ordinary Lists pack automatically it is wrong: Inverse of a packed
 * 3x3 answered NDArray[{{...}}] where the same value as a plain List answered
 * {{...}}, which is a visible change of head from an invisible optimisation.
 *
 * `a` and `b` are the call's operands (b may be NULL for a one-argument head);
 * nd_present_src2 applies the same visible-dominates-packed rule the elementwise
 * paths use, so NDArray[m] . packedVector stays visible.
 *
 * ONLY a top-level array is stamped. A head that answers with a LIST of arrays
 * (QRDecomposition, LUDecomposition, SingularValueDecomposition, Eigenvectors)
 * must leave them materialised: List is deliberately not packed-aware -- that is
 * what makes the transparency gate an O(argc) top-level scan instead of an
 * O(tree) walk -- so a buffer sitting inside a plain List would be reachable by
 * an unaware head with no gate between them. Those heads already return plain
 * Lists and are left alone.
 */
static Expr* nd_result_like(const Expr* a, const Expr* b, Expr* out)
{
    if (!out || !is_ndarray(out)) return out;
    const Expr* src = nd_present_src2(a, b);
    if (src) out->data.ndarray.present_as = src->data.ndarray.present_as;
    return out;
}

Expr* linalg_delist_and_reeval(Expr* res)
{
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    Expr** args = (n > 0) ? (Expr**)malloc(sizeof(Expr*) * n) : NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* a = res->data.function.args[i];
        args[i] = is_ndarray(a) ? ndarray_to_nested_list(a) : expr_copy(a);
    }
    Expr* call = expr_new_function(expr_copy(res->data.function.head), args, n);
    if (args) free(args);
    /* The rebuilt call has no NDArray arguments, so evaluating it cannot
     * re-enter a linalg_call_has_ndarray guard: no recursion. */
    return eval_and_free(call);
}

/* Does the sole matrix/vector argument carry a complex dtype? */
static bool nd_arg_is_complex(const Expr* e)
{
    return is_ndarray(e) && ndt_is_complex(e->data.ndarray.dtype);
}

/* ------------------------------------------------------------------ */
/*  In-house real double LU (column-major, partial pivoting)           */
/* ------------------------------------------------------------------ */

/* Factor the n x n column-major matrix A in place. piv[k] is the row swapped
 * into position k. *sign accumulates the permutation sign. Returns the number
 * of nonzero pivots (== n for a nonsingular matrix); a zero pivot leaves the
 * remaining column untouched and does not increment the count. */
static int nd_lu_real(double* A, int n, int* piv, int* sign)
{
    *sign = 1;
    int rank = 0;
    for (int k = 0; k < n; k++) {
        int p = k;
        double mx = fabs(A[k + (size_t)k * n]);
        for (int i = k + 1; i < n; i++) {
            double v = fabs(A[i + (size_t)k * n]);
            if (v > mx) { mx = v; p = i; }
        }
        piv[k] = p;
        if (mx == 0.0) continue;              /* singular column */
        rank++;
        if (p != k) {
            for (int j = 0; j < n; j++) {
                double t = A[k + (size_t)j * n];
                A[k + (size_t)j * n] = A[p + (size_t)j * n];
                A[p + (size_t)j * n] = t;
            }
            *sign = -*sign;
        }
        double akk = A[k + (size_t)k * n];
        for (int i = k + 1; i < n; i++) A[i + (size_t)k * n] /= akk;
        for (int j = k + 1; j < n; j++) {
            double akj = A[k + (size_t)j * n];
            if (akj != 0.0)
                for (int i = k + 1; i < n; i++)
                    A[i + (size_t)j * n] -= A[i + (size_t)k * n] * akj;
        }
    }
    return rank;
}

/* Solve L U x = P b in place on one column-major rhs of length n, using the
 * factors produced by nd_lu_real. */
static void nd_lu_solve_col(const double* LU, int n, const int* piv, double* b)
{
    for (int k = 0; k < n; k++) {                 /* apply P */
        if (piv[k] != k) { double t = b[k]; b[k] = b[piv[k]]; b[piv[k]] = t; }
    }
    for (int i = 0; i < n; i++)                    /* forward: unit lower L */
        for (int k = 0; k < i; k++)
            b[i] -= LU[i + (size_t)k * n] * b[k];
    for (int i = n - 1; i >= 0; i--) {             /* back: upper U */
        for (int k = i + 1; k < n; k++)
            b[i] -= LU[i + (size_t)k * n] * b[k];
        b[i] /= LU[i + (size_t)i * n];
    }
}

/* ------------------------------------------------------------------ */
/*  Det                                                                */
/* ------------------------------------------------------------------ */

/* Build a REAL determinant from the LU U-diagonal.
 *
 * det = sign * prod(U[k,k]).  Accumulating that product in a raw double
 * overflows to +/-inf when the true |det| exceeds DBL_MAX (a 200x200
 * RandomReal[{-10,10}] matrix has |det| ~ 1e340), and underflows to 0.0
 * mid-product even when the final value is representable.  So:
 *   - an exact-zero pivot  -> the matrix is singular          -> Real 0.0
 *   - a finite nonzero product -> the common case             -> Real
 *   - otherwise (overflow / underflow of a non-singular product) -> re-multiply
 *     the SAME diagonal in a 53-bit-mantissa mpfr_t, whose exponent range is
 *     effectively unbounded, and return an arbitrary-precision real.  The
 *     mantissa stays 53 bits because the input is machine precision; only the
 *     exponent is widened, so the answer is exactly what the double product
 *     would be if doubles had no exponent limit.
 * Without USE_MPFR the (possibly inf) double result is returned unchanged.
 *
 * LU is the column-major factored buffer; the U diagonal is LU[k + k*n]. */
static Expr* nd_real_det_result(const double* LU, int n, int sign)
{
    double det = (double)sign;
    bool zero_pivot = false;
    for (int k = 0; k < n; k++) {
        double u = LU[k + (size_t)k * n];
        if (u == 0.0) zero_pivot = true;
        det *= u;
    }
    if (zero_pivot) return expr_new_real(0.0);          /* singular */
    if (isfinite(det) && det != 0.0) return expr_new_real(det);
#ifdef USE_MPFR
    mpfr_t acc;
    mpfr_init2(acc, 53);
    mpfr_set_si(acc, sign, MPFR_RNDN);
    for (int k = 0; k < n; k++)
        mpfr_mul_d(acc, acc, LU[k + (size_t)k * n], MPFR_RNDN);
    return expr_new_mpfr_move(acc);                      /* moves acc */
#else
    return expr_new_real(det);
#endif
}

/* Complex analogue.  The U diagonal is stored interleaved: A[2*(k + k*n)]
 * (real) and A[2*(k + k*n)+1] (imag).  A finite running product is returned as
 * a machine complex; on overflow/underflow of either component the product is
 * recomputed in mpfr and returned as Complex[mpfr, mpfr] (evaluated, so a zero
 * imaginary part collapses to a real). */
static Expr* nd_complex_det_result(const double* A, int n, int sign)
{
    double dr = (double)sign, di = 0.0;
    for (int k = 0; k < n; k++) {
        double ur = A[2 * ((size_t)k + (size_t)k * n)];
        double ui = A[2 * ((size_t)k + (size_t)k * n) + 1];
        double nr = dr * ur - di * ui;
        double ni = dr * ui + di * ur;
        dr = nr; di = ni;
    }
    if (isfinite(dr) && isfinite(di)) return na_scalar(dr, di);
#ifdef USE_MPFR
    mpfr_t ar, ai, ur, ui, nr, ni, t;
    mpfr_inits2(53, ar, ai, ur, ui, nr, ni, t, (mpfr_ptr)0);
    mpfr_set_si(ar, sign, MPFR_RNDN);
    mpfr_set_zero(ai, 1);
    for (int k = 0; k < n; k++) {
        mpfr_set_d(ur, A[2 * ((size_t)k + (size_t)k * n)], MPFR_RNDN);
        mpfr_set_d(ui, A[2 * ((size_t)k + (size_t)k * n) + 1], MPFR_RNDN);
        /* (ar + I ai)(ur + I ui) = (ar*ur - ai*ui) + I (ar*ui + ai*ur) */
        mpfr_mul(nr, ar, ur, MPFR_RNDN); mpfr_mul(t, ai, ui, MPFR_RNDN);
        mpfr_sub(nr, nr, t, MPFR_RNDN);
        mpfr_mul(ni, ar, ui, MPFR_RNDN); mpfr_mul(t, ai, ur, MPFR_RNDN);
        mpfr_add(ni, ni, t, MPFR_RNDN);
        mpfr_set(ar, nr, MPFR_RNDN); mpfr_set(ai, ni, MPFR_RNDN);
    }
    Expr* re = expr_new_mpfr_copy(ar);
    Expr* im = expr_new_mpfr_copy(ai);
    mpfr_clears(ar, ai, ur, ui, nr, ni, t, (mpfr_ptr)0);
    Expr* c = expr_new_function(expr_new_symbol(SYM_Complex),
                                (Expr*[]){ re, im }, 2);
    return eval_and_free(c);
#else
    return na_scalar(dr, di);
#endif
}

Expr* ndla_det(Expr* res)
{
    if (res->data.function.arg_count != 1) return linalg_delist_and_reeval(res);
    Expr* arg = res->data.function.args[0];
    if (!is_ndarray(arg) || arg->data.ndarray.rank != 2
        || arg->data.ndarray.dims[0] != arg->data.ndarray.dims[1]
        || arg->data.ndarray.dims[0] == 0) {
        /* non-square / wrong rank: let the real builtin emit Det::matsq */
        return linalg_delist_and_reeval(res);
    }

    /* Exact-integer buffer: the determinant of an integer matrix is an exact
     * integer that routinely exceeds int64 and the 53-bit double mantissa (a
     * 2x2 with 1e9 entries already has det = 999999999999999999, which a double
     * rounds to 1e18).  Delist to the exact FLINT path, which promotes to a GMP
     * bignum -- rather than computing a lossy float LU. */
    if (arg->data.ndarray.dtype == NDT_INT64) {
        return linalg_delist_and_reeval(res);
    }

    bool cplx = nd_arg_is_complex(arg);

    if (!cplx) {
        int n, cc; double* A = NULL;
        if (!na_load_matrix(arg, false, true, &n, &cc, &A))
            return linalg_delist_and_reeval(res);
        int* piv = (int*)malloc(sizeof(int) * (size_t)n);
        int sign = 1;
        if (mathilda_lapack_probe()) {
            /* LAPACK dgetrf (blocked, vectorized) — det = sign * prod(diag U).
             * The product itself is formed by nd_real_det_result, which keeps a
             * finite answer when the raw double product would over/underflow. */
            int info = mat_lapack_dgetrf(n, n, A, n, piv);
            if (info < 0) { free(piv); free(A); return linalg_delist_and_reeval(res); }
            for (int k = 0; k < n; k++)
                if (piv[k] != k + 1) sign = -sign;    /* ipiv is 1-indexed */
        } else {
            nd_lu_real(A, n, piv, &sign);
        }
        Expr* out = nd_real_det_result(A, n, sign);
        free(piv); free(A);
        return out;
    }

    /* Complex: use LAPACK zgetrf when available, else defer. */
    if (!mathilda_lapack_probe()) return linalg_delist_and_reeval(res);
    int n, cc; double* A = NULL;
    if (!na_load_matrix(arg, true, true, &n, &cc, &A))
        return linalg_delist_and_reeval(res);
    int* piv = (int*)malloc(sizeof(int) * (size_t)n);
    int info = mat_lapack_zgetrf(n, n, A, n, piv);
    if (info < 0) { free(piv); free(A); return linalg_delist_and_reeval(res); }
    int sign = 1;
    for (int k = 0; k < n; k++)
        if (piv[k] != k + 1) sign = -sign;            /* ipiv is 1-indexed */
    Expr* out = nd_complex_det_result(A, n, sign);
    free(piv); free(A);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Inverse                                                            */
/* ------------------------------------------------------------------ */

Expr* ndla_inverse(Expr* res)
{
    if (res->data.function.arg_count != 1) return linalg_delist_and_reeval(res);
    Expr* arg = res->data.function.args[0];
    if (!is_ndarray(arg) || arg->data.ndarray.rank != 2
        || arg->data.ndarray.dims[0] != arg->data.ndarray.dims[1]
        || arg->data.ndarray.dims[0] == 0)
        return linalg_delist_and_reeval(res);
    bool cplx = nd_arg_is_complex(arg);

    if (!cplx) {
        int n, cc; double* A = NULL;
        if (!na_load_matrix(arg, false, true, &n, &cc, &A))
            return linalg_delist_and_reeval(res);
        int* piv = (int*)malloc(sizeof(int) * (size_t)n);
        if (mathilda_lapack_probe()) {
            /* LAPACK dgetrf + dgetri (blocked, vectorized). */
            int info = mat_lapack_dgetrf(n, n, A, n, piv);
            if (info != 0) {                          /* <0 bad arg, >0 singular */
                free(piv); free(A); return linalg_delist_and_reeval(res);
            }
            info = mat_lapack_dgetri(n, A, n, piv);
            free(piv);
            if (info != 0) { free(A); return linalg_delist_and_reeval(res); }
            Expr* out = nd_result_like(arg, NULL,
                                       na_build_matrix(A, n, n, false, true));
            free(A);
            return out;
        }
        int sign, rank = nd_lu_real(A, n, piv, &sign);
        if (rank < n) {                              /* singular: defer for Inverse::sing */
            free(piv); free(A);
            return linalg_delist_and_reeval(res);
        }
        /* Solve A X = I column by column (col-major result). */
        double* X = (double*)malloc(sizeof(double) * (size_t)n * (size_t)n);
        double* col = (double*)malloc(sizeof(double) * (size_t)n);
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) col[i] = (i == j) ? 1.0 : 0.0;
            nd_lu_solve_col(A, n, piv, col);
            for (int i = 0; i < n; i++) X[i + (size_t)j * n] = col[i];
        }
        Expr* out = nd_result_like(arg, NULL,
                                   na_build_matrix(X, n, n, false, true));
        free(col); free(X); free(piv); free(A);
        return out;
    }

    if (!mathilda_lapack_probe()) return linalg_delist_and_reeval(res);
    int n, cc; double* A = NULL;
    if (!na_load_matrix(arg, true, true, &n, &cc, &A))
        return linalg_delist_and_reeval(res);
    int* piv = (int*)malloc(sizeof(int) * (size_t)n);
    int info = mat_lapack_zgetrf(n, n, A, n, piv);
    if (info != 0) { free(piv); free(A); return linalg_delist_and_reeval(res); }
    info = mat_lapack_zgetri(n, A, n, piv);
    free(piv);
    if (info != 0) { free(A); return linalg_delist_and_reeval(res); }
    Expr* out = nd_result_like(arg, NULL,
                               na_build_matrix(A, n, n, true, true));
    free(A);
    return out;
}

/* ------------------------------------------------------------------ */
/*  LinearSolve                                                        */
/* ------------------------------------------------------------------ */

/* Load the rhs (vector or matrix) of a solve into a col-major buffer with n
 * rows. Returns false (and defers via the caller) on shape/dtype mismatch.
 * *is_vec records whether the rhs was rank-1 (so the result is a vector). */
static bool nd_load_rhs(const Expr* b, bool want_complex, int n,
                        int* nrhs, bool* is_vec, double** buf)
{
    int bn = 0, bc = 0;
    if (is_ndarray(b) && b->data.ndarray.rank == 1) {
        if (!na_load_vector(b, want_complex, &bn, buf)) return false;
        if (bn != n) { free(*buf); *buf = NULL; return false; }
        *nrhs = 1; *is_vec = true; return true;
    }
    if (is_ndarray(b) && b->data.ndarray.rank == 2) {
        if (!na_load_matrix(b, want_complex, true, &bn, &bc, buf)) return false;
        if (bn != n) { free(*buf); *buf = NULL; return false; }
        *nrhs = bc; *is_vec = false; return true;
    }
    /*
     * A PLAIN LIST rhs. This used to defer "to keep behaviour uniform", and
     * that deferral was a crash: automatic packing keys on total element count,
     * so LinearSolve[m, b] on a 90x90 system arrives with m packed (8100
     * elements) and b NOT packed (90) -- the single commonest shape there is.
     * Deferring sent a machine-real solve to linalg_delist_and_reeval, whose
     * ordinary path is fraction-free Bareiss elimination with a polynomial GCD
     * per pivot, and it STACK-OVERFLOWED.
     *
     * na_load_vector / na_load_matrix already accept an ordinary numeric List
     * (see numarray.h) and fail cleanly on a symbolic or non-machine leaf, so
     * accepting one here needs no new marshalling and cannot widen what the
     * fast path claims. Vector is tried first: a nested List presents its rows
     * as non-numeric leaves, so na_load_vector rejects it, which discriminates
     * rank 1 from rank 2 without a separate shape probe.
     */
    if (na_load_vector(b, want_complex, &bn, buf)) {
        if (bn != n) { free(*buf); *buf = NULL; return false; }
        *nrhs = 1; *is_vec = true; return true;
    }
    if (na_load_matrix(b, want_complex, true, &bn, &bc, buf)) {
        if (bn != n) { free(*buf); *buf = NULL; return false; }
        *nrhs = bc; *is_vec = false; return true;
    }
    return false;
}

Expr* ndla_linearsolve(Expr* res)
{
    if (res->data.function.arg_count != 2) return linalg_delist_and_reeval(res);
    Expr* m = res->data.function.args[0];
    Expr* b = res->data.function.args[1];
    if (!is_ndarray(m) || m->data.ndarray.rank != 2
        || m->data.ndarray.dims[0] != m->data.ndarray.dims[1]
        || m->data.ndarray.dims[0] == 0)
        return linalg_delist_and_reeval(res);
    bool cplx = nd_arg_is_complex(m) || nd_arg_is_complex(b);

    if (!cplx) {
        int n, cc; double* A = NULL;
        if (!na_load_matrix(m, false, true, &n, &cc, &A))
            return linalg_delist_and_reeval(res);
        int nrhs; bool is_vec; double* B = NULL;
        if (!nd_load_rhs(b, false, n, &nrhs, &is_vec, &B)) {
            free(A); return linalg_delist_and_reeval(res);
        }
        int* piv = (int*)malloc(sizeof(int) * (size_t)n);
        if (mathilda_lapack_probe()) {
            /* LAPACK dgesv (blocked LU + solve). */
            int info = mat_lapack_dgesv(n, nrhs, A, n, piv, B, n);
            free(piv); free(A);
            if (info != 0) { free(B); return linalg_delist_and_reeval(res); }
            Expr* out = nd_result_like(m, b, is_vec
                            ? na_build_vector(B, n, false)
                            : na_build_matrix(B, n, nrhs, false, true));
            free(B);
            return out;
        }
        int sign, rank = nd_lu_real(A, n, piv, &sign);
        if (rank < n) {                              /* singular: defer */
            free(piv); free(A); free(B);
            return linalg_delist_and_reeval(res);
        }
        for (int j = 0; j < nrhs; j++) nd_lu_solve_col(A, n, piv, B + (size_t)j * n);
        Expr* out = nd_result_like(m, b, is_vec
                        ? na_build_vector(B, n, false)
                        : na_build_matrix(B, n, nrhs, false, true));
        free(piv); free(A); free(B);
        return out;
    }

    if (!mathilda_lapack_probe()) return linalg_delist_and_reeval(res);
    int n, cc; double* A = NULL;
    if (!na_load_matrix(m, true, true, &n, &cc, &A))
        return linalg_delist_and_reeval(res);
    int nrhs; bool is_vec; double* B = NULL;
    if (!nd_load_rhs(b, true, n, &nrhs, &is_vec, &B)) {
        free(A); return linalg_delist_and_reeval(res);
    }
    int* piv = (int*)malloc(sizeof(int) * (size_t)n);
    int info = mat_lapack_zgesv(n, nrhs, A, n, piv, B, n);
    free(piv); free(A);
    if (info != 0) { free(B); return linalg_delist_and_reeval(res); }
    Expr* out = nd_result_like(m, b, is_vec
                    ? na_build_vector(B, n, true)
                    : na_build_matrix(B, n, nrhs, true, true));
    free(B);
    return out;
}

/* ------------------------------------------------------------------ */
/*  MatrixRank (real; complex defers)                                  */
/* ------------------------------------------------------------------ */

Expr* ndla_matrixrank(Expr* res)
{
    if (res->data.function.arg_count != 1) return linalg_delist_and_reeval(res);
    Expr* arg = res->data.function.args[0];
    if (!is_ndarray(arg) || arg->data.ndarray.rank != 2 || nd_arg_is_complex(arg))
        return linalg_delist_and_reeval(res);

    int r, c; double* A = NULL;                     /* row-major r x c */
    if (!na_load_matrix(arg, false, false, &r, &c, &A))
        return linalg_delist_and_reeval(res);

    double maxabs = 0.0;
    for (size_t i = 0; i < (size_t)r * (size_t)c; i++) {
        double v = fabs(A[i]); if (v > maxabs) maxabs = v;
    }
    double tol = (double)((r > c) ? r : c) * 2.220446049250313e-16 * maxabs;

    int rank = 0;
    for (int col = 0, row = 0; col < c && row < r; col++) {
        int piv = row; double mx = fabs(A[(size_t)row * c + col]);
        for (int i = row + 1; i < r; i++) {
            double v = fabs(A[(size_t)i * c + col]);
            if (v > mx) { mx = v; piv = i; }
        }
        if (mx <= tol) continue;                     /* no pivot in this column */
        if (piv != row)
            for (int j = 0; j < c; j++) {
                double t = A[(size_t)row * c + j];
                A[(size_t)row * c + j] = A[(size_t)piv * c + j];
                A[(size_t)piv * c + j] = t;
            }
        double d = A[(size_t)row * c + col];
        for (int i = row + 1; i < r; i++) {
            double f = A[(size_t)i * c + col] / d;
            if (f != 0.0)
                for (int j = col; j < c; j++)
                    A[(size_t)i * c + j] -= f * A[(size_t)row * c + j];
        }
        rank++; row++;
    }
    free(A);
    return expr_new_integer(rank);
}

/* ------------------------------------------------------------------ */
/*  Tr (2-D, default Plus over the diagonal)                           */
/* ------------------------------------------------------------------ */

Expr* ndla_tr(Expr* res)
{
    if (res->data.function.arg_count != 1) return linalg_delist_and_reeval(res);
    Expr* arg = res->data.function.args[0];
    if (!is_ndarray(arg) || arg->data.ndarray.rank != 2)
        return linalg_delist_and_reeval(res);
    int r = (int)arg->data.ndarray.dims[0];
    int c = (int)arg->data.ndarray.dims[1];
    int d = (r < c) ? r : c;
    NDType dt = arg->data.ndarray.dtype;
    double sr = 0.0, si = 0.0;
    for (int i = 0; i < d; i++) {
        double re, im;
        ndt_get(arg->data.ndarray.data, (size_t)i * (size_t)c + (size_t)i, dt, &re, &im);
        sr += re; si += im;
    }
    return na_scalar(sr, si);
}

/* ------------------------------------------------------------------ */
/*  Norm (vectors: default 2-norm, integer/real p, or Infinity)        */
/* ------------------------------------------------------------------ */

/* Matrix norm of a rank-2 NDArray. Default (or p==2) is the spectral norm (the
 * largest singular value, via LAPACK gesdd); p==1 / Infinity / "Frobenius" go
 * through LAPACK lange. Returns NULL to signal "defer" (no LAPACK, or an
 * unsupported p). */
Expr* ndla_matrix_norm_direct(Expr* v, Expr* pe)
{
    char kind = 'S';                                  /* S = spectral (2-norm) */
    if (pe) {
        if (pe->type == EXPR_INTEGER && pe->data.integer == 2) kind = 'S';
        else if (pe->type == EXPR_INTEGER && pe->data.integer == 1) kind = '1';
        else if (pe->type == EXPR_SYMBOL && pe->data.symbol.name == SYM_Infinity) kind = 'I';
        else if (pe->type == EXPR_STRING && pe->data.string &&
                 strcmp(pe->data.string, "Frobenius") == 0) kind = 'F';
        else return NULL;
    }
    if (!mathilda_lapack_probe()) return NULL;
    bool cplx = nd_arg_is_complex(v);
    int r, c; double* A = NULL;                        /* column-major r x c */
    if (!na_load_matrix(v, cplx, true, &r, &c, &A)) return NULL;

    double out;
    if (kind == 'S') {
        int mn = r < c ? r : c;
        double* S = (double*)malloc(sizeof(double) * (size_t)mn);
        int info = cplx ? mat_lapack_zgesdd('N', r, c, A, r, S, NULL, 1, NULL, 1)
                        : mat_lapack_dgesdd('N', r, c, A, r, S, NULL, 1, NULL, 1);
        out = (info == 0) ? S[0] : -1.0;               /* gesdd sorts S descending */
        free(S);
        if (info != 0) { free(A); return NULL; }
    } else {
        out = cplx ? mat_lapack_zlange(kind, r, c, A, r)
                   : mat_lapack_dlange(kind, r, c, A, r);
        if (out < 0.0) { free(A); return NULL; }
    }
    free(A);
    return expr_new_real(out);
}

/* Euclidean 2-norm of a rank-1 buffer without spurious overflow/underflow.
 * The naive sqrt(sum x^2) overflows once any |x_i| exceeds sqrt(DBL_MAX)
 * (~1.3e154), and underflows every term to 0 for tiny inputs, even when the
 * true norm is representable -- so Norm[{1e200, 1e200}] returned inf and
 * Normalize of the same vector returned the zero vector.  This is LAPACK
 * dnrm2's running-scale recurrence: track the largest |component| as the scale
 * and accumulate the sum of squares relative to it.  For a complex buffer the
 * interleaved (re, im) reals are summed as two reals, which is exactly
 * sum |z_i|^2 = sum (re^2 + im^2). */
static double nd_vec_2norm(const double* buf, int n, bool cplx)
{
    int comps = cplx ? 2 * n : n;
    double scale = 0.0, ssq = 1.0;
    for (int i = 0; i < comps; i++) {
        double a = fabs(buf[i]);
        if (a != 0.0) {
            if (scale < a) {
                double r = scale / a;
                ssq = 1.0 + ssq * r * r;
                scale = a;
            } else {
                double r = a / scale;
                ssq += r * r;
            }
        }
    }
    return scale * sqrt(ssq);
}

Expr* ndla_norm(Expr* res)
{
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return linalg_delist_and_reeval(res);
    Expr* v = res->data.function.args[0];
    if (!is_ndarray(v)) return linalg_delist_and_reeval(res);

    if (v->data.ndarray.rank == 2) {                  /* matrix norm via LAPACK */
        Expr* r = ndla_matrix_norm_direct(v, argc == 2 ? res->data.function.args[1] : NULL);
        return r ? r : linalg_delist_and_reeval(res);
    }
    if (v->data.ndarray.rank != 1)
        return linalg_delist_and_reeval(res);

    /* Resolve p: default 2, positive Integer/Real, or Infinity. */
    double p = 2.0; bool p_inf = false;
    if (argc == 2) {
        Expr* pe = res->data.function.args[1];
        if (pe->type == EXPR_SYMBOL && pe->data.symbol.name == SYM_Infinity) p_inf = true;
        else if (pe->type == EXPR_INTEGER && pe->data.integer > 0) p = (double)pe->data.integer;
        else if (pe->type == EXPR_REAL && pe->data.real > 0.0) p = pe->data.real;
        else return linalg_delist_and_reeval(res);
    }

    int n; double* buf = NULL;
    bool cplx = nd_arg_is_complex(v);
    if (!na_load_vector(v, cplx, &n, &buf))
        return linalg_delist_and_reeval(res);

    double out;
    if (p_inf) {
        out = 0.0;
        for (int i = 0; i < n; i++) {
            double mag = cplx ? hypot(buf[2 * i], buf[2 * i + 1]) : fabs(buf[i]);
            if (mag > out) out = mag;
        }
    } else if (p == 2.0) {
        out = nd_vec_2norm(buf, n, cplx);
    } else {
        /* General p-norm, scaled by the max |component| so the running sum
         * neither overflows (large components) nor underflows to 0 (tiny
         * components) when the true norm is representable. */
        double mx = 0.0;
        for (int i = 0; i < n; i++) {
            double mag = cplx ? hypot(buf[2 * i], buf[2 * i + 1]) : fabs(buf[i]);
            if (mag > mx) mx = mag;
        }
        if (mx == 0.0) {
            out = 0.0;
        } else {
            double s = 0.0;
            for (int i = 0; i < n; i++) {
                double mag = cplx ? hypot(buf[2 * i], buf[2 * i + 1]) : fabs(buf[i]);
                s += pow(mag / mx, p);
            }
            out = mx * pow(s, 1.0 / p);
        }
    }
    free(buf);
    return expr_new_real(out);
}

/* ------------------------------------------------------------------ */
/*  Normalize (vector / its 2-norm)                                    */
/* ------------------------------------------------------------------ */

Expr* ndla_normalize(Expr* res)
{
    if (res->data.function.arg_count != 1) return linalg_delist_and_reeval(res);
    Expr* v = res->data.function.args[0];
    if (!is_ndarray(v) || v->data.ndarray.rank != 1)
        return linalg_delist_and_reeval(res);
    bool cplx = nd_arg_is_complex(v);
    int n; double* buf = NULL;
    if (!na_load_vector(v, cplx, &n, &buf))
        return linalg_delist_and_reeval(res);

    double nrm = nd_vec_2norm(buf, n, cplx);
    if (nrm > 0.0) {
        int comps = cplx ? 2 * n : n;
        for (int i = 0; i < comps; i++) buf[i] /= nrm;
    }
    Expr* out = nd_result_like(v, NULL, na_build_vector(buf, n, cplx));
    free(buf);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Cross (binary, length-3 vectors)                                   */
/* ------------------------------------------------------------------ */

Expr* ndla_cross(Expr* res)
{
    if (res->data.function.arg_count != 2) return linalg_delist_and_reeval(res);
    Expr* u = res->data.function.args[0];
    Expr* v = res->data.function.args[1];
    if (!is_ndarray(u) || !is_ndarray(v)
        || u->data.ndarray.rank != 1 || v->data.ndarray.rank != 1
        || u->data.ndarray.dims[0] != 3 || v->data.ndarray.dims[0] != 3)
        return linalg_delist_and_reeval(res);
    bool cplx = nd_arg_is_complex(u) || nd_arg_is_complex(v);
    int un, vn; double* U = NULL; double* V = NULL;
    if (!na_load_vector(u, cplx, &un, &U)) return linalg_delist_and_reeval(res);
    if (!na_load_vector(v, cplx, &vn, &V)) { free(U); return linalg_delist_and_reeval(res); }

    double out[6];
    if (!cplx) {
        out[0] = U[1] * V[2] - U[2] * V[1];
        out[1] = U[2] * V[0] - U[0] * V[2];
        out[2] = U[0] * V[1] - U[1] * V[0];
    } else {
        /* complex cross via interleaved (re,im); (a*b) is complex multiply */
        #define CR(i) U[2*(i)]
        #define CI(i) U[2*(i)+1]
        #define DR(i) V[2*(i)]
        #define DI(i) V[2*(i)+1]
        /* component k = U[a]*V[b] - U[b]*V[a] */
        int a[3] = {1, 2, 0}, bb[3] = {2, 0, 1};
        for (int k = 0; k < 3; k++) {
            int i = a[k], j = bb[k];
            double r1 = CR(i) * DR(j) - CI(i) * DI(j);
            double m1 = CR(i) * DI(j) + CI(i) * DR(j);
            double r2 = CR(j) * DR(i) - CI(j) * DI(i);
            double m2 = CR(j) * DI(i) + CI(j) * DR(i);
            out[2 * k]     = r1 - r2;
            out[2 * k + 1] = m1 - m2;
        }
        #undef CR
        #undef CI
        #undef DR
        #undef DI
    }
    free(U); free(V);
    return nd_result_like(u, v, na_build_vector(out, 3, cplx));
}

/* ------------------------------------------------------------------ */
/*  PseudoInverse / LeastSquares (real, from one gesdd)                */
/* ------------------------------------------------------------------ */

/*
 * Both heads were on pack.c's AWARE list already, so the transparency gate did
 * NOT materialise for them -- and both were still unusable on a machine matrix,
 * because each opened with linalg_delist_and_reeval and then ran the exact
 * pipeline in inv.c: rationalise every entry, row-reduce over Q, invert two
 * Gram matrices, numericalise. PseudoInverse[A300] and LeastSquares[A500,b500]
 * did not finish in 180 s. That is a hang under ordinary input, not a slow row.
 *
 * The shared kernel is one thin SVD. A = U S V^T with U m x k and V^T k x n
 * (k = min(m,n)); then
 *
 *     A^+   = V S^+ U^T                 (n x m)
 *     x     = V S^+ U^T b               (the minimum-norm least-squares x)
 *
 * where S^+ inverts the singular values above the cutoff and zeroes the rest.
 * LeastSquares does NOT form A^+ and multiply -- it applies the three factors
 * to b in order, which is O(k(m + n)) per right-hand side instead of O(mn) plus
 * an n x m product, and is the same value.
 *
 * Neither routine may defer by re-evaluating the call: builtin_pseudoinverse
 * dispatches here, so a linalg_delist_and_reeval decline would re-enter it.
 * NULL means decline, and the caller carries on with its own arguments.
 */

/* Singular-value cutoff. Automatic is LAPACK's own rank criterion, which is
 * also NumPy's default rcond for pinv/lstsq; an explicit Tolerance is a
 * fraction of the largest singular value, which is what Mathematica documents.
 * smax is S[0] -- gesdd returns the singular values in descending order. */
static double nd_svd_cutoff(int m, int n, double smax,
                            bool tol_automatic, double tol_value)
{
    if (!tol_automatic) return tol_value * smax;
    return (double)((m > n) ? m : n) * 2.220446049250313e-16 * smax;
}

/* A real rank-2 float64 buffer is the entire domain of both routines: the
 * exact pipeline is the RIGHT answer for an exact matrix, and an int64 buffer
 * cannot arrive because neither head is on pack.c's INT64_OK list. */
static bool nd_svd_arg_ok(const Expr* a)
{
    return a && is_ndarray(a) && a->data.ndarray.rank == 2
             && a->data.ndarray.dims[0] > 0 && a->data.ndarray.dims[1] > 0
             && !nd_arg_is_complex(a) && mathilda_lapack_probe();
}

/* Thin SVD of `a`. On success the caller owns *S, *U and *VT and must free all
 * three; U is column-major m x k, VT column-major k x n, both as gesdd leaves
 * them. Returns false on a load failure or a gesdd that did not converge. */
static bool nd_thin_svd(const Expr* a, int* m_out, int* n_out, int* k_out,
                        double** S_out, double** U_out, double** VT_out)
{
    int m, n; double* A = NULL;                       /* column-major m x n */
    if (!na_load_matrix(a, false, true, &m, &n, &A)) return false;
    int k = (m < n) ? m : n;

    double* S  = (double*)malloc(sizeof(double) * (size_t)k);
    double* U  = (double*)malloc(sizeof(double) * (size_t)m * (size_t)k);
    double* VT = (double*)malloc(sizeof(double) * (size_t)k * (size_t)n);
    if (!S || !U || !VT) {
        free(A); free(S); free(U); free(VT);
        return false;
    }
    int info = mat_lapack_dgesdd('S', m, n, A, m, S, U, m, VT, k);
    free(A);                                          /* gesdd destroys A */
    if (info != 0) { free(S); free(U); free(VT); return false; }

    *m_out = m; *n_out = n; *k_out = k;
    *S_out = S; *U_out = U; *VT_out = VT;
    return true;
}

Expr* ndla_pseudoinverse_direct(const Expr* a, bool tol_automatic, double tol_value)
{
    if (!nd_svd_arg_ok(a)) return NULL;

    int m, n, k; double *S = NULL, *U = NULL, *VT = NULL;
    if (!nd_thin_svd(a, &m, &n, &k, &S, &U, &VT)) return NULL;

    double cutoff = nd_svd_cutoff(m, n, S[0], tol_automatic, tol_value);

    /* Scale V^T's rows by 1/sigma in place, so the product below is a plain
     * V * (S^+ U^T). A singular value at or below the cutoff zeroes its row,
     * which is exactly the Moore-Penrose truncation. */
    for (int t = 0; t < k; t++) {
        double f = (S[t] > cutoff) ? 1.0 / S[t] : 0.0;
        for (int j = 0; j < n; j++) VT[t + (size_t)j * k] *= f;
    }

    /* P = (V^T)^T . U^T, i.e. n x m, column-major. */
    double* P = (double*)malloc(sizeof(double) * (size_t)n * (size_t)m);
    if (!P) { free(S); free(U); free(VT); return NULL; }
#ifdef USE_LAPACK
    /* Column-major: dgemm with CblasColMajor, C = A^T B^T for A = VT (k x n)
     * and B = U (m x k) gives P[i + j*n] = sum_t VT[t + i*k] * U[j + t*m]. */
    cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans,
                n, m, k, 1.0, VT, k, U, m, 0.0, P, n);
#else
    for (size_t z = 0; z < (size_t)n * (size_t)m; z++) P[z] = 0.0;
    for (int j = 0; j < m; j++)
        for (int t = 0; t < k; t++) {
            double c = U[j + (size_t)t * m];
            if (c == 0.0) continue;
            for (int i = 0; i < n; i++)
                P[i + (size_t)j * n] += VT[t + (size_t)i * k] * c;
        }
#endif
    free(S); free(U); free(VT);

    Expr* out = na_build_matrix(P, n, m, false, true);
    free(P);
    return nd_result_like(a, NULL, out);
}

Expr* ndla_leastsquares_direct(const Expr* a, const Expr* b,
                               bool tol_automatic, double tol_value)
{
    if (!nd_svd_arg_ok(a)) return NULL;
    int rows = (int)a->data.ndarray.dims[0];

    /* The right-hand side may be a vector or a matrix, packed or plain. A plain
     * List is accepted for the reason ndla_linearsolve accepts one: packing
     * keys on element count, so LeastSquares[A500, b90] arrives with the matrix
     * packed and the vector not, and deferring on that would send a machine
     * solve back to the pipeline this function exists to avoid. `a` being a
     * float64 buffer already makes the answer inexact, so loading an exact b
     * through doubles loses nothing the result would have kept. */
    int nrhs = 1; bool is_vec = true; int bn = 0, bc = 0; double* B = NULL;
    if (na_load_vector(b, false, &bn, &B)) {
        if (bn != rows) { free(B); return NULL; }
    } else if (na_load_matrix(b, false, true, &bn, &bc, &B)) {
        if (bn != rows) { free(B); return NULL; }
        nrhs = bc; is_vec = false;
    } else {
        return NULL;
    }

    int m, n, k; double *S = NULL, *U = NULL, *VT = NULL;
    if (!nd_thin_svd(a, &m, &n, &k, &S, &U, &VT)) { free(B); return NULL; }

    double cutoff = nd_svd_cutoff(m, n, S[0], tol_automatic, tol_value);

    /* X (n x nrhs, column-major) = V S^+ U^T B, applied factor by factor. */
    double* X = (double*)malloc(sizeof(double) * (size_t)n * (size_t)nrhs);
    double* t = (double*)malloc(sizeof(double) * (size_t)k);
    if (!X || !t) { free(X); free(t); free(B); free(S); free(U); free(VT); return NULL; }

    for (int j = 0; j < nrhs; j++) {
        const double* bj = B + (size_t)j * rows;
        for (int q = 0; q < k; q++) {                 /* t = S^+ U^T b */
            double s = 0.0;
            const double* uq = U + (size_t)q * m;     /* U column q */
            for (int i = 0; i < m; i++) s += uq[i] * bj[i];
            t[q] = (S[q] > cutoff) ? s / S[q] : 0.0;
        }
        double* xj = X + (size_t)j * n;
        for (int i = 0; i < n; i++) {                 /* x = V t */
            double s = 0.0;
            const double* vi = VT + (size_t)i * k;    /* V^T column i == V row i */
            for (int q = 0; q < k; q++) s += vi[q] * t[q];
            xj[i] = s;
        }
    }
    free(t); free(B); free(S); free(U); free(VT);

    Expr* out = is_vec ? na_build_vector(X, n, false)
                       : na_build_matrix(X, n, nrhs, false, true);
    free(X);
    return nd_result_like(a, is_ndarray(b) ? b : NULL, out);
}
