/*
 * schurdecomp_mpfr.c -- arbitrary-precision SchurDecomposition[m].
 *
 * Handles the standard, real, block-diagonal, un-pivoted case at MPFR
 * precision by reusing the in-house Hessenberg + Francis QR from eigen_direct.c
 * (eigen_schur_real_mpfr), which leaves the real Schur form T in the matrix and
 * the orthogonal Schur vectors Q alongside, with m == Q * T * Q^T.
 *
 * Complex, generalized, RealBlockDiagonalForm -> False, and Pivoting inputs are
 * NOT handled here (a documented gap -- they would need a complex MPFR QR / an
 * MPFR QZ); the dispatcher only calls this kernel for the supported case, and
 * on any failure it falls back to the machine LAPACK kernel.
 *
 * Memory contract: standard builtin ownership (SPEC.md §4).  Never frees the
 * input; every mpfr_t initialised here is cleared on every exit path.
 */

#include "schurdecomp.h"
#include "schurdecomp_internal.h"
#include "eigen.h"
#include "numeric.h"
#include "sym_names.h"
#include "expr.h"

#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

#ifndef USE_MPFR

/* Without MPFR there is no arbitrary-precision Schur kernel: schurdecomp.c only
 * references schur_mpfr_standard_real under #ifdef USE_MPFR, so nothing here is
 * declared or called.  A typedef keeps the translation unit non-empty. */
typedef int schurdecomp_mpfr_no_mpfr_build;

#else

/* Allocate / free a flat array of `count` mpfr_t cells at `bits` precision. */
static mpfr_t* schur_mpfr_alloc(size_t count, mpfr_prec_t bits) {
    mpfr_t* a = (mpfr_t*)malloc(sizeof(mpfr_t) * count);
    for (size_t i = 0; i < count; i++) mpfr_init2(a[i], bits);
    return a;
}

static void schur_mpfr_free(mpfr_t* a, size_t count) {
    for (size_t i = 0; i < count; i++) mpfr_clear(a[i]);
    free(a);
}

/* Load the n x n real matrix (a boxed List-of-Lists -- an NDArray never reaches
 * this arbitrary-precision path) into a row-major mpfr buffer at `bits`.  The
 * imaginary part is discarded (the caller asserts realness via schur_all_real).
 * Returns false, leaking nothing, on a ragged shape or a non-numeric leaf. */
static bool schur_mpfr_load(const Expr* m, int n, mpfr_prec_t bits, mpfr_t* A) {
    if (m->type != EXPR_FUNCTION) return false;
    mpfr_t tmp_im;
    mpfr_init2(tmp_im, bits);
    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (row->type != EXPR_FUNCTION || (int)row->data.function.arg_count != n) {
            mpfr_clear(tmp_im);
            return false;
        }
        for (int j = 0; j < n; j++) {
            bool inexact = false;
            if (!get_approx_mpfr(row->data.function.args[j],
                                 A[(size_t)i * (size_t)n + (size_t)j],
                                 tmp_im, &inexact)) {
                mpfr_clear(tmp_im);
                return false;
            }
        }
    }
    mpfr_clear(tmp_im);
    return true;
}

/* Wrap a row-major n x n mpfr buffer as a nested List of EXPR_MPFR entries. */
static Expr* schur_mpfr_wrap(const mpfr_t* buf, int n) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int j = 0; j < n; j++)
            elems[j] = expr_new_mpfr_copy(buf[(size_t)i * (size_t)n + (size_t)j]);
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), elems, (size_t)n);
        free(elems);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    return out;
}

Expr* schur_mpfr_standard_real(const Expr* m, int n, mpfr_prec_t bits,
                               const SchurOpts* opts) {
    (void)opts;  /* only reached for real / block-diagonal / un-pivoted input */

    size_t nn = (size_t)n * (size_t)n;
    mpfr_t* A = schur_mpfr_alloc(nn, bits);  /* -> Schur form T */
    if (!schur_mpfr_load(m, n, bits, A)) {
        schur_mpfr_free(A, nn);
        return NULL;
    }

    mpfr_t* Q = schur_mpfr_alloc(nn, bits);  /* -> Schur vectors q */
    int status = eigen_schur_real_mpfr(A, n, bits, Q);

    Expr* result = NULL;
    if (status == 0) {
        Expr* q = schur_mpfr_wrap(Q, n);
        Expr* t = schur_mpfr_wrap(A, n);
        Expr* items[2] = { q, t };
        result = expr_new_function(expr_new_symbol(SYM_List), items, 2);
    }

    schur_mpfr_free(A, nn);
    schur_mpfr_free(Q, nn);
    return result;
}

#endif /* USE_MPFR */
