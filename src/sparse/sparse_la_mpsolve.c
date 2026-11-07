/* Mathilda -- sparse linear algebra: arbitrary-precision LinearSolve on a
 * SparseArray by mixed-precision iterative refinement.
 *
 * LinearSolve[s, b] with s a square SparseArray of arbitrary-precision reals
 * (N[s, 60]) used to be declined to the dense path: the gate densified s
 * (impossible beyond the densify cap -- SparseArray::dense, unevaluated) and
 * the dense MPFR LU then spent n^3/3 multi-limb operations on a matrix with a
 * handful of entries a row (n = 10000 at 60 digits: hours).
 *
 * The sparsity is used instead, the classical way (Wilkinson; Moler, "Iterative
 * refinement in floating point", JACM 14, 1967; for a solver that is not
 * exact in its own precision, Higham, Accuracy and Stability, 2nd ed., ch. 12):
 *
 *   1. A0 = s rounded to machine numbers, factored ONCE by the machine sparse
 *      backend (SpFactor: band / dense LU, UMFPACK / CHOLMOD where linked) --
 *      or, when no direct factor fits, the Krylov solver;
 *   2. x = 0; repeat
 *        r = b - s x           in MPFR at the working precision + 64 bits,
 *                              over the nonzeros only (O(nnz) multi-limb ops)
 *        A0 d = r / 2^e        in machine arithmetic (e scales r to O(1), so
 *                              r never underflows a double at any precision)
 *        x = x + 2^e d         in MPFR
 *      until |d| <= 2^-(p+4) |x|  (p the working precision in bits).
 *
 * Each pass multiplies the error by about cond(A) eps, so 60 digits take
 * four or five passes when cond(A) is moderate; every pass costs one machine
 * solve and one sparse MPFR product.
 *
 * The result is certified by the iteration itself: the corrections must
 * contract geometrically (|d_k| <= |d_{k-1}| / 2 from the second pass), and
 * then the error left after the last correction is at most that correction,
 * which is below 2^-(p+4) |x|.  When they do not contract -- A0 singular, or
 * cond(A) eps >= ~1/2, where machine corrections carry no information -- the
 * solve declines and the dense arbitrary-precision path answers as before.
 * The entries come back at the working precision p, as the dense MPFR LU
 * gives them.
 *
 * Scope: a real matrix (exact or MPFR entries, at least one inexact operand
 * above machine precision), a zero background, a vector right-hand side.
 * Complex systems, a matrix right-hand side and options go to the dense
 * path.  A TimeConstrained deadline is polled once a pass.
 */
#include "sparse_la.h"
#include "sparse_internal.h"
#include "sym_names.h"
#include "common.h"
#include "core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPFR
#include <mpfr.h>
#include "numeric.h"

/* Machine-real view: the pattern of `sa` with doubles (owned x / ci). */
static bool mp_machine_view(const SparseArray* sa, const mpfr_t* v, SpMat* A) {
    memset(A, 0, sizeof *A);
    const SparsePattern* p = sa->pat;
    A->m = sa->dims[0];
    A->n = sa->dims[1];
    A->nnz = p->nnz;
    A->rp = p->rowptr;
    size_t nz = (size_t)(p->nnz > 0 ? p->nnz : 1);
    if (p->idx32) {
        A->owned_ci = malloc(sizeof(int64_t) * nz);
        if (!A->owned_ci) return false;
        const int32_t* c32 = p->colidx;
        for (int64_t k = 0; k < p->nnz; k++) A->owned_ci[k] = c32[k];
        A->ci = A->owned_ci;
    } else {
        A->ci = p->colidx;
    }
    A->owned_x = malloc(sizeof(double) * nz);
    if (!A->owned_x) { free(A->owned_ci); memset(A, 0, sizeof *A); return false; }
    for (int64_t k = 0; k < p->nnz; k++) {
        double d = mpfr_get_d(v[k], MPFR_RNDN);
        if (!isfinite(d)) { spla_view_free(A); return false; }
        A->owned_x[k] = d;
    }
    A->x = A->owned_x;
    return true;
}

static mpfr_t* mp_vec_new(size_t n, mpfr_prec_t bits) {
    mpfr_t* v = malloc(sizeof(mpfr_t) * (n > 0 ? n : 1));
    if (!v) return NULL;
    for (size_t i = 0; i < n; i++) mpfr_init2(v[i], bits);
    return v;
}

static void mp_vec_free(mpfr_t* v, size_t n) {
    if (!v) return;
    for (size_t i = 0; i < n; i++) mpfr_clear(v[i]);
    free(v);
}

/* Real number leaf (exact or inexact) into `out`; false for anything else. */
static bool mp_read_real(const Expr* e, mpfr_t out, mpfr_t scratch_im) {
    bool inexact = false;
    if (!get_approx_mpfr(e, out, scratch_im, &inexact)) return false;
    return mpfr_zero_p(scratch_im) != 0 && mpfr_number_p(out);
}

/* The right-hand side as a vector of n reals, or NULL. */
static mpfr_t* mp_load_rhs(const Expr* b, int64_t n, mpfr_prec_t bits) {
    Expr* dense = NULL;
    const Expr* src = b;
    if (b->type == EXPR_SPARSE) {
        const SparseArray* sb = sp_of(b);
        if (sb->rank != 1) return NULL;
        dense = sp_normal(sb);
        if (!dense) return NULL;
        src = dense;
    }
    mpfr_t* v = NULL;
    if (src->type == EXPR_FUNCTION && src->data.function.head->type == EXPR_SYMBOL &&
        src->data.function.head->data.symbol.name == SYM_List &&
        (int64_t)src->data.function.arg_count == n) {
        v = mp_vec_new((size_t)n, bits);
        mpfr_t im;
        mpfr_init2(im, bits);
        for (int64_t i = 0; v && i < n; i++)
            if (!mp_read_real(src->data.function.args[i], v[i], im)) { mp_vec_free(v, (size_t)n); v = NULL; }
        mpfr_clear(im);
    }
    expr_free(dense);
    return v;
}

static void mp_norm_inf(const mpfr_t* v, int64_t n, mpfr_t out) {
    mpfr_set_zero(out, 1);
    for (int64_t i = 0; i < n; i++)
        if (mpfr_cmpabs(v[i], out) > 0) mpfr_abs(out, v[i], MPFR_RNDU);
}

Expr* spla_mp_solve(const Expr* Ae, const Expr* Be, long bits) {
    const SparseArray* sa = sp_of(Ae);
    if (!sa || sa->rank != 2 || sa->dims[0] != sa->dims[1] || sa->dims[0] < 2) return NULL;
    if (!spla_is_zero(sa->bg) || bits <= 53) return NULL;
    if (sa->kind != SPV_EXPR && sa->kind != SPV_INT64) return NULL;
    int64_t n = sa->dims[0];
    int64_t nnz = sa->pat->nnz;
    mpfr_prec_t p = (mpfr_prec_t)bits;
    mpfr_prec_t wp = p + 64;                      /* the residual's precision */

    /* the values at the working precision (exact entries converted) */
    mpfr_t* av = mp_vec_new((size_t)nnz, wp);
    if (!av) return NULL;
    {
        mpfr_t im;
        mpfr_init2(im, wp);
        bool ok = true;
        if (sa->kind == SPV_INT64) {
            const int64_t* iv = sa->vals;
            for (int64_t k = 0; k < nnz; k++) {
                /* exactly, in two 32-bit halves (no 64-bit long assumed) */
                int64_t hi = iv[k] / ((int64_t)1 << 32), lo = iv[k] - hi * ((int64_t)1 << 32);
                mpfr_set_d(av[k], (double)hi, MPFR_RNDN);
                mpfr_mul_2ui(av[k], av[k], 32, MPFR_RNDN);
                mpfr_add_d(av[k], av[k], (double)lo, MPFR_RNDN);
            }
        } else {
            Expr* const* ev = sa->vals;
            for (int64_t k = 0; k < nnz && ok; k++) ok = mp_read_real(ev[k], av[k], im);
        }
        mpfr_clear(im);
        if (!ok) { mp_vec_free(av, (size_t)nnz); return NULL; }
    }
    mpfr_t* bv = mp_load_rhs(Be, n, wp);
    if (!bv) { mp_vec_free(av, (size_t)nnz); return NULL; }

    SpMat A0;
    if (!mp_machine_view(sa, av, &A0)) { mp_vec_free(av, (size_t)nnz); mp_vec_free(bv, (size_t)n); return NULL; }
    SpfStatus st;
    SpFactor* f = spf_factor(&A0, n <= 128 ? SPF_WANT_DENSE : SPF_WANT_AUTO, &st);
    KryOpts ko;
    kry_defaults(&ko);
    bool use_krylov = !f && st == SPF_TOOBIG;

    Expr* result = NULL;
    mpfr_t* xv = mp_vec_new((size_t)n, wp);
    mpfr_t* rv = mp_vec_new((size_t)n, wp);
    double* rd = malloc(sizeof(double) * (size_t)n);
    double* dd = malloc(sizeof(double) * (size_t)n);
    mpfr_t rn, xn, dn, prev, t;
    mpfr_inits2(wp, rn, xn, dn, prev, t, (mpfr_ptr)0);
    const int64_t* ci64 = A0.ci;
    bool converged = false, timed_out = false;
    if ((f || use_krylov) && xv && rv && rd && dd) {
        for (int64_t i = 0; i < n; i++) mpfr_set_zero(xv[i], 1);
        mpfr_set_inf(prev, 1);
        int maxpass = (int)(p / 10) + 8;
        for (int pass = 0; pass < maxpass; pass++) {
            if (tc_deadline_passed()) { timed_out = true; break; }
            /* r = b - A x over the nonzeros */
            for (int64_t i = 0; i < n; i++) {
                mpfr_set_zero(rv[i], 1);
                for (int64_t k = A0.rp[i]; k < A0.rp[i + 1]; k++)
                    if (!mpfr_zero_p(xv[ci64[k]]))
                        mpfr_fma(rv[i], av[k], xv[ci64[k]], rv[i], MPFR_RNDN);
                mpfr_sub(rv[i], bv[i], rv[i], MPFR_RNDN);
            }
            mp_norm_inf(rv, n, rn);
            if (mpfr_zero_p(rn)) { converged = true; break; }
            /* scale r to O(1) and solve in machine arithmetic */
            long e = mpfr_get_exp(rn);
            for (int64_t i = 0; i < n; i++) {
                mpfr_mul_2si(t, rv[i], -e, MPFR_RNDN);
                rd[i] = mpfr_get_d(t, MPFR_RNDN);
            }
            bool ok;
            if (use_krylov) {
                memset(dd, 0, sizeof(double) * (size_t)n);
                int64_t it = 0;
                KryStatus ks = kry_solve(&A0, rd, dd, &ko, &it);
                ok = ks == KRY_CONVERGED || ks == KRY_MAXIT;
            } else {
                ok = spf_solve(f, rd, dd, 1, false);
            }
            if (!ok) break;
            /* x += 2^e d, and |d| for the contraction test */
            mpfr_set_zero(dn, 1);
            bool finite = true;
            for (int64_t i = 0; i < n; i++) {
                if (!isfinite(dd[i])) { finite = false; break; }
                mpfr_set_d(t, dd[i], MPFR_RNDN);
                mpfr_mul_2si(t, t, e, MPFR_RNDN);
                if (mpfr_cmpabs(t, dn) > 0) mpfr_abs(dn, t, MPFR_RNDU);
                mpfr_add(xv[i], xv[i], t, MPFR_RNDN);
            }
            if (!finite) break;
            /* contraction: the certificate of the refinement */
            if (pass > 0) {
                mpfr_div_2ui(t, prev, 1, MPFR_RNDN);
                if (mpfr_cmp(dn, t) > 0) break;
            }
            mpfr_set(prev, dn, MPFR_RNDN);
            mp_norm_inf(xv, n, xn);
            mpfr_mul_2si(t, xn, -(long)p - 4, MPFR_RNDN);
            if (mpfr_cmp(dn, t) <= 0) { converged = true; break; }
        }
    }
    if (converged) {
        Expr** items = malloc(sizeof(Expr*) * (size_t)n);
        if (items) {
            for (int64_t i = 0; i < n; i++) {
                Expr* x = expr_new_mpfr_bits(p);
                mpfr_set(x->data.mpfr, xv[i], MPFR_RNDN);
                items[i] = x;
            }
            result = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
            free(items);
        }
    }
    mpfr_clears(rn, xn, dn, prev, t, (mpfr_ptr)0);
    free(rd);
    free(dd);
    mp_vec_free(xv, xv ? (size_t)n : 0);
    mp_vec_free(rv, rv ? (size_t)n : 0);
    spf_free(f);
    spla_view_free(&A0);
    mp_vec_free(av, (size_t)nnz);
    mp_vec_free(bv, (size_t)n);
    if (timed_out) tc_check_deadline();           /* everything is freed */
    return result;
}

#else

Expr* spla_mp_solve(const Expr* Ae, const Expr* Be, long bits) {
    (void)Ae; (void)Be; (void)bits;
    return NULL;
}

#endif
