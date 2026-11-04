/*
 * src/linalg/hessenberg.c -- HessenbergDecomposition[m].
 *
 * Gives {p, h} with p unitary, h upper Hessenberg and
 * m == p . h . ConjugateTranspose[p].  Pivoting -> True also returns the
 * balancing matrix d (LAPACK gebal, job 'B'), the decomposition then being of
 * the balanced matrix: m . d == d . p . h . ConjugateTranspose[p].
 *
 * The reference treats this as a numerical decomposition: an exact matrix
 * gives ::schurf ("received a matrix with infinite precision") and a symbolic
 * one ::schurn, both unevaluated.  So:
 *
 *   machine  -> hess_machine: LAPACK gehrd + orghr (dgehrd/dorghr, or
 *               zgehrd/zunghr for a complex matrix).  Packed factors carrying
 *               the input's presentation.  The ONE kernel behind List, packed,
 *               NDArray and Compile[] (which lowers it as a rank-3 2 x n x n
 *               array {p, h}).
 *   MPFR     -> hess_mpfr: Householder reduction at the input precision with
 *               LAPACK's reflector convention (zlarfg: a real beta of sign
 *               opposite to Re x0), so its signs match the machine answer.
 *               Entries a Householder reduction never touches -- p's first row
 *               and column, h below the subdiagonal -- are exact 0, as in the
 *               reference.
 */

#include "la_dense.h"
#include "ndlinalg.h"
#include "numarray.h"
#include "lapack.h"
#include "eigen_internal.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "ndarray.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "message.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

#define HESS_HEAD "HessenbergDecomposition"

/* ------------------------------------------------------------------------ *
 *  Machine kernel                                                          *
 * ------------------------------------------------------------------------ */

/* Raw machine Hessenberg of an NDArray / boxed machine matrix of order n.
 * On success *P and *H are column-major buffers (complex: interleaved) and
 * *D the balancing matrix when `pivot`; returns false when the matrix is not
 * machine-numeric or LAPACK fails. */
static bool hess_machine_raw(const Expr* m, int n, bool pivot, bool* cplx_out,
                             double** P, double** H, double** D) {
    *P = *H = NULL; if (D) *D = NULL;
    int r = 0, c = 0;
    double* A = NULL;
    bool cplx = false;
    if (!na_load_matrix((Expr*)m, false, true, &r, &c, &A)) {
        free(A); A = NULL;
        if (!na_load_matrix((Expr*)m, true, true, &r, &c, &A)) { free(A); return false; }
        cplx = true;
    }
    if (r != n || c != n) { free(A); return false; }
    size_t comps = cplx ? 2 : 1, nn = (size_t)n * (size_t)n;
    int ilo = 1, ihi = n;
    if (pivot) {
        double* scale = (double*)malloc(sizeof(double) * (size_t)n);
        int info = cplx ? mat_lapack_zgebal('B', n, A, n, &ilo, &ihi, scale)
                        : mat_lapack_dgebal('B', n, A, n, &ilo, &ihi, scale);
        if (info != 0) { free(scale); free(A); return false; }
        double* Dm = (double*)calloc(comps * nn, sizeof(double));
        for (int i = 0; i < n; i++) Dm[comps * ((size_t)i + (size_t)i * n)] = 1.0;
        if (cplx) mat_lapack_zgebak('B', 'R', n, ilo, ihi, scale, n, Dm, n);
        else      mat_lapack_dgebak('B', 'R', n, ilo, ihi, scale, n, Dm, n);
        free(scale);
        *D = Dm;
    }
    double* tau = (double*)calloc(comps * (size_t)(n > 1 ? n - 1 : 1), sizeof(double));
    int info = 0;
    if (n > 1) info = cplx ? mat_lapack_zgehrd(n, ilo, ihi, A, n, tau)
                           : mat_lapack_dgehrd(n, ilo, ihi, A, n, tau);
    if (info != 0) { free(tau); free(A); if (D) { free(*D); *D = NULL; } return false; }
    /* h: A on and above the first subdiagonal */
    double* Hm = (double*)calloc(comps * nn, sizeof(double));
    for (int j = 0; j < n; j++)
        for (int i = 0; i <= j + 1 && i < n; i++) {
            size_t o = (size_t)i + (size_t)j * n;
            for (size_t q = 0; q < comps; q++) Hm[comps * o + q] = A[comps * o + q];
        }
    if (n > 1) info = cplx ? mat_lapack_zunghr(n, ilo, ihi, A, n, tau)
                           : mat_lapack_dorghr(n, ilo, ihi, A, n, tau);
    else { A[0] = 1.0; if (cplx) A[1] = 0.0; }
    free(tau);
    if (info != 0) { free(Hm); free(A); if (D) { free(*D); *D = NULL; } return false; }
    *P = A; *H = Hm; *cplx_out = cplx;
    return true;
}

/* {p, h} (or {p, h, d}) as packed factors, or NULL. */
static Expr* hess_machine(const Expr* m, int n, bool pivot) {
    double *P, *H, *D = NULL;
    bool cplx = false;
    if (!hess_machine_raw(m, n, pivot, &cplx, &P, &H, pivot ? &D : NULL)) return NULL;
    NDPresentation pres = na_result_presentation(m);
    Expr* items[3];
    items[0] = na_build_matrix_as(P, n, n, cplx, true, pres);
    items[1] = na_build_matrix_as(H, n, n, cplx, true, pres);
    size_t k = 2;
    if (pivot) items[k++] = na_build_matrix_as(D, n, n, cplx, true, pres);
    free(P); free(H); free(D);
    return expr_new_function(expr_new_symbol(SYM_List), items, k);
}

/* Compile[]: the {p, h} pair as ONE rank-3 2 x n x n machine array, from the
 * same raw kernel, or NULL (silently). */
Expr* la_hessenberg_machine(const Expr* m) {
    int r = 0, c = 0;
    if (!is_ndarray(m) || !la_matrix_dims(m, &r, &c) || r < 1 || r != c) return NULL;
    double *P, *H;
    bool cplx = false;
    if (!hess_machine_raw(m, r, false, &cplx, &P, &H, NULL)) return NULL;
    size_t comps = cplx ? 2 : 1, nn = (size_t)r * (size_t)r;
    double* buf = (double*)malloc(sizeof(double) * comps * 2 * nn);
    /* column-major -> row-major, p then h */
    for (int s = 0; s < 2; s++) {
        const double* src = s == 0 ? P : H;
        for (int i = 0; i < r; i++)
            for (int j = 0; j < r; j++)
                for (size_t q = 0; q < comps; q++)
                    buf[comps * ((size_t)s * nn + (size_t)i * r + j) + q] =
                        src[comps * ((size_t)i + (size_t)j * r) + q];
    }
    free(P); free(H);
    int64_t dims[3] = { 2, r, r };
    return expr_new_ndarray_like(m, 3, dims, buf, cplx ? NDT_COMPLEX64 : NDT_FLOAT64);
}

/* ------------------------------------------------------------------------ *
 *  MPFR kernel                                                             *
 * ------------------------------------------------------------------------ */

#ifdef USE_MPFR
/* complex helpers over separate re/im mpfr_t */
static Expr* hm_value(const mpfr_t re, const mpfr_t im, bool cplx) {
    if (!cplx) return expr_new_mpfr_copy(re);
    return la_ev2(SYM_Complex, expr_new_mpfr_copy(re), expr_new_mpfr_copy(im));
}

static Expr* hess_mpfr(const Expr* m, int n, mpfr_prec_t bits) {
    MatM A;
    if (!matM_load((Expr*)m, (size_t)n, bits, &A)) return NULL;
    bool cplx = A.is_complex != 0;
    size_t nn = (size_t)n * (size_t)n;
    mpfr_t* ar = A.re;
    mpfr_t* ai = cplx ? A.im : mpfr_array_alloc(nn, bits);
    if (!cplx) for (size_t k = 0; k < nn; k++) mpfr_set_zero(ai[k], 1);
    mpfr_t* qr = mpfr_array_alloc(nn, bits);
    mpfr_t* qi = mpfr_array_alloc(nn, bits);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            mpfr_set_si(qr[(size_t)i * n + j], (i == j) ? 1L : 0L, MPFR_RNDN);
            mpfr_set_zero(qi[(size_t)i * n + j], 1);
        }
    mpfr_t* vr = mpfr_array_alloc((size_t)n, bits);
    mpfr_t* vi = mpfr_array_alloc((size_t)n, bits);
    mpfr_t xn, beta, t, sr, si, taur, taui, dr, di, den, u;
    mpfr_inits2(bits, xn, beta, t, sr, si, taur, taui, dr, di, den, u, (mpfr_ptr)0);

#define AR(i, j) ar[(size_t)(i) * n + (j)]
#define AI(i, j) ai[(size_t)(i) * n + (j)]
    /* n - 1 reflectors, as zgehrd: the last one only makes a complex
     * subdiagonal entry real (it is the identity for a real matrix) */
    for (int k = 0; k + 1 < n; k++) {
        int p0 = k + 1;
        /* xnorm = || x[1:] ||, x = A[p0.., k] */
        mpfr_set_zero(xn, 1);
        for (int i = p0 + 1; i < n; i++) {
            mpfr_sqr(t, AR(i, k), MPFR_RNDN); mpfr_add(xn, xn, t, MPFR_RNDN);
            mpfr_sqr(t, AI(i, k), MPFR_RNDN); mpfr_add(xn, xn, t, MPFR_RNDN);
        }
        if (mpfr_zero_p(xn) && mpfr_zero_p(AI(p0, k))) continue;     /* H = I */
        /* beta = -sign(Re alpha) sqrt(|alpha|^2 + xnorm^2) */
        mpfr_sqr(t, AR(p0, k), MPFR_RNDN); mpfr_add(beta, xn, t, MPFR_RNDN);
        mpfr_sqr(t, AI(p0, k), MPFR_RNDN); mpfr_add(beta, beta, t, MPFR_RNDN);
        mpfr_sqrt(beta, beta, MPFR_RNDN);
        if (mpfr_sgn(AR(p0, k)) >= 0) mpfr_neg(beta, beta, MPFR_RNDN);
        /* tau = (beta - alpha) / beta */
        mpfr_sub(taur, beta, AR(p0, k), MPFR_RNDN); mpfr_div(taur, taur, beta, MPFR_RNDN);
        mpfr_neg(taui, AI(p0, k), MPFR_RNDN);       mpfr_div(taui, taui, beta, MPFR_RNDN);
        /* v = [1; x[1:] / (alpha - beta)] */
        mpfr_sub(dr, AR(p0, k), beta, MPFR_RNDN);
        mpfr_set(di, AI(p0, k), MPFR_RNDN);
        mpfr_sqr(den, dr, MPFR_RNDN); mpfr_sqr(t, di, MPFR_RNDN); mpfr_add(den, den, t, MPFR_RNDN);
        mpfr_set_ui(vr[p0], 1, MPFR_RNDN); mpfr_set_zero(vi[p0], 1);
        for (int i = p0 + 1; i < n; i++) {
            /* (xr + i xi)(dr - i di) / den */
            mpfr_mul(sr, AR(i, k), dr, MPFR_RNDN); mpfr_mul(t, AI(i, k), di, MPFR_RNDN);
            mpfr_add(sr, sr, t, MPFR_RNDN);        mpfr_div(vr[i], sr, den, MPFR_RNDN);
            mpfr_mul(si, AI(i, k), dr, MPFR_RNDN); mpfr_mul(t, AR(i, k), di, MPFR_RNDN);
            mpfr_sub(si, si, t, MPFR_RNDN);        mpfr_div(vi[i], si, den, MPFR_RNDN);
        }
        /* Left: A[p0.., c] -= conj(tau) v (v^H A[p0.., c]) for every column c */
        for (int c = 0; c < n; c++) {
            mpfr_set_zero(sr, 1); mpfr_set_zero(si, 1);   /* s = sum conj(v_s) A[s][c] */
            for (int s = p0; s < n; s++) {
                mpfr_mul(t, vr[s], AR(s, c), MPFR_RNDN); mpfr_add(sr, sr, t, MPFR_RNDN);
                mpfr_mul(t, vi[s], AI(s, c), MPFR_RNDN); mpfr_add(sr, sr, t, MPFR_RNDN);
                mpfr_mul(t, vr[s], AI(s, c), MPFR_RNDN); mpfr_add(si, si, t, MPFR_RNDN);
                mpfr_mul(t, vi[s], AR(s, c), MPFR_RNDN); mpfr_sub(si, si, t, MPFR_RNDN);
            }
            /* w = conj(tau) * s */
            mpfr_mul(dr, taur, sr, MPFR_RNDN); mpfr_mul(t, taui, si, MPFR_RNDN); mpfr_add(dr, dr, t, MPFR_RNDN);
            mpfr_mul(di, taur, si, MPFR_RNDN); mpfr_mul(t, taui, sr, MPFR_RNDN); mpfr_sub(di, di, t, MPFR_RNDN);
            for (int r = p0; r < n; r++) {   /* A[r][c] -= v_r * w */
                mpfr_mul(t, vr[r], dr, MPFR_RNDN); mpfr_sub(AR(r, c), AR(r, c), t, MPFR_RNDN);
                mpfr_mul(t, vi[r], di, MPFR_RNDN); mpfr_add(AR(r, c), AR(r, c), t, MPFR_RNDN);
                mpfr_mul(t, vr[r], di, MPFR_RNDN); mpfr_sub(AI(r, c), AI(r, c), t, MPFR_RNDN);
                mpfr_mul(t, vi[r], dr, MPFR_RNDN); mpfr_sub(AI(r, c), AI(r, c), t, MPFR_RNDN);
            }
        }
        /* Right on A and on Q: M[r, p0..] -= tau (M[r, p0..] v) v^H */
        for (int pass = 0; pass < 2; pass++) {
            mpfr_t* Mr = pass == 0 ? ar : qr;
            mpfr_t* Mi = pass == 0 ? ai : qi;
            for (int r = 0; r < n; r++) {
                mpfr_set_zero(sr, 1); mpfr_set_zero(si, 1);   /* s = sum M[r][s] v_s */
                for (int s = p0; s < n; s++) {
                    size_t o = (size_t)r * n + s;
                    mpfr_mul(t, Mr[o], vr[s], MPFR_RNDN); mpfr_add(sr, sr, t, MPFR_RNDN);
                    mpfr_mul(t, Mi[o], vi[s], MPFR_RNDN); mpfr_sub(sr, sr, t, MPFR_RNDN);
                    mpfr_mul(t, Mr[o], vi[s], MPFR_RNDN); mpfr_add(si, si, t, MPFR_RNDN);
                    mpfr_mul(t, Mi[o], vr[s], MPFR_RNDN); mpfr_add(si, si, t, MPFR_RNDN);
                }
                /* w = tau * s */
                mpfr_mul(dr, taur, sr, MPFR_RNDN); mpfr_mul(t, taui, si, MPFR_RNDN); mpfr_sub(dr, dr, t, MPFR_RNDN);
                mpfr_mul(di, taur, si, MPFR_RNDN); mpfr_mul(t, taui, sr, MPFR_RNDN); mpfr_add(di, di, t, MPFR_RNDN);
                for (int c = p0; c < n; c++) {   /* M[r][c] -= w conj(v_c) */
                    size_t o = (size_t)r * n + c;
                    mpfr_mul(t, dr, vr[c], MPFR_RNDN); mpfr_sub(Mr[o], Mr[o], t, MPFR_RNDN);
                    mpfr_mul(t, di, vi[c], MPFR_RNDN); mpfr_sub(Mr[o], Mr[o], t, MPFR_RNDN);
                    mpfr_mul(t, di, vr[c], MPFR_RNDN); mpfr_sub(Mi[o], Mi[o], t, MPFR_RNDN);
                    mpfr_mul(t, dr, vi[c], MPFR_RNDN); mpfr_add(Mi[o], Mi[o], t, MPFR_RNDN);
                }
            }
        }
        /* the reduced column: beta on the subdiagonal, exact zeros below */
        mpfr_set(AR(p0, k), beta, MPFR_RNDN); mpfr_set_zero(AI(p0, k), 1);
        for (int i = p0 + 1; i < n; i++) { mpfr_set_zero(AR(i, k), 1); mpfr_set_zero(AI(i, k), 1); }
    }
#undef AR
#undef AI

    Expr** pc = (Expr**)malloc(sizeof(Expr*) * nn);
    Expr** hc = (Expr**)malloc(sizeof(Expr*) * nn);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            size_t o = (size_t)i * n + j;
            bool pz = (i == 0 || j == 0) && i != j;      /* never touched */
            pc[o] = pz ? expr_new_integer(0) : hm_value(qr[o], qi[o], cplx);
            hc[o] = (i > j + 1) ? expr_new_integer(0) : hm_value(ar[o], ai[o], cplx);
        }
    Expr* items[2] = { la_matrix_from(pc, n, n), la_matrix_from(hc, n, n) };
    free(pc); free(hc);
    mpfr_clears(xn, beta, t, sr, si, taur, taui, dr, di, den, u, (mpfr_ptr)0);
    mpfr_array_free(vr, (size_t)n); mpfr_array_free(vi, (size_t)n);
    mpfr_array_free(qr, nn); mpfr_array_free(qi, nn);
    if (!cplx) mpfr_array_free(ai, nn);
    matM_free(&A);
    return expr_new_function(expr_new_symbol(SYM_List), items, 2);
}
#endif

/* ------------------------------------------------------------------------ *
 *  Entry points                                                            *
 * ------------------------------------------------------------------------ */

static Expr* hess_run(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return builtin_arg_error(HESS_HEAD, argc, 1, 1);
    bool pivot = false;
    for (size_t i = 1; i < argc; i++) {
        const Expr* o = res->data.function.args[i];
        bool rule = o->type == EXPR_FUNCTION && o->data.function.arg_count == 2
                 && (head_is(o, SYM_Rule) || head_is(o, SYM_RuleDelayed));
        if (!rule) return builtin_arg_error(HESS_HEAD, argc, 1, 1);
        const Expr* key = o->data.function.args[0];
        const Expr* rhs = o->data.function.args[1];
        const char* kn = key->type == EXPR_SYMBOL ? key->data.symbol.name : NULL;
        if (kn == SYM_Pivoting) {
            if (rhs->type != EXPR_SYMBOL
                || (rhs->data.symbol.name != SYM_True && rhs->data.symbol.name != SYM_False))
                return NULL;
            pivot = rhs->data.symbol.name == SYM_True;
        } else if (kn == SYM_TargetStructure) {
            continue;
        } else {
            char* k = la_str(key);
            char* call = la_str(res);
            mth_message(HESS_HEAD, "optx", "Unknown option %s in %s.", k, call);
            free(k); free(call);
            return NULL;
        }
    }
    const Expr* m = res->data.function.args[0];
    int r = 0, c = 0;
    if (!la_matrix_dims(m, &r, &c) || r < 1 || r != c) {
        char* s = la_str(m);
        mth_message(HESS_HEAD, "matsq", "Argument %s at position 1 is not a non-empty square matrix.", s);
        free(s);
        return NULL;
    }
    long bits = 53;
    LaPrecClass pc = la_prec_class(m, &bits);
    if (pc == LA_EXACT) {
        if (la_all_numeric(m))
            mth_message(HESS_HEAD, "schurf", "HessenbergDecomposition has received a matrix "
                        "with infinite precision.");
        else
            mth_message(HESS_HEAD, "schurn", "HessenbergDecomposition has received a matrix "
                        "with non-numeric elements.");
        return NULL;
    }
    if (!la_all_numeric(m)) {
        mth_message(HESS_HEAD, "schurn", "HessenbergDecomposition has received a matrix "
                    "with non-numeric elements.");
        return NULL;
    }
#ifdef USE_MPFR
    if (pc == LA_MPFR && !pivot) {
        Expr* boxed = la_boxed_copy(m);
        Expr* out = hess_mpfr(boxed, r, (mpfr_prec_t)bits);
        expr_free(boxed);
        if (out) return out;
    }
#endif
    Expr* out = hess_machine(m, r, pivot);
#ifdef USE_MPFR
    if (!out && !pivot) {
        /* no LAPACK: the Householder kernel at 53 bits, rounded to machine */
        Expr* boxed = la_boxed_copy(m);
        Expr* mp = hess_mpfr(boxed, r, 53);
        expr_free(boxed);
        if (mp) { out = common_numericalize_result(mp, 53); expr_free(mp); }
    }
#endif
    return out;
}

Expr* ndla_hessenbergdecomposition(Expr* res) { return hess_run(res); }

Expr* builtin_hessenbergdecomposition(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (linalg_call_has_ndarray(res)) return ndla_hessenbergdecomposition(res);
    return hess_run(res);
}

void hessenberg_init(void) {
    symtab_add_builtin(HESS_HEAD, builtin_hessenbergdecomposition);
    symtab_get_def(HESS_HEAD)->attributes |= ATTR_PROTECTED;
    Expr* r[2] = {
        la_call2(SYM_Rule, expr_new_symbol(SYM_Pivoting), expr_new_symbol(SYM_False)),
        la_call2(SYM_Rule, expr_new_symbol(SYM_TargetStructure), expr_new_symbol(SYM_Automatic)),
    };
    symtab_set_options(HESS_HEAD, expr_new_function(expr_new_symbol(SYM_List), r, 2));
}
