/* Mathilda — NDSolve arbitrary-precision integrator.
 *
 * A genuine MPFR path: the state vector Y, the independent variable t and the
 * step size h are all carried as mpfr_t at a guard-padded working precision, so
 * non-autonomous right-hand sides f(t, Y) are evaluated at full precision.  The
 * adaptive step-size control (a scalar error ratio) is done in double, which is
 * ample for step selection.  Reuses the same tableaux as the machine path
 * (ndsolve_tableau.h).  Only the explicit methods are provided here; a
 * high-precision request for an implicit/multistep method uses the adaptive
 * DOPRI5 integrator (a note is emitted by the front-end). */
#include "ndsolve_common.h"
#include "ndsolve_tableau.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../common.h"
#include "../arithmetic.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>

#ifdef USE_MPFR

/* ---- mpfr vector helpers ---- */
static mpfr_t* mp_vec(size_t n, long bits) {
    mpfr_t* v = malloc(sizeof(mpfr_t) * n);
    for (size_t i = 0; i < n; i++) mpfr_init2(v[i], bits);
    return v;
}
static void mp_vec_free(mpfr_t* v, size_t n) {
    if (!v) return;
    for (size_t i = 0; i < n; i++) mpfr_clear(v[i]);
    free(v);
}

/* Evaluate the reduced RHS at (t, Y) in MPFR: out[i] = f_i(t, Y). */
static bool nd_rhs_mpfr(NdProblem* P, const mpfr_t t, mpfr_t* Y, mpfr_t* out, long bits) {
    size_t d = P->d;
    nd_bind_set(&P->bind_t, expr_new_mpfr_copy(t));
    for (size_t i = 0; i < d; i++) nd_bind_set(&P->bind_y[i], expr_new_mpfr_copy(Y[i]));
    eval_clock_bump();
    if (P->eval_monitor) { Expr* m = eval_and_free(expr_copy(P->eval_monitor)); expr_free(m); }
    mpfr_t im; mpfr_init2(im, bits);
    bool ok = true;
    for (size_t i = 0; i < d && ok; i++) {
        arith_warnings_mute_push();
        Expr* raw = eval_and_free(expr_copy(P->f[i]));
        arith_warnings_mute_pop();
        if (!raw) { ok = false; break; }
        Expr* num = numericalize(raw, P->spec);
        expr_free(raw);
        bool inexact;
        if (!num || !get_approx_mpfr(num, out[i], im, &inexact) || !mpfr_number_p(out[i]))
            ok = false;
        expr_free(num);
    }
    mpfr_clear(im);
    return ok;
}

/* Scalar WRMS error from mpfr e against a double-scale sc_i = atol+rtol*max|y|. */
static double mp_wrms(size_t d, mpfr_t* e, mpfr_t* y, mpfr_t* yn, NdTol tol) {
    double sum = 0.0;
    for (size_t i = 0; i < d; i++) {
        double ay = fabs(mpfr_get_d(y[i], MPFR_RNDN));
        double an = fabs(mpfr_get_d(yn[i], MPFR_RNDN));
        double sc = tol.atol + tol.rtol * (ay > an ? ay : an);
        if (sc <= 0.0) sc = 1e-300;
        double r = mpfr_get_d(e[i], MPFR_RNDN) / sc;
        sum += r * r;
    }
    return sqrt(sum / (double)d);
}

/* One DOPRI5 step in MPFR from (t, Y) by h.  Writes Ynew, and (adaptive) the
 * error estimate into Yerr.  k must hold 7*d mpfr slots. */
static bool dopri5_mpfr(NdProblem* P, const mpfr_t t, mpfr_t* Y, const mpfr_t h,
                        mpfr_t* Ynew, mpfr_t* Yerr, mpfr_t* k, long bits) {
    size_t d = P->d;
    mpfr_t* tmp = mp_vec(d, bits);
    mpfr_t ts, acc, prod;
    mpfr_init2(ts, bits); mpfr_init2(acc, bits); mpfr_init2(prod, bits);
    bool ok = nd_rhs_mpfr(P, t, Y, &k[0], bits);
    for (int s = 1; s < 7 && ok; s++) {
        for (size_t i = 0; i < d; i++) {
            mpfr_set_zero(acc, 1);
            for (int j = 0; j < s; j++) {
                mpfr_mul_d(prod, k[j*d + i], DP_A[s][j], MPFR_RNDN);
                mpfr_add(acc, acc, prod, MPFR_RNDN);
            }
            mpfr_mul(acc, acc, h, MPFR_RNDN);       /* h * sum a_sj k_j */
            mpfr_add(tmp[i], Y[i], acc, MPFR_RNDN);
        }
        mpfr_mul_d(ts, h, DP_C[s], MPFR_RNDN);
        mpfr_add(ts, ts, t, MPFR_RNDN);
        ok = nd_rhs_mpfr(P, ts, tmp, &k[s*d], bits);
    }
    if (ok) {
        for (size_t i = 0; i < d; i++) {
            mpfr_t sol, err; mpfr_init2(sol, bits); mpfr_init2(err, bits);
            mpfr_set_zero(sol, 1); mpfr_set_zero(err, 1);
            for (int s = 0; s < 7; s++) {
                mpfr_mul_d(prod, k[s*d + i], DP_B[s], MPFR_RNDN); mpfr_add(sol, sol, prod, MPFR_RNDN);
                mpfr_mul_d(prod, k[s*d + i], DP_E[s], MPFR_RNDN); mpfr_add(err, err, prod, MPFR_RNDN);
            }
            mpfr_mul(sol, sol, h, MPFR_RNDN); mpfr_add(Ynew[i], Y[i], sol, MPFR_RNDN);
            if (Yerr) mpfr_mul(Yerr[i], err, h, MPFR_RNDN);
            mpfr_clear(sol); mpfr_clear(err);
        }
    }
    mp_vec_free(tmp, d);
    mpfr_clear(ts); mpfr_clear(acc); mpfr_clear(prod);
    return ok;
}

/* One classical RK4 step in MPFR (fixed order 4). */
static bool rk4_mpfr(NdProblem* P, const mpfr_t t, mpfr_t* Y, const mpfr_t h,
                     mpfr_t* Ynew, mpfr_t* k, long bits) {
    size_t d = P->d;
    mpfr_t* tmp = mp_vec(d, bits);
    mpfr_t th, hh, prod; mpfr_init2(th, bits); mpfr_init2(hh, bits); mpfr_init2(prod, bits);
    mpfr_mul_d(hh, h, 0.5, MPFR_RNDN);
    bool ok = nd_rhs_mpfr(P, t, Y, &k[0], bits);
    if (ok) { for (size_t i = 0; i < d; i++) { mpfr_mul(prod, hh, k[0*d+i], MPFR_RNDN); mpfr_add(tmp[i], Y[i], prod, MPFR_RNDN); }
              mpfr_add(th, t, hh, MPFR_RNDN); ok = nd_rhs_mpfr(P, th, tmp, &k[1*d], bits); }
    if (ok) { for (size_t i = 0; i < d; i++) { mpfr_mul(prod, hh, k[1*d+i], MPFR_RNDN); mpfr_add(tmp[i], Y[i], prod, MPFR_RNDN); }
              ok = nd_rhs_mpfr(P, th, tmp, &k[2*d], bits); }
    if (ok) { for (size_t i = 0; i < d; i++) { mpfr_mul(prod, h, k[2*d+i], MPFR_RNDN); mpfr_add(tmp[i], Y[i], prod, MPFR_RNDN); }
              mpfr_add(th, t, h, MPFR_RNDN); ok = nd_rhs_mpfr(P, th, tmp, &k[3*d], bits); }
    if (ok) for (size_t i = 0; i < d; i++) {
        mpfr_t s; mpfr_init2(s, bits); mpfr_set_zero(s, 1);
        mpfr_add(s, s, k[0*d+i], MPFR_RNDN);
        mpfr_mul_d(prod, k[1*d+i], 2.0, MPFR_RNDN); mpfr_add(s, s, prod, MPFR_RNDN);
        mpfr_mul_d(prod, k[2*d+i], 2.0, MPFR_RNDN); mpfr_add(s, s, prod, MPFR_RNDN);
        mpfr_add(s, s, k[3*d+i], MPFR_RNDN);
        mpfr_mul(s, s, h, MPFR_RNDN); mpfr_div_d(s, s, 6.0, MPFR_RNDN);
        mpfr_add(Ynew[i], Y[i], s, MPFR_RNDN); mpfr_clear(s);
    }
    mp_vec_free(tmp, d);
    mpfr_clear(th); mpfr_clear(hh); mpfr_clear(prod);
    return ok;
}

/* ---- mpfr solution store ---- */
typedef struct { size_t n, cap, d; long bits; mpfr_t* t; mpfr_t* Y; mpfr_t* dY; } MpSol;
static void mpsol_init(MpSol* s, size_t d, long bits) { s->n=s->cap=0; s->d=d; s->bits=bits; s->t=NULL; s->Y=NULL; s->dY=NULL; }
static void mpsol_free(MpSol* s) {
    for (size_t i = 0; i < s->n; i++) mpfr_clear(s->t[i]);
    for (size_t i = 0; i < s->n * s->d; i++) { mpfr_clear(s->Y[i]); mpfr_clear(s->dY[i]); }
    free(s->t); free(s->Y); free(s->dY);
}
static void mpsol_push(MpSol* s, const mpfr_t t, mpfr_t* Y, mpfr_t* dY) {
    if (s->n == s->cap) {
        size_t nc = s->cap ? s->cap*2 : 64;
        s->t  = realloc(s->t,  sizeof(mpfr_t) * nc);
        s->Y  = realloc(s->Y,  sizeof(mpfr_t) * nc * s->d);
        s->dY = realloc(s->dY, sizeof(mpfr_t) * nc * s->d);
        s->cap = nc;
    }
    mpfr_init2(s->t[s->n], s->bits); mpfr_set(s->t[s->n], t, MPFR_RNDN);
    for (size_t i = 0; i < s->d; i++) {
        mpfr_init2(s->Y[s->n*s->d + i],  s->bits); mpfr_set(s->Y[s->n*s->d + i],  Y[i],  MPFR_RNDN);
        mpfr_init2(s->dY[s->n*s->d + i], s->bits); mpfr_set(s->dY[s->n*s->d + i], dY[i], MPFR_RNDN);
    }
    s->n++;
}
/* insertion sort by ascending t + drop duplicates */
static void mpsol_sort(MpSol* s) {
    size_t d = s->d;
    for (size_t i = 1; i < s->n; i++) {
        for (size_t j = i; j > 0 && mpfr_cmp(s->t[j-1], s->t[j]) > 0; j--) {
            mpfr_swap(s->t[j-1], s->t[j]);
            for (size_t c = 0; c < d; c++) { mpfr_swap(s->Y[(j-1)*d+c], s->Y[j*d+c]); mpfr_swap(s->dY[(j-1)*d+c], s->dY[j*d+c]); }
        }
    }
}

/* ================================================================= *
 *  MPFR implicit (stiff) BDF — variable step, variable order 1-5     *
 * ================================================================= *
 *
 * The state, node times and BDF coefficients are carried in MPFR (accuracy-
 * critical, since the coefficients scale the state into the theta-method base);
 * the predictor, error estimate, order and step control are done in double
 * (heuristics).  The Newton residual and update are MPFR, but the iteration
 * matrix I - coef*J reuses the DOUBLE Jacobian (nd_jacobian_real) and double
 * dense solve: an inexact Jacobian only changes Newton's convergence rate, not
 * the root, which the MPFR residual drives to full working precision.  This is
 * what lets stiff problems (where the explicit MPFR integrator would need an
 * impractically tiny step) run at arbitrary precision. */
#define MP_QMAX 5

static double mp_grow_cap(int q) {
    switch (q) { case 1: case 2: return 5.0; case 3: return 4.0; case 4: return 3.0; default: return 2.0; }
}

/* BDF-q coefficients c[0..q] in MPFR from node times: x0=t1 (unknown),
 * x_{1..q}=th[0..q-1].  c_j = L_j'(t1). */
static void mp_bdf_coeffs(const mpfr_t t1, mpfr_t* th, int q, mpfr_t* c, long bits) {
    mpfr_t x[MP_QMAX + 1], diff, num, den;
    for (int j = 0; j <= q; j++) mpfr_init2(x[j], bits);
    mpfr_init2(diff, bits); mpfr_init2(num, bits); mpfr_init2(den, bits);
    mpfr_set(x[0], t1, MPFR_RNDN);
    for (int j = 1; j <= q; j++) mpfr_set(x[j], th[j - 1], MPFR_RNDN);
    mpfr_set_zero(c[0], 1);
    for (int m = 1; m <= q; m++) { mpfr_sub(diff, x[0], x[m], MPFR_RNDN);
        mpfr_ui_div(diff, 1, diff, MPFR_RNDN); mpfr_add(c[0], c[0], diff, MPFR_RNDN); }
    for (int j = 1; j <= q; j++) {
        mpfr_set_ui(num, 1, MPFR_RNDN); mpfr_set_ui(den, 1, MPFR_RNDN);
        for (int m = 0; m <= q; m++) {
            if (m == j) continue;
            mpfr_sub(diff, x[j], x[m], MPFR_RNDN); mpfr_mul(den, den, diff, MPFR_RNDN);
            if (m == 0) continue;
            mpfr_sub(diff, x[0], x[m], MPFR_RNDN); mpfr_mul(num, num, diff, MPFR_RNDN);
        }
        mpfr_div(c[j], num, den, MPFR_RNDN);
    }
    for (int j = 0; j <= q; j++) mpfr_clear(x[j]);
    mpfr_clear(diff); mpfr_clear(num); mpfr_clear(den);
}

/* Degree-q polynomial predictor in MPFR through (th[0..q], yh rows), evaluated at
 * t1.  MPFR is required: for tight tolerances the corrector-predictor difference
 * (the local-error estimate) is far below double roundoff of the O(1) states, so
 * a double predictor would drown the estimate in cancellation noise. */
static void mp_predict(const mpfr_t t1, mpfr_t* th, mpfr_t* yh, int q, size_t d,
                       mpfr_t* pred, long bits) {
    mpfr_t w[MP_QMAX + 1], nu, de, diff, term;
    for (int j = 0; j <= q; j++) mpfr_init2(w[j], bits);
    mpfr_init2(nu, bits); mpfr_init2(de, bits); mpfr_init2(diff, bits); mpfr_init2(term, bits);
    for (int j = 0; j <= q; j++) {
        mpfr_set_ui(nu, 1, MPFR_RNDN); mpfr_set_ui(de, 1, MPFR_RNDN);
        for (int mm = 0; mm <= q; mm++) {
            if (mm == j) continue;
            mpfr_sub(diff, t1, th[mm], MPFR_RNDN); mpfr_mul(nu, nu, diff, MPFR_RNDN);
            mpfr_sub(diff, th[j], th[mm], MPFR_RNDN); mpfr_mul(de, de, diff, MPFR_RNDN);
        }
        mpfr_div(w[j], nu, de, MPFR_RNDN);
    }
    for (size_t i = 0; i < d; i++) {
        mpfr_set_zero(pred[i], 1);
        for (int j = 0; j <= q; j++) { mpfr_mul(term, w[j], yh[(size_t)j * d + i], MPFR_RNDN);
            mpfr_add(pred[i], pred[i], term, MPFR_RNDN); }
    }
    for (int j = 0; j <= q; j++) mpfr_clear(w[j]);
    mpfr_clear(nu); mpfr_clear(de); mpfr_clear(diff); mpfr_clear(term);
}

/* MPFR dense Gaussian elimination with partial pivoting: solve M x = b in place
 * (M destroyed, b overwritten with x).  Returns false on a singular pivot. */
static bool mp_dense_solve(size_t n, mpfr_t* M, mpfr_t* b, long bits) {
    mpfr_t factor, t; mpfr_init2(factor, bits); mpfr_init2(t, bits);
    bool ok = true;
    for (size_t col = 0; col < n && ok; col++) {
        size_t piv = col;
        for (size_t r = col + 1; r < n; r++)
            if (mpfr_cmpabs(M[r*n + col], M[piv*n + col]) > 0) piv = r;
        if (mpfr_zero_p(M[piv*n + col])) { ok = false; break; }
        if (piv != col) {
            for (size_t j = col; j < n; j++) mpfr_swap(M[col*n + j], M[piv*n + j]);
            mpfr_swap(b[col], b[piv]);
        }
        for (size_t r = col + 1; r < n; r++) {
            mpfr_div(factor, M[r*n + col], M[col*n + col], MPFR_RNDN);
            if (mpfr_zero_p(factor)) continue;
            for (size_t j = col; j < n; j++) {
                mpfr_mul(t, factor, M[col*n + j], MPFR_RNDN);
                mpfr_sub(M[r*n + j], M[r*n + j], t, MPFR_RNDN);
            }
            mpfr_mul(t, factor, b[col], MPFR_RNDN);
            mpfr_sub(b[r], b[r], t, MPFR_RNDN);
        }
    }
    for (size_t i = n; ok && i-- > 0; ) {
        for (size_t j = i + 1; j < n; j++) {
            mpfr_mul(t, M[i*n + j], b[j], MPFR_RNDN);
            mpfr_sub(b[i], b[i], t, MPFR_RNDN);
        }
        mpfr_div(b[i], b[i], M[i*n + i], MPFR_RNDN);
    }
    mpfr_clear(factor); mpfr_clear(t);
    return ok;
}

/* MPFR Newton for  Z = Ybase + coef * f(t1, Z)  (coef = 1/c0).  The iteration
 * matrix I - coef*J is built from the DOUBLE Jacobian (its inaccuracy only slows
 * convergence, not the root), but the residual G and the linear SOLVE are MPFR,
 * so the correction dZ — and hence the root — reaches full working precision
 * (a double solve would floor dZ at ~1e-16 and never satisfy a tight goal). */
static bool mp_newton_bdf(NdProblem* P, const mpfr_t t1, mpfr_t* Ybase, const mpfr_t coef,
                          mpfr_t* guess, mpfr_t* Ynew, NdTol tol, long bits) {
    size_t d = P->d;
    mpfr_t* Z = mp_vec(d, bits);
    mpfr_t* f = mp_vec(d, bits);
    mpfr_t* G = mp_vec(d, bits);
    mpfr_t* M = mp_vec(d * d, bits);
    mpfr_t  tmp; mpfr_init2(tmp, bits);
    double* Jd = malloc(sizeof(double) * d * d);
    double* Zd = malloc(sizeof(double) * d);
    double* dZd = malloc(sizeof(double) * d);
    for (size_t i = 0; i < d; i++) mpfr_set(Z[i], guess[i], MPFR_RNDN);
    double t1d = mpfr_get_d(t1, MPFR_RNDN), coefd = mpfr_get_d(coef, MPFR_RNDN);
    bool conv = false;
    for (int it = 0; it < 20; it++) {
        if (!nd_rhs_mpfr(P, t1, Z, f, bits)) break;
        for (size_t i = 0; i < d; i++) {                 /* G = Z - Ybase - coef f  */
            mpfr_mul(tmp, coef, f[i], MPFR_RNDN);
            mpfr_sub(tmp, Z[i], tmp, MPFR_RNDN);
            mpfr_sub(G[i], tmp, Ybase[i], MPFR_RNDN);
            Zd[i] = mpfr_get_d(Z[i], MPFR_RNDN);
        }
        if (!nd_jacobian_real(P, t1d, Zd, Jd)) break;    /* double Jacobian        */
        for (size_t i = 0; i < d; i++)
            for (size_t j = 0; j < d; j++)
                mpfr_set_d(M[i*d + j], (i == j ? 1.0 : 0.0) - coefd * Jd[i*d + j], MPFR_RNDN);
        if (!mp_dense_solve(d, M, G, bits)) break;        /* G := dZ (MPFR)         */
        for (size_t i = 0; i < d; i++) {
            mpfr_sub(Z[i], Z[i], G[i], MPFR_RNDN);
            dZd[i] = mpfr_get_d(G[i], MPFR_RNDN);
        }
        double nrm = nd_wrms_norm(d, dZd, Zd, NULL, tol);
        if (nrm <= 1e-3) { conv = true; break; }
    }
    if (conv) for (size_t i = 0; i < d; i++) mpfr_set(Ynew[i], Z[i], MPFR_RNDN);
    mp_vec_free(Z, d); mp_vec_free(f, d); mp_vec_free(G, d); mp_vec_free(M, d * d);
    mpfr_clear(tmp); free(Jd); free(Zd); free(dZd);
    return conv;
}

/* Variable-step variable-order MPFR BDF in one direction (mirrors the machine
 * bdf_dir: MPFR state + double control). */
static NdStatus mpfr_bdf_dir(NdProblem* P, const NdOpts* o, MpSol* sol, NdTol tol,
                             double target, long bits, int64_t* budget) {
    size_t d = P->d;
    double dir = (target > P->t0) ? 1.0 : -1.0;
    double span = fabs(P->tmax - P->tmin);
    if (fabs(target - P->t0) <= 1e-15 * (span + 1.0)) return ND_OK;
    const int CAP = MP_QMAX + 1;

    mpfr_t* yh   = mp_vec((size_t)CAP * d, bits);   /* history states (row 0 = newest) */
    mpfr_t* th   = mp_vec((size_t)CAP, bits);       /* history node times              */
    mpfr_t* cc   = mp_vec(MP_QMAX + 1, bits);
    mpfr_t* base = mp_vec(d, bits);
    mpfr_t* ynext= mp_vec(d, bits);
    mpfr_t* fnext= mp_vec(d, bits);
    mpfr_t* predm= mp_vec(d, bits);            /* order-q predictor (Newton guess) */
    mpfr_t* predm2=mp_vec(d, bits);            /* order-(q-1) predictor (down-check)*/
    mpfr_t* diffm= mp_vec(d, bits);            /* corrector - predictor (MPFR)     */
    mpfr_t  t, t1m, coef, acc, tmp;
    mpfr_init2(t, bits); mpfr_init2(t1m, bits); mpfr_init2(coef, bits);
    mpfr_init2(acc, bits); mpfr_init2(tmp, bits);
    double* y0d    = malloc(sizeof(double) * d);
    double* fcur_d = malloc(sizeof(double) * d);

    mpfr_set_d(t, P->t0, MPFR_RNDN);
    mpfr_set(th[0], t, MPFR_RNDN);
    for (size_t i = 0; i < d; i++) { mpfr_set_d(yh[i], P->Y0[i], MPFR_RNDN); y0d[i] = P->Y0[i]; }
    int m = 1;
    NdStatus st = ND_OK;
    if (!nd_rhs_real(P, P->t0, y0d, fcur_d)) { st = ND_ERR_SAMPLE; goto done; }

    double h = (o->starting_step > 0.0) ? dir * o->starting_step
                                        : nd_initial_step(P, o, tol, P->t0, y0d, fcur_d, 1, dir);
    if (h == 0.0) h = dir * span / 100.0;
    int qtarget = 1, hold = 0;
    bool no_grow = false;
    double h_cap = (o->max_step_size > 0.0) ? o->max_step_size : HUGE_VAL;
    double frac_cap = (o->max_step_fraction > 0.0) ? o->max_step_fraction * span : HUGE_VAL;

    while (dir * (target - mpfr_get_d(t, MPFR_RNDN)) > 1e-14 * (fabs(mpfr_get_d(t, MPFR_RNDN)) + 1.0)) {
        double td = mpfr_get_d(t, MPFR_RNDN);
        double hmag = fabs(h);
        if (hmag > h_cap) hmag = h_cap;
        if (hmag > frac_cap) hmag = frac_cap;
        double hs = dir * hmag;
        if ((td + hs - target) * dir > 0.0) hs = target - td;
        if (fabs(hs) < 1e-14 * (fabs(td) + 1.0)) { st = ND_ERR_STEPSIZE; break; }
        if (--(*budget) < 0) { st = ND_ERR_MAXSTEPS; break; }
        mpfr_add_d(t1m, t, hs, MPFR_RNDN);

        int q;
        bool ok;
        if (m == 1) {
            q = 1;                                           /* backward Euler start */
            for (size_t i = 0; i < d; i++) { mpfr_set_d(tmp, hs * fcur_d[i], MPFR_RNDN);
                mpfr_add(predm[i], yh[i], tmp, MPFR_RNDN); }
            mpfr_set_d(coef, hs, MPFR_RNDN);                 /* c0=1/hs -> coef=hs */
            for (size_t i = 0; i < d; i++) mpfr_set(base[i], yh[i], MPFR_RNDN);
            ok = mp_newton_bdf(P, t1m, base, coef, predm, ynext, tol, bits);
        } else {
            q = qtarget; if (q > m - 1) q = m - 1; if (q > MP_QMAX) q = MP_QMAX; if (q < 1) q = 1;
            mp_bdf_coeffs(t1m, th, q, cc, bits);
            for (size_t i = 0; i < d; i++) {                 /* base = -(1/c0) sum c_j y_{n+1-j} */
                mpfr_set_zero(acc, 1);
                for (int j = 1; j <= q; j++) { mpfr_mul(tmp, cc[j], yh[(size_t)(j-1)*d + i], MPFR_RNDN); mpfr_add(acc, acc, tmp, MPFR_RNDN); }
                mpfr_div(acc, acc, cc[0], MPFR_RNDN); mpfr_neg(base[i], acc, MPFR_RNDN);
            }
            mpfr_ui_div(coef, 1, cc[0], MPFR_RNDN);          /* coef = 1/c0        */
            mp_predict(t1m, th, yh, q, d, predm, bits);
            ok = mp_newton_bdf(P, t1m, base, coef, predm, ynext, tol, bits);
        }

        if (!ok) { no_grow = true; if (qtarget > 1) qtarget--; h = 0.5 * hs;
                   if (fabs(h) < 1e-14 * (fabs(td) + 1.0)) { st = ND_ERR_NONCONV; break; } continue; }

        /* MPFR local-error estimate: corrector - predictor (no cancellation). */
        for (size_t i = 0; i < d; i++) mpfr_sub(diffm[i], ynext[i], predm[i], MPFR_RNDN);
        double err = mp_wrms(d, diffm, yh, ynext, tol);
        double qexp = (double)q + 1.0;

        if (err <= 1.0) {
            mpfr_set(t, t1m, MPFR_RNDN);
            if (!nd_rhs_mpfr(P, t, ynext, fnext, bits)) { st = ND_ERR_SAMPLE; break; }
            mpsol_push(sol, t, ynext, fnext);
            if (o->step_monitor) { Expr* mo = eval_and_free(expr_copy(o->step_monitor)); expr_free(mo); }
            int qnew = q; double fac = 1.0;
            if (no_grow) { fac = 1.0; hold = q + 1; }
            else if (hold > 0) { fac = 1.0; hold--; }
            else {
                double fac_same = pow(1.0 / (err > 1e-300 ? err : 1e-300), 1.0 / qexp);
                if (q > 1) {
                    mp_predict(t1m, th, yh, q - 1, d, predm2, bits);
                    for (size_t i = 0; i < d; i++) mpfr_sub(diffm[i], ynext[i], predm2[i], MPFR_RNDN);
                    double err_dn = mp_wrms(d, diffm, yh, ynext, tol);
                    double fac_dn = pow(1.0 / (err_dn > 1e-300 ? err_dn : 1e-300), 1.0 / (double)q);
                    if (fac_dn > 1.5 * fac_same) qnew = q - 1;
                }
                if (qnew == q && q < MP_QMAX && q <= m - 1) qnew = q + 1;
                fac = 0.9 * fac_same;
                double cap = mp_grow_cap(qnew);
                if (fac < 0.2) fac = 0.2;
                if (fac > cap) fac = cap;
                hold = qnew + 1;
            }
            /* shift history ring */
            int keep = (m < CAP) ? m : CAP - 1;
            for (int j = keep; j >= 1; j--) {
                mpfr_set(th[j], th[j-1], MPFR_RNDN);
                for (size_t i = 0; i < d; i++) mpfr_set(yh[(size_t)j*d+i], yh[(size_t)(j-1)*d+i], MPFR_RNDN);
            }
            mpfr_set(th[0], t, MPFR_RNDN);
            for (size_t i = 0; i < d; i++) mpfr_set(yh[i], ynext[i], MPFR_RNDN);
            if (m < CAP) m++;
            qtarget = qnew; no_grow = false;
            h = dir * fabs(hs) * fac;
        } else {
            no_grow = true; hold = 0;
            if (qtarget > 1) qtarget--;
            double fac = 0.9 * pow(1.0 / err, 1.0 / qexp);
            if (fac < 0.2) fac = 0.2;
            if (fac > 1.0) fac = 1.0;
            h = dir * fabs(hs) * fac;
        }
    }
done:
    mp_vec_free(yh, (size_t)CAP * d); mp_vec_free(th, (size_t)CAP); mp_vec_free(cc, MP_QMAX + 1);
    mp_vec_free(base, d); mp_vec_free(ynext, d); mp_vec_free(fnext, d);
    mp_vec_free(predm, d); mp_vec_free(predm2, d); mp_vec_free(diffm, d);
    mpfr_clear(t); mpfr_clear(t1m); mpfr_clear(coef); mpfr_clear(acc); mpfr_clear(tmp);
    free(y0d); free(fcur_d);
    return st;
}

/* Integrate one direction toward target with the chosen explicit mpfr method. */
static NdStatus mpfr_dir(NdProblem* P, const NdOpts* o, MpSol* sol, NdTol tol,
                         bool adaptive, double target, long bits, int64_t* budget) {
    size_t d = P->d;
    double dir = (target > P->t0) ? 1.0 : -1.0;
    double span = fabs(P->tmax - P->tmin);
    if (fabs(target - P->t0) <= 1e-15 * (span + 1.0)) return ND_OK;

    mpfr_t t, h, tgt, hmag, tmpr;
    mpfr_init2(t, bits); mpfr_init2(h, bits); mpfr_init2(tgt, bits);
    mpfr_init2(hmag, bits); mpfr_init2(tmpr, bits);
    mpfr_set_d(t, P->t0, MPFR_RNDN); mpfr_set_d(tgt, target, MPFR_RNDN);
    double h0 = (o->starting_step > 0.0) ? o->starting_step : nd_fixed_step(P, o, tol, 1.0);
    mpfr_set_d(h, dir * fabs(h0), MPFR_RNDN);

    mpfr_t* Y  = mp_vec(d, bits);
    mpfr_t* Yn = mp_vec(d, bits);
    mpfr_t* Ye = mp_vec(d, bits);
    mpfr_t* fn = mp_vec(d, bits);
    mpfr_t* k  = mp_vec(7 * d, bits);
    for (size_t i = 0; i < d; i++) mpfr_set_d(Y[i], P->Y0[i], MPFR_RNDN);

    double h_cap = (o->max_step_size > 0.0) ? o->max_step_size : HUGE_VAL;
    double frac_cap = (o->max_step_fraction > 0.0) ? o->max_step_fraction * span : HUGE_VAL;
    NdStatus st = ND_OK;

    while (dir * (target - mpfr_get_d(t, MPFR_RNDN)) > 1e-14 * (fabs(mpfr_get_d(t, MPFR_RNDN)) + 1.0)) {
        double hm = fabs(mpfr_get_d(h, MPFR_RNDN));
        if (hm > h_cap) hm = h_cap;
        if (hm > frac_cap) hm = frac_cap;
        mpfr_set_d(h, dir * hm, MPFR_RNDN);
        mpfr_sub(tmpr, tgt, t, MPFR_RNDN);           /* remaining */
        if (dir * (mpfr_get_d(h, MPFR_RNDN) - mpfr_get_d(tmpr, MPFR_RNDN)) > 0.0) mpfr_set(h, tmpr, MPFR_RNDN);
        if (fabs(mpfr_get_d(h, MPFR_RNDN)) < 1e-300) { st = ND_ERR_STEPSIZE; break; }

        bool ok; double err = 0.0;
        if (adaptive) ok = dopri5_mpfr(P, t, Y, h, Yn, Ye, k, bits);
        else          ok = rk4_mpfr(P, t, Y, h, Yn, k, bits);
        if (--(*budget) < 0) { st = ND_ERR_MAXSTEPS; break; }
        if (!ok) { mpfr_mul_d(h, h, 0.5, MPFR_RNDN);
                   if (fabs(mpfr_get_d(h, MPFR_RNDN)) < 1e-300) { st = ND_ERR_STEPSIZE; break; } continue; }
        if (adaptive) err = mp_wrms(d, Ye, Y, Yn, tol);
        int q = adaptive ? 5 : 5;
        if (!adaptive || err <= 1.0) {
            mpfr_add(t, t, h, MPFR_RNDN);
            for (size_t i = 0; i < d; i++) mpfr_set(Y[i], Yn[i], MPFR_RNDN);
            if (!nd_rhs_mpfr(P, t, Y, fn, bits)) { st = ND_ERR_SAMPLE; break; }
            mpsol_push(sol, t, Y, fn);
            if (o->step_monitor) { Expr* m = eval_and_free(expr_copy(o->step_monitor)); expr_free(m); }
            double fac = adaptive ? (err > 0.0 ? 0.9 * pow(1.0/err, 1.0/q) : 5.0) : 1.0;
            if (fac < 0.2) fac = 0.2;
            if (fac > 5.0) fac = 5.0;
            mpfr_mul_d(h, h, fac, MPFR_RNDN);
        } else {
            double fac = 0.9 * pow(1.0/err, 1.0/q);
            if (fac < 0.2) fac = 0.2;
            if (fac > 1.0) fac = 1.0;
            mpfr_mul_d(h, h, fac, MPFR_RNDN);
        }
    }
    mp_vec_free(Y, d); mp_vec_free(Yn, d); mp_vec_free(Ye, d); mp_vec_free(fn, d); mp_vec_free(k, 7*d);
    mpfr_clear(t); mpfr_clear(h); mpfr_clear(tgt); mpfr_clear(hmag); mpfr_clear(tmpr);
    return st;
}

/* Build the mpfr InterpolatingFunction rule list from the node store. */
static Expr* mpfr_build_result(NdProblem* P, const MpSol* sol, long out_bits) {
    if (sol->n < 2) return NULL;
    size_t n = sol->n, d = sol->d;
    Expr** rules = malloc(sizeof(Expr*) * P->nfun);
    size_t nr = 0;
    for (size_t kf = 0; kf < P->nfun; kf++) {
        size_t comp = P->fun_state0[kf];
        Expr** entries = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            mpfr_t rt, rv, rd; mpfr_init2(rt, out_bits); mpfr_init2(rv, out_bits); mpfr_init2(rd, out_bits);
            mpfr_set(rt, sol->t[i], MPFR_RNDN);
            mpfr_set(rv, sol->Y[i*d + comp], MPFR_RNDN);
            mpfr_set(rd, sol->dY[i*d + comp], MPFR_RNDN);
            Expr* coord_el = expr_new_mpfr_copy(rt);
            Expr* coord = expr_new_function(expr_new_symbol(SYM_List), &coord_el, 1);
            Expr* trip[3] = { coord, expr_new_mpfr_copy(rv), expr_new_mpfr_copy(rd) };
            entries[i] = expr_new_function(expr_new_symbol(SYM_List), trip, 3);
            mpfr_clear(rt); mpfr_clear(rv); mpfr_clear(rd);
        }
        Expr* data = expr_new_function(expr_new_symbol(SYM_List), entries, n);
        free(entries);
        Expr* ifun = eval_and_free(expr_new_function(expr_new_symbol(SYM_Interpolation), &data, 1));
        if (!head_is(ifun, SYM_InterpolatingFunction)) {
            expr_free(ifun);
            for (size_t q = 0; q < nr; q++) expr_free(rules[q]);
            free(rules); return NULL;
        }
        Expr* lhs;
        if (P->fun_applied) {
            Expr* xa[1] = { expr_new_symbol(P->tvar) };
            lhs = expr_new_function(expr_new_symbol(P->fun_names[kf]), xa, 1);
            Expr* xa2[1] = { expr_new_symbol(P->tvar) };
            ifun = expr_new_function(ifun, xa2, 1);
        } else lhs = expr_new_symbol(P->fun_names[kf]);
        Expr* rargs[2] = { lhs, ifun };
        rules[nr++] = expr_new_function(expr_new_symbol(SYM_Rule), rargs, 2);
    }
    Expr* inner = expr_new_function(expr_new_symbol(SYM_List), rules, nr);
    free(rules);
    Expr* outer = expr_new_function(expr_new_symbol(SYM_List), &inner, 1);
    return outer;
}

/* Evaluate a (bound) boundary expression to mpfr at the current binding. */
static bool mp_eval_expr(NdProblem* P, Expr* e, mpfr_t out, long bits) {
    arith_warnings_mute_push();
    Expr* raw = eval_and_free(expr_copy(e));
    arith_warnings_mute_pop();
    if (!raw) return false;
    Expr* num = numericalize(raw, P->spec);
    expr_free(raw);
    mpfr_t im; mpfr_init2(im, bits);
    bool inexact;
    bool ok = num && get_approx_mpfr(num, out, im, &inexact) && mpfr_number_p(out);
    mpfr_clear(im);
    expr_free(num);
    return ok;
}

/* Run the shared mpfr integration of P into `sol`; returns the status. */
static NdStatus mpfr_integrate(NdProblem* P, const NdOpts* o, const NdStepper* S,
                               MpSol* sol, long bits) {
    NdTol tol = nd_resolve_tol(o);
    /* Stiff (implicit / multistep) methods use the MPFR BDF; explicit methods use
     * the DOPRI5/RK4 path.  RK4 alone is the fixed-step explicit request. */
    bool implicit = S && (S->flags & (ND_IMPLICIT | ND_MULTISTEP));
    bool adaptive = !(S && S->name && strcmp(S->name, "RK4") == 0);
    mpfr_t t0m; mpfr_init2(t0m, bits); mpfr_set_d(t0m, P->t0, MPFR_RNDN);
    mpfr_t* Y0 = mp_vec(P->d, bits);
    mpfr_t* f0 = mp_vec(P->d, bits);
    for (size_t i = 0; i < P->d; i++) mpfr_set_d(Y0[i], P->Y0[i], MPFR_RNDN);
    NdStatus st = ND_OK;
    if (!nd_rhs_mpfr(P, t0m, Y0, f0, bits)) st = ND_ERR_SAMPLE;
    else {
        P->tol.rtol = tol.rtol; P->tol.atol = tol.atol;
        mpsol_push(sol, t0m, Y0, f0);
        int64_t budget = (o->max_steps > 0) ? o->max_steps : 10000;
        if (P->tmax > P->t0)
            st = implicit ? mpfr_bdf_dir(P, o, sol, tol, P->tmax, bits, &budget)
                          : mpfr_dir(P, o, sol, tol, adaptive, P->tmax, bits, &budget);
        if (P->tmin < P->t0) {
            NdStatus s2 = implicit ? mpfr_bdf_dir(P, o, sol, tol, P->tmin, bits, &budget)
                                   : mpfr_dir(P, o, sol, tol, adaptive, P->tmin, bits, &budget);
            if (st == ND_OK) st = s2;
        }
    }
    mp_vec_free(Y0, P->d); mp_vec_free(f0, P->d); mpfr_clear(t0m);
    mpsol_sort(sol);
    return st;
}

Expr* nd_solve_mpfr_mol(NdProblem* P, const NdOpts* o, const NdStepper* S,
                        const NdMolGrid* g) {
    long out_bits = o->wp_bits > 53 ? o->wp_bits : 53;
    long bits = out_bits + 64;
    MpSol sol; mpsol_init(&sol, P->d, bits);
    NdStatus st = mpfr_integrate(P, o, S, &sol, bits);
    (void)st;
    size_t n = sol.n, d = sol.d, nx = g->nx;
    Expr* result = NULL;
    if (n >= 2) {
        size_t npts = n * nx;
        Expr** entries = malloc(sizeof(Expr*) * npts);
        size_t p = 0;
        mpfr_t xg, val, hb; mpfr_init2(xg, out_bits); mpfr_init2(val, out_bits); mpfr_init2(hb, out_bits);
        mpfr_set_d(hb, g->h, MPFR_RNDN);
        for (size_t i = 0; i < n; i++) {
            /* bind t and all states for boundary-value evaluation */
            nd_bind_set(&P->bind_t, expr_new_mpfr_copy(sol.t[i]));
            for (size_t k = 0; k < d; k++) nd_bind_set(&P->bind_y[k], expr_new_mpfr_copy(sol.Y[i*d + k]));
            for (size_t gi = 0; gi < nx; gi++) {
                mpfr_mul_ui(xg, hb, (unsigned long)gi, MPFR_RNDN);
                mpfr_add_d(xg, xg, g->xmin, MPFR_RNDN);
                if (g->periodic) {
                    size_t gg = (gi == nx - 1) ? 0 : gi;
                    mpfr_set(val, sol.Y[i*d + gg * (size_t)g->torder + 0], MPFR_RNDN);
                } else if (gi == 0) {
                    if (!mp_eval_expr(P, g->bc_left, val, bits)) mpfr_set_zero(val, 1);
                } else if (gi == nx - 1) {
                    if (!mp_eval_expr(P, g->bc_right, val, bits)) mpfr_set_zero(val, 1);
                } else {
                    mpfr_set(val, sol.Y[i*d + (gi - 1) * (size_t)g->torder + 0], MPFR_RNDN);
                }
                Expr* coord[2] = { expr_new_mpfr_copy(sol.t[i]), expr_new_mpfr_copy(xg) };
                Expr* coordL = expr_new_function(expr_new_symbol(SYM_List), coord, 2);
                Expr* pair[2] = { coordL, expr_new_mpfr_copy(val) };
                entries[p++] = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
            }
        }
        mpfr_clear(xg); mpfr_clear(val); mpfr_clear(hb);
        Expr* data = expr_new_function(expr_new_symbol(SYM_List), entries, npts);
        free(entries);
        Expr* ifun = eval_and_free(expr_new_function(expr_new_symbol(SYM_Interpolation), &data, 1));
        if (head_is(ifun, SYM_InterpolatingFunction)) {
            Expr* lhs;
            if (g->applied) {
                Expr* ar[2] = { expr_new_symbol(g->tvar), expr_new_symbol(g->xvar) };
                lhs = expr_new_function(expr_new_symbol(g->fname), ar, 2);
                Expr* ar2[2] = { expr_new_symbol(g->tvar), expr_new_symbol(g->xvar) };
                ifun = expr_new_function(ifun, ar2, 2);
            } else lhs = expr_new_symbol(g->fname);
            Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                                           (Expr*[]){ lhs, ifun }, 2);
            Expr* inner = expr_new_function(expr_new_symbol(SYM_List), &rule, 1);
            result = expr_new_function(expr_new_symbol(SYM_List), &inner, 1);
        } else expr_free(ifun);
    }
    mpsol_free(&sol);
    return result;
}

Expr* nd_solve_mpfr(NdProblem* P, const NdOpts* o, const NdStepper* S) {
    long out_bits = o->wp_bits > 53 ? o->wp_bits : 53;
    long bits = out_bits + 64;               /* guard digits */
    MpSol sol; mpsol_init(&sol, P->d, bits);
    NdStatus st = mpfr_integrate(P, o, S, &sol, bits);   /* routes explicit/implicit */
    (void)st;
    Expr* result = mpfr_build_result(P, &sol, out_bits);
    mpsol_free(&sol);
    return result;
}

#else
Expr* nd_solve_mpfr(NdProblem* P, const NdOpts* o, const NdStepper* S) {
    (void)P; (void)o; (void)S; return NULL;
}
#endif
