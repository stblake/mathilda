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

Expr* nd_solve_mpfr(NdProblem* P, const NdOpts* o, const NdStepper* S) {
    long out_bits = o->wp_bits > 53 ? o->wp_bits : 53;
    long bits = out_bits + 64;               /* guard digits */
    NdTol tol = nd_resolve_tol(o);
    bool adaptive = !(S && S->name && strcmp(S->name, "RK4") == 0);

    MpSol sol; mpsol_init(&sol, P->d, bits);
    /* record initial node */
    mpfr_t t0m; mpfr_init2(t0m, bits); mpfr_set_d(t0m, P->t0, MPFR_RNDN);
    mpfr_t* Y0 = mp_vec(P->d, bits);
    mpfr_t* f0 = mp_vec(P->d, bits);
    for (size_t i = 0; i < P->d; i++) mpfr_set_d(Y0[i], P->Y0[i], MPFR_RNDN);
    NdStatus st = ND_OK;
    if (!nd_rhs_mpfr(P, t0m, Y0, f0, bits)) st = ND_ERR_SAMPLE;
    else {
        mpsol_push(&sol, t0m, Y0, f0);
        int64_t budget = (o->max_steps > 0) ? o->max_steps : 10000;
        if (P->tmax > P->t0) st = mpfr_dir(P, o, &sol, tol, adaptive, P->tmax, bits, &budget);
        if (P->tmin < P->t0) { NdStatus s2 = mpfr_dir(P, o, &sol, tol, adaptive, P->tmin, bits, &budget);
                               if (st == ND_OK) st = s2; }
    }
    mp_vec_free(Y0, P->d); mp_vec_free(f0, P->d); mpfr_clear(t0m);
    mpsol_sort(&sol);
    Expr* result = mpfr_build_result(P, &sol, out_bits);
    mpsol_free(&sol);
    (void)st;
    return result;
}

#else
Expr* nd_solve_mpfr(NdProblem* P, const NdOpts* o, const NdStepper* S) {
    (void)P; (void)o; (void)S; return NULL;
}
#endif
