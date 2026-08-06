/* Mathilda — NDSolve Adams predictor–corrector (order 2 PECE: AB2 predictor +
 * trapezoidal corrector), self-started with one RK4 step.  Adaptive variable
 * step size: the (corrector - predictor) difference is a free Milne local-error
 * estimate driving accept/reject.  Adams is an explicit multistep family;
 * higher orders are the documented extension. */
#include "ndsolve_common.h"
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>

/* One RK4 step (self-start), local to this module. */
static bool adams_rk4(NdProblem* P, double t, const double* Y, double h, double* Ynew) {
    size_t d = P->d;
    double* k1 = malloc(sizeof(double) * d);
    double* k2 = malloc(sizeof(double) * d);
    double* k3 = malloc(sizeof(double) * d);
    double* k4 = malloc(sizeof(double) * d);
    double* tmp = malloc(sizeof(double) * d);
    bool ok = nd_rhs_real(P, t, Y, k1);
    if (ok) { for (size_t i = 0; i < d; i++) tmp[i] = Y[i] + 0.5*h*k1[i]; ok = nd_rhs_real(P, t+0.5*h, tmp, k2); }
    if (ok) { for (size_t i = 0; i < d; i++) tmp[i] = Y[i] + 0.5*h*k2[i]; ok = nd_rhs_real(P, t+0.5*h, tmp, k3); }
    if (ok) { for (size_t i = 0; i < d; i++) tmp[i] = Y[i] + h*k3[i];     ok = nd_rhs_real(P, t+h, tmp, k4); }
    if (ok) for (size_t i = 0; i < d; i++) Ynew[i] = Y[i] + (h/6.0)*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
    free(k1); free(k2); free(k3); free(k4); free(tmp);
    return ok;
}

/* Adaptive PECE loop.  Once history exists, the AB2 predictor and trapezoidal
 * corrector are both order-2, so their difference (corrector - predictor) is a
 * ready-made Milne local-error estimate — free, since both are already formed.
 * A WRMS norm <= 1 accepts the step; otherwise the step is shrunk and retried.
 * The RK4 self-start step is accepted as-is (no cheap same-order estimate). */
static NdStatus adams_dir(NdProblem* P, const NdOpts* o, NdSolution* sol, NdTol tol,
                          double target, int64_t* budget) {
    size_t d = P->d;
    double dir = (target > P->t0) ? 1.0 : -1.0;
    double span = fabs(P->tmax - P->tmin);
    if (fabs(target - P->t0) <= 16.0 * DBL_EPSILON * (span + 1.0)) return ND_OK;

    double* ycur  = malloc(sizeof(double) * d);
    double* ynext = malloc(sizeof(double) * d);
    double* fprev = malloc(sizeof(double) * d);   /* f_{n-1} */
    double* fcur  = malloc(sizeof(double) * d);   /* f_n     */
    double* pred  = malloc(sizeof(double) * d);
    double* fp    = malloc(sizeof(double) * d);
    memcpy(ycur, P->Y0, sizeof(double) * d);
    double t = P->t0;
    NdStatus st = ND_OK;
    bool have_hist = false;

    if (!nd_rhs_real(P, t, ycur, fcur)) { st = ND_ERR_SAMPLE; goto done; }

    double h = (o->starting_step > 0.0) ? dir * o->starting_step
                                        : nd_initial_step(P, o, tol, t, ycur, fcur, 4, dir);
    if (h == 0.0) h = dir * span / 100.0;
    double h_prev = 0.0;        /* last accepted step, for the variable-step AB2 */
    bool   no_grow = false;
    double h_cap = (o->max_step_size > 0.0) ? o->max_step_size : HUGE_VAL;
    double frac_cap = (o->max_step_fraction > 0.0) ? o->max_step_fraction * span : HUGE_VAL;

    while ((target - t) * dir > 16.0 * DBL_EPSILON * (fabs(t) + 1.0)) {
        double hmag = fabs(h);
        if (hmag > h_cap) hmag = h_cap;
        if (hmag > frac_cap) hmag = frac_cap;
        double hs = dir * hmag;
        if ((t + hs - target) * dir > 0.0) hs = target - t;
        if (fabs(hs) < 16.0 * DBL_EPSILON * (fabs(t) + 1.0)) { st = ND_ERR_STEPSIZE; break; }
        if (--(*budget) < 0) { st = ND_ERR_MAXSTEPS; break; }

        if (!have_hist) {
            /* RK4 self-start for the very first step (accepted unconditionally). */
            if (!adams_rk4(P, t, ycur, hs, ynext)) { st = ND_ERR_SAMPLE; break; }
            t += hs;
            double* fnext = fp;
            if (!nd_rhs_real(P, t, ynext, fnext)) { st = ND_ERR_SAMPLE; break; }
            nd_solution_push(sol, t, ynext, fnext);
            if (o->step_monitor) { Expr* m = eval_and_free(expr_copy(o->step_monitor)); expr_free(m); }
            memcpy(fprev, fcur, sizeof(double) * d);
            memcpy(fcur, fnext, sizeof(double) * d);
            memcpy(ycur, ynext, sizeof(double) * d);
            h_prev = hs; have_hist = true; no_grow = false;
            continue;
        }

        /* Variable-step AB2 predictor: y_n + h[(1 + w/2) f_n - (w/2) f_{n-1}],
         * w = h/h_prev (w=1 recovers 3/2, -1/2).  Trapezoidal corrector.
         * fcur is maintained as f(t, ycur) across accept/reject, so no re-eval. */
        double w = hs / h_prev;
        for (size_t i = 0; i < d; i++)
            pred[i] = ycur[i] + hs*((1.0 + 0.5*w)*fcur[i] - (0.5*w)*fprev[i]);
        if (!nd_rhs_real(P, t + hs, pred, fp)) { st = ND_ERR_SAMPLE; break; }
        for (size_t i = 0; i < d; i++) ynext[i] = ycur[i] + 0.5*hs*(fcur[i] + fp[i]);

        /* Milne error estimate = corrector - predictor (both order 2). */
        double err = 0.0;
        { double* e = pred;   /* reuse pred as the difference buffer */
          for (size_t i = 0; i < d; i++) e[i] = ynext[i] - pred[i];
          err = nd_wrms_norm(d, e, ycur, ynext, tol); }

        if (err > 1.0) {
            /* reject: shrink and retry (do not advance). */
            double fac = 0.9 * pow(1.0/err, 1.0/3.0);
            if (fac < 0.2) fac = 0.2;
            if (fac > 1.0) fac = 1.0;
            h = dir * fabs(hs) * fac; no_grow = true;
            continue;
        }

        /* accept */
        t += hs;
        double* fnext = fp;
        if (!nd_rhs_real(P, t, ynext, fnext)) { st = ND_ERR_SAMPLE; break; }
        nd_solution_push(sol, t, ynext, fnext);
        if (o->step_monitor) { Expr* m = eval_and_free(expr_copy(o->step_monitor)); expr_free(m); }
        memcpy(fprev, fcur, sizeof(double) * d);
        memcpy(fcur, fnext, sizeof(double) * d);
        memcpy(ycur, ynext, sizeof(double) * d);
        h_prev = hs;
        double fac = (err > 0.0) ? 0.9 * pow(1.0/err, 1.0/3.0) : 5.0;
        if (fac < 0.2) fac = 0.2;
        if (fac > (no_grow ? 1.0 : 5.0)) fac = no_grow ? 1.0 : 5.0;
        no_grow = false;
        h = dir * fabs(hs) * fac;
    }
done:
    free(ycur); free(ynext); free(fprev); free(fcur); free(pred); free(fp);
    return st;
}

NdStatus nd_multistep_adams(NdProblem* P, const NdOpts* o, NdSolution* sol) {
    size_t d = P->d;
    NdTol tol = nd_resolve_tol(o);
    P->tol.rtol = tol.rtol; P->tol.atol = tol.atol;
    double* f0 = malloc(sizeof(double) * d);
    if (!nd_rhs_real(P, P->t0, P->Y0, f0)) { free(f0); return ND_ERR_SAMPLE; }
    nd_solution_push(sol, P->t0, P->Y0, f0);
    free(f0);
    int64_t budget = (o->max_steps > 0) ? o->max_steps : ND_AUTO_MAX_STEPS;
    NdStatus st = ND_OK;
    if (P->tmax > P->t0) st = adams_dir(P, o, sol, tol, P->tmax, &budget);
    if (P->tmin < P->t0) { NdStatus s2 = adams_dir(P, o, sol, tol, P->tmin, &budget);
                           if (st == ND_OK) st = s2; }
    return st;
}

const NdStepper nd_stepper_adams = {
    "Adams", ND_MULTISTEP, 2, 1, 0, NULL,
    "NDSolve`Adams[eqns, u, {x, xmin, xmax}]\n"
    "\tAdams predictor–corrector multistep method (order 2 PECE: AB2 predictor,\n"
    "\ttrapezoidal corrector), self-started with RK4.  Adaptive variable step\n"
    "\tsize via the predictor-corrector local error estimate."
};
