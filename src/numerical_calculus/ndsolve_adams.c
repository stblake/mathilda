/* Mathilda — NDSolve Adams predictor–corrector (order 2 PECE: AB2 predictor +
 * trapezoidal corrector), self-started with one RK4 step.  Adams is an explicit
 * multistep family; higher orders are the documented extension. */
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

static NdStatus adams_dir(NdProblem* P, const NdOpts* o, NdSolution* sol, NdTol tol,
                          double target, int64_t* budget) {
    size_t d = P->d;
    double dir = (target > P->t0) ? 1.0 : -1.0;
    double span = fabs(P->tmax - P->tmin);
    if (fabs(target - P->t0) <= 16.0 * DBL_EPSILON * (span + 1.0)) return ND_OK;
    double h = nd_fixed_step(P, o, tol, dir);
    if (h == 0.0) h = dir * span / 100.0;

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

    while ((target - t) * dir > 16.0 * DBL_EPSILON * (fabs(t) + 1.0)) {
        double hs = h;
        if ((t + hs - target) * dir > 0.0) hs = target - t;
        if (--(*budget) < 0) { st = ND_ERR_MAXSTEPS; break; }
        if (!have_hist) {
            /* RK4 self-start for the very first step. */
            if (!adams_rk4(P, t, ycur, hs, ynext)) { st = ND_ERR_SAMPLE; break; }
        } else {
            /* AB2 predictor then trapezoidal corrector (PECE). */
            for (size_t i = 0; i < d; i++) pred[i] = ycur[i] + hs*(1.5*fcur[i] - 0.5*fprev[i]);
            if (!nd_rhs_real(P, t + hs, pred, fp)) { st = ND_ERR_SAMPLE; break; }
            for (size_t i = 0; i < d; i++) ynext[i] = ycur[i] + 0.5*hs*(fcur[i] + fp[i]);
        }
        t += hs;
        double* fnext = fp;   /* reuse buffer */
        if (!nd_rhs_real(P, t, ynext, fnext)) { st = ND_ERR_SAMPLE; break; }
        nd_solution_push(sol, t, ynext, fnext);
        if (o->step_monitor) { Expr* m = eval_and_free(expr_copy(o->step_monitor)); expr_free(m); }
        memcpy(fprev, fcur, sizeof(double) * d);
        memcpy(fcur, fnext, sizeof(double) * d);
        memcpy(ycur, ynext, sizeof(double) * d);
        have_hist = true;
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
    int64_t budget = (o->max_steps > 0) ? o->max_steps : 10000;
    NdStatus st = ND_OK;
    if (P->tmax > P->t0) st = adams_dir(P, o, sol, tol, P->tmax, &budget);
    if (P->tmin < P->t0) { NdStatus s2 = adams_dir(P, o, sol, tol, P->tmin, &budget);
                           if (st == ND_OK) st = s2; }
    return st;
}

const NdStepper nd_stepper_adams = {
    "Adams", ND_MULTISTEP, 2, 0, 0, NULL,
    "NDSolve`Adams[eqns, u, {x, xmin, xmax}]\n"
    "\tAdams predictor–corrector multistep method (order 2 PECE: AB2 predictor,\n"
    "\ttrapezoidal corrector), self-started with RK4."
};
