/* Mathilda — NDSolve implicit / stiff steppers: BackwardEuler, ImplicitTrapezoid
 * (one-step, via the shared theta-method Newton solve) and BDF (variable-step,
 * order-2 backward differentiation, self-started with backward Euler). */
#include "ndsolve_common.h"
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>

static NdTol prob_tol(NdProblem* P) { NdTol t; t.rtol = P->tol.rtol; t.atol = P->tol.atol; return t; }

/* Backward Euler: Ynew = Yn + h f(t+h, Ynew).  theta = 1, Ybase = Yn. */
static bool backward_euler_step(const NdStepper* S, NdProblem* P, double t, const double* Y,
                                double h, double* Ynew, double* Yerr, double* K) {
    (void)S; (void)Yerr; (void)K;
    size_t d = P->d;
    double* f0 = malloc(sizeof(double) * d);
    double* guess = malloc(sizeof(double) * d);
    bool ok = nd_rhs_real(P, t, Y, f0);
    if (ok) { for (size_t i = 0; i < d; i++) guess[i] = Y[i] + h * f0[i];
              ok = nd_newton_theta(P, t + h, Y, h, 1.0, NULL, guess, Ynew, prob_tol(P)); }
    free(f0); free(guess);
    return ok;
}

/* Implicit trapezoid: Ynew = Yn + h/2 (f(tn,Yn) + f(t+h,Ynew)).
 * theta = 1/2, Ybase = Yn, rhs_const = h/2 f(tn,Yn). */
static bool implicit_trapezoid_step(const NdStepper* S, NdProblem* P, double t, const double* Y,
                                    double h, double* Ynew, double* Yerr, double* K) {
    (void)S; (void)Yerr; (void)K;
    size_t d = P->d;
    double* f0 = malloc(sizeof(double) * d);
    double* rc = malloc(sizeof(double) * d);
    double* guess = malloc(sizeof(double) * d);
    bool ok = nd_rhs_real(P, t, Y, f0);
    if (ok) {
        for (size_t i = 0; i < d; i++) { rc[i] = 0.5 * h * f0[i]; guess[i] = Y[i] + h * f0[i]; }
        ok = nd_newton_theta(P, t + h, Y, h, 0.5, rc, guess, Ynew, prob_tol(P));
    }
    free(f0); free(rc); free(guess);
    return ok;
}

/* -------------------- BDF multistep (order 2, fixed step) -------------------- */

/* Integrate one direction with BDF2, self-started by one backward-Euler step. */
static NdStatus bdf_dir(NdProblem* P, const NdOpts* o, NdSolution* sol, NdTol tol,
                        double target, int64_t* budget) {
    size_t d = P->d;
    double dir = (target > P->t0) ? 1.0 : -1.0;
    double span = fabs(P->tmax - P->tmin);
    if (fabs(target - P->t0) <= 16.0 * DBL_EPSILON * (span + 1.0)) return ND_OK;
    double h = nd_fixed_step(P, o, tol, dir);
    if (h == 0.0) h = dir * span / 100.0;

    double* yprev = malloc(sizeof(double) * d);   /* y_{n-1} */
    double* ycur  = malloc(sizeof(double) * d);   /* y_n     */
    double* ynext = malloc(sizeof(double) * d);
    double* base  = malloc(sizeof(double) * d);
    double* guess = malloc(sizeof(double) * d);
    double* f     = malloc(sizeof(double) * d);
    double* f0    = malloc(sizeof(double) * d);
    memcpy(ycur, P->Y0, sizeof(double) * d);
    memcpy(yprev, P->Y0, sizeof(double) * d);
    double t = P->t0;
    NdStatus st = ND_OK;
    bool first = true;

    while ((target - t) * dir > 16.0 * DBL_EPSILON * (fabs(t) + 1.0)) {
        double hs = h;
        if ((t + hs - target) * dir > 0.0) hs = target - t;   /* land on target */
        if (--(*budget) < 0) { st = ND_ERR_MAXSTEPS; break; }
        bool ok;
        if (first) {
            /* backward Euler startup */
            if (!nd_rhs_real(P, t, ycur, f0)) { st = ND_ERR_SAMPLE; break; }
            for (size_t i = 0; i < d; i++) guess[i] = ycur[i] + hs * f0[i];
            ok = nd_newton_theta(P, t + hs, ycur, hs, 1.0, NULL, guess, ynext, tol);
        } else {
            /* BDF2: ynext = (4/3)y_n - (1/3)y_{n-1} + (2/3) h f(t1, ynext) */
            for (size_t i = 0; i < d; i++) {
                base[i] = (4.0/3.0) * ycur[i] - (1.0/3.0) * yprev[i];
                guess[i] = ycur[i];
            }
            ok = nd_newton_theta(P, t + hs, base, hs, 2.0/3.0, NULL, guess, ynext, tol);
        }
        if (!ok) { st = ND_ERR_NONCONV; break; }
        t += hs;
        if (!nd_rhs_real(P, t, ynext, f)) { st = ND_ERR_SAMPLE; break; }
        nd_solution_push(sol, t, ynext, f);
        if (o->step_monitor) { Expr* m = eval_and_free(expr_copy(o->step_monitor)); expr_free(m); }
        memcpy(yprev, ycur, sizeof(double) * d);
        memcpy(ycur, ynext, sizeof(double) * d);
        first = false;
    }
    free(yprev); free(ycur); free(ynext); free(base); free(guess); free(f); free(f0);
    return st;
}

NdStatus nd_multistep_bdf(NdProblem* P, const NdOpts* o, NdSolution* sol) {
    size_t d = P->d;
    NdTol tol = nd_resolve_tol(o);
    P->tol.rtol = tol.rtol; P->tol.atol = tol.atol;
    double* f0 = malloc(sizeof(double) * d);
    if (!nd_rhs_real(P, P->t0, P->Y0, f0)) { free(f0); return ND_ERR_SAMPLE; }
    nd_solution_push(sol, P->t0, P->Y0, f0);
    free(f0);
    int64_t budget = (o->max_steps > 0) ? o->max_steps : 10000;
    NdStatus st = ND_OK;
    if (P->tmax > P->t0) st = bdf_dir(P, o, sol, tol, P->tmax, &budget);
    if (P->tmin < P->t0) { NdStatus s2 = bdf_dir(P, o, sol, tol, P->tmin, &budget);
                           if (st == ND_OK) st = s2; }
    return st;
}

const NdStepper nd_stepper_backward_euler = {
    "BackwardEuler", ND_IMPLICIT, 1, 0, 0, backward_euler_step,
    "NDSolve`BackwardEuler[eqns, u, {x, xmin, xmax}]\n"
    "\tImplicit (backward) Euler method, order 1; A-stable, robust for stiff\n"
    "\tproblems.  Each step solves a nonlinear system by Newton iteration."
};

const NdStepper nd_stepper_implicit_trapezoid = {
    "ImplicitTrapezoid", ND_IMPLICIT, 2, 0, 0, implicit_trapezoid_step,
    "NDSolve`ImplicitTrapezoid[eqns, u, {x, xmin, xmax}]\n"
    "\tImplicit trapezoidal rule, order 2, A-stable."
};

const NdStepper nd_stepper_bdf = {
    "BDF", ND_IMPLICIT | ND_MULTISTEP, 2, 0, 0, NULL,
    "NDSolve`BDF[eqns, u, {x, xmin, xmax}]\n"
    "\tBackward differentiation formula (order 2), a stiff multistep method\n"
    "\tself-started with backward Euler."
};
