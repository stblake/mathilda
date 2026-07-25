/* Mathilda — NDSolve implicit / stiff steppers: BackwardEuler, ImplicitTrapezoid
 * (one-step, via the shared theta-method Newton solve) and BDF (adaptive
 * variable-step backward differentiation, orders 1-2, with local error control
 * and Newton-failure step recovery; self-started at order 1). */
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

/* ------------------ BDF multistep (variable step, orders 1-2) ---------------- *
 *
 * Adaptive variable-step-size BDF, self-started at order 1 (backward Euler) and
 * rising to order 2 once a second node exists.  Three ingredients make it robust
 * on the stiff / awkward problems fixed-step BDF cannot handle:
 *
 *   1. Local error control (Milne device).  Each step forms an explicit
 *      predictor of the SAME order as the BDF corrector; the (corrector -
 *      predictor) difference is a cheap local truncation-error estimate (no
 *      extra RHS/Newton work).  A WRMS norm <= 1 accepts the step; otherwise it
 *      is rejected and the step shrunk.
 *        order 1:  predictor = explicit Euler  y_n + h f_n
 *        order 2:  predictor = Hermite quadratic through (t_{n-1}, y_{n-1}) and
 *                  (t_n, y_n, f_n), extrapolated to t_{n+1}.
 *   2. Newton-failure recovery.  When nd_newton_theta diverges (the classic
 *      failure at incompatible IC/BC corners, `ndcf`), the step is halved and
 *      retried instead of aborting; after repeated trouble the order is dropped
 *      to the L-stable backward Euler.  Only a genuine step-size collapse is
 *      terminal.
 *   3. Variable-step BDF2 coefficients.  With ratio w = h_n / h_{n-1},
 *        a*y_{n+1} - b*y_n + c*y_{n-1} = h_n f(t_{n+1}, y_{n+1}),
 *        a=(1+2w)/(1+w),  b=(1+w),  c=w^2/(1+w),
 *      solved as the theta-method  y = Ybase + (h_n/a) f  with
 *        Ybase=(b/a)y_n-(c/a)y_{n-1},  theta=1/a  (w=1 recovers 4/3,1/3,2/3).
 */
static NdStatus bdf_dir(NdProblem* P, const NdOpts* o, NdSolution* sol, NdTol tol,
                        double target, int64_t* budget) {
    size_t d = P->d;
    double dir = (target > P->t0) ? 1.0 : -1.0;
    double span = fabs(P->tmax - P->tmin);
    if (fabs(target - P->t0) <= 16.0 * DBL_EPSILON * (span + 1.0)) return ND_OK;

    double* yprev = malloc(sizeof(double) * d);   /* y_{n-1}                    */
    double* ycur  = malloc(sizeof(double) * d);   /* y_n                        */
    double* ynext = malloc(sizeof(double) * d);   /* trial y_{n+1}             */
    double* base  = malloc(sizeof(double) * d);   /* theta-method Ybase        */
    double* pred  = malloc(sizeof(double) * d);   /* explicit predictor        */
    double* fcur  = malloc(sizeof(double) * d);   /* f(t_n, y_n)               */
    double* fnext = malloc(sizeof(double) * d);   /* f(t_{n+1}, y_{n+1})       */

    memcpy(ycur, P->Y0, sizeof(double) * d);
    memcpy(yprev, P->Y0, sizeof(double) * d);
    double t = P->t0;
    NdStatus st = ND_OK;
    if (!nd_rhs_real(P, t, ycur, fcur)) { st = ND_ERR_SAMPLE; goto done; }

    /* starting step: Hairer heuristic (order 1), then adapted */
    double h = (o->starting_step > 0.0) ? dir * o->starting_step
                                        : nd_initial_step(P, o, tol, t, ycur, fcur, 1, dir);
    if (h == 0.0) h = dir * span / 100.0;

    double h_prev = 0.0;        /* size of the last accepted step (for w)       */
    bool   have_prev = false;   /* is y_{n-1} a genuine second history point?    */
    bool   no_grow = false;     /* just recovered from a reject/failure          */
    int    trouble = 0;         /* consecutive rejects/failures at this t        */
    double h_cap = (o->max_step_size > 0.0) ? o->max_step_size : HUGE_VAL;
    double frac_cap = (o->max_step_fraction > 0.0) ? o->max_step_fraction * span : HUGE_VAL;

    while ((target - t) * dir > 16.0 * DBL_EPSILON * (fabs(t) + 1.0)) {
        /* clamp magnitude, then land exactly on the target */
        double hmag = fabs(h);
        if (hmag > h_cap) hmag = h_cap;
        if (hmag > frac_cap) hmag = frac_cap;
        double hs = dir * hmag;
        if ((t + hs - target) * dir > 0.0) hs = target - t;
        if (fabs(hs) < 16.0 * DBL_EPSILON * (fabs(t) + 1.0)) { st = ND_ERR_STEPSIZE; break; }
        if (--(*budget) < 0) { st = ND_ERR_MAXSTEPS; break; }

        /* Order 2 only with a valid second point; after repeated trouble drop to
         * the very robust backward Euler for one attempt. */
        bool order2 = have_prev && trouble < 2;
        double q;               /* step-control exponent = order + 1            */
        bool ok;
        if (order2) {
            double w = hs / h_prev;
            double a = (1.0 + 2.0*w) / (1.0 + w);
            double b = 1.0 + w;
            double c = (w*w) / (1.0 + w);
            for (size_t i = 0; i < d; i++) {
                base[i] = (b/a) * ycur[i] - (c/a) * yprev[i];
                /* Hermite quadratic predictor through y_{n-1}, (y_n, f_n) */
                double A = (yprev[i] - ycur[i] + fcur[i]*h_prev) / (h_prev*h_prev);
                pred[i] = ycur[i] + fcur[i]*hs + A*hs*hs;
            }
            ok = nd_newton_theta(P, t + hs, base, hs, 1.0/a, NULL, pred, ynext, tol);
            q = 3.0;
        } else {
            for (size_t i = 0; i < d; i++) pred[i] = ycur[i] + hs*fcur[i]; /* Euler */
            ok = nd_newton_theta(P, t + hs, ycur, hs, 1.0, NULL, pred, ynext, tol);
            q = 2.0;
        }

        if (!ok) {
            /* Newton diverged: shrink and retry; escalate the order-1 fallback. */
            trouble++; no_grow = true;
            h = 0.5 * hs;
            if (fabs(h) < 16.0 * DBL_EPSILON * (fabs(t) + 1.0)) { st = ND_ERR_NONCONV; break; }
            continue;
        }

        /* Milne local-error estimate: corrector - predictor. */
        for (size_t i = 0; i < d; i++) pred[i] = ynext[i] - pred[i];
        double err = nd_wrms_norm(d, pred, ycur, ynext, tol);

        if (err <= 1.0) {
            /* accept */
            t += hs;
            if (!nd_rhs_real(P, t, ynext, fnext)) { st = ND_ERR_SAMPLE; break; }
            nd_solution_push(sol, t, ynext, fnext);
            if (o->step_monitor) { Expr* m = eval_and_free(expr_copy(o->step_monitor)); expr_free(m); }
            memcpy(yprev, ycur, sizeof(double) * d);
            memcpy(ycur,  ynext, sizeof(double) * d);
            memcpy(fcur,  fnext, sizeof(double) * d);
            h_prev = hs; have_prev = true; trouble = 0;
            double fac = (err > 0.0) ? 0.9 * pow(1.0/err, 1.0/q) : 5.0;
            if (fac < 0.2) fac = 0.2;
            if (fac > (no_grow ? 1.0 : 5.0)) fac = no_grow ? 1.0 : 5.0;
            no_grow = false;
            h = dir * fabs(hs) * fac;
        } else {
            /* reject: shrink, never grow */
            trouble++; no_grow = true;
            double fac = 0.9 * pow(1.0/err, 1.0/q);
            if (fac < 0.2) fac = 0.2;
            if (fac > 1.0) fac = 1.0;
            h = dir * fabs(hs) * fac;
        }
    }
done:
    free(yprev); free(ycur); free(ynext); free(base);
    free(pred); free(fcur); free(fnext);
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
    "BDF", ND_IMPLICIT | ND_MULTISTEP, 2, 1, 0, NULL,
    "NDSolve`BDF[eqns, u, {x, xmin, xmax}]\n"
    "\tBackward differentiation formula, a stiff implicit multistep method.\n"
    "\tAdaptive variable step size (orders 1-2) with predictor-corrector local\n"
    "\terror control and Newton-failure step recovery; self-started at order 1."
};
