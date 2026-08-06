/* Mathilda — NDSolve implicit / stiff steppers: BackwardEuler, ImplicitTrapezoid
 * (one-step, via the shared theta-method Newton solve) and BDF (variable-step
 * variable-order backward differentiation, orders 1-5, with exact nonuniform-mesh
 * coefficients, local error control, order ramping and Newton-failure step/order
 * recovery; self-started at order 1). */
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

/* ------------- BDF multistep: variable step, variable order 1-5 -------------- *
 *
 * A variable-step variable-order (VSVO) BDF integrator — the standard stiff
 * workhorse (LSODE/CVODE/DASSL family), giving high accuracy at few steps.
 *
 *  * Exact variable-step coefficients.  BDF-q enforces  L'(t_{n+1}) = f, where L
 *    is the degree-q polynomial interpolating (t_{n+1},y_{n+1}) and the q most
 *    recent nodes.  Differentiating the Lagrange basis at t_{n+1} gives the
 *    coefficients  c_j = L_j'(t_{n+1})  directly from the actual (nonuniform)
 *    node times — exact for any step distribution and any order (bdf_coeffs).
 *    The implicit relation  c_0 y_{n+1} + sum_{j>=1} c_j y_{n+1-j} = f  is solved
 *    as the shared theta-method  y = Ybase + (1/c_0) f  (Ybase = -(1/c_0) sum
 *    c_j y_{n+1-j}, theta = 1/(c_0 h)).  (q=1 recovers backward Euler; q=2, equal
 *    steps, recovers 4/3, 1/3, 2/3.)
 *  * Milne local-error estimate.  The predictor is the degree-q polynomial
 *    through the q+1 most recent nodes, extrapolated to t_{n+1} (lagrange_predict)
 *    — same order as the corrector, so (corrector - predictor) estimates the
 *    local truncation error with no extra RHS/Newton work.  WRMS <= 1 accepts.
 *  * Order control.  Order ramps up by one per successful step (up to the number
 *    of available history points and BDF_QMAX) and drops on trouble.  Step growth
 *    is capped tighter at higher order to respect variable-step BDF's zero-
 *    stability step-ratio bounds, and small step changes are suppressed (held at
 *    the current h) to keep the mesh smooth — both stabilize high-order BDF.
 *  * Newton-failure recovery.  A diverging Newton iteration (the `ndcf` failure at
 *    incompatible IC/BC corners) halves the step, drops the order toward the
 *    L-stable backward Euler, and retries; only a step-size collapse is terminal.
 */
#define BDF_QMAX 5

/* BDF-q coefficients c[0..q]: c_j is the derivative at t1 of the Lagrange basis
 * for node x_j, where x_0 = t1 (the unknown y_{n+1}) and x_{1..q} = th[0..q-1]
 * (the q most recent nodes, th[0] most recent).  Then
 *   c_0 y_{n+1} + sum_{j=1}^{q} c_j y_{n+1-j} = f(t_{n+1}, y_{n+1}). */
static void bdf_coeffs(double t1, const double* th, int q, double* c) {
    double x[BDF_QMAX + 1];
    x[0] = t1;
    for (int j = 1; j <= q; j++) x[j] = th[j - 1];
    double s = 0.0;
    for (int m = 1; m <= q; m++) s += 1.0 / (x[0] - x[m]);
    c[0] = s;                                     /* L_0'(t1) = sum 1/(x0-xm) */
    for (int j = 1; j <= q; j++) {
        double num = 1.0, den = 1.0;
        for (int m = 0; m <= q; m++) {
            if (m == j) continue;
            den *= (x[j] - x[m]);                 /* all m != j                */
            if (m == 0) continue;
            num *= (x[0] - x[m]);                 /* m != j and m != 0         */
        }
        c[j] = num / den;
    }
}

/* Degree-q polynomial through the q+1 most recent nodes (th[0..q], yh rows),
 * extrapolated to t1 -> pred[0..d).  th[0] is the most recent node. */
static void lagrange_predict(double t1, const double* th, const double* yh,
                             int q, size_t d, double* pred) {
    double w[BDF_QMAX + 1];
    for (int j = 0; j <= q; j++) {
        double num = 1.0, den = 1.0;
        for (int m = 0; m <= q; m++) {
            if (m == j) continue;
            num *= (t1 - th[m]);
            den *= (th[j] - th[m]);
        }
        w[j] = num / den;
    }
    for (size_t i = 0; i < d; i++) {
        double acc = 0.0;
        for (int j = 0; j <= q; j++) acc += w[j] * yh[(size_t)j * d + i];
        pred[i] = acc;
    }
}

/* Per-reassessment cap on a step jump.  The step is held constant between
 * reassessments (uniform mesh -> BDF zero-stable at every order <= 6), so the cap
 * only bounds the size of an occasional jump; it tightens with order because a
 * large ratio is riskier at high order.  Index by order 1..5. */
static double bdf_grow_cap(int q) {
    switch (q) {
        case 1: case 2: return 5.0;
        case 3:         return 4.0;
        case 4:         return 3.0;
        default:        return 2.0;   /* q >= 5 */
    }
}

static NdStatus bdf_dir(NdProblem* P, const NdOpts* o, NdSolution* sol, NdTol tol,
                        double target, int64_t* budget) {
    size_t d = P->d;
    double dir = (target > P->t0) ? 1.0 : -1.0;
    double span = fabs(P->tmax - P->tmin);
    if (fabs(target - P->t0) <= 16.0 * DBL_EPSILON * (span + 1.0)) return ND_OK;

    /* history ring: th[0]/yh row 0 = most recent node, up to BDF_QMAX+1 nodes */
    const int CAP = BDF_QMAX + 1;
    double  th[BDF_QMAX + 1];
    double* yh    = malloc(sizeof(double) * (size_t)CAP * d);
    double  cc[BDF_QMAX + 1];
    double* ynext = malloc(sizeof(double) * d);
    double* base  = malloc(sizeof(double) * d);
    double* pred  = malloc(sizeof(double) * d);
    double* fcur  = malloc(sizeof(double) * d);   /* f(t_n, y_n) (order-1 start) */
    double* fnext = malloc(sizeof(double) * d);

    th[0] = P->t0;
    memcpy(yh, P->Y0, sizeof(double) * d);
    int    m = 1;               /* number of stored history nodes               */
    double t = P->t0;
    NdStatus st = ND_OK;
    if (!nd_rhs_real(P, t, yh, fcur)) { st = ND_ERR_SAMPLE; goto done; }

    double h = (o->starting_step > 0.0) ? dir * o->starting_step
                                        : nd_initial_step(P, o, tol, t, yh, fcur, 1, dir);
    if (h == 0.0) h = dir * span / 100.0;

    int    qtarget = 1;         /* desired order                                 */
    int    hold = 0;            /* steps to hold (order, step) before reassessing */
    bool   no_grow = false;     /* just recovered from a reject/failure          */
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

        double t1 = t + hs;
        int q;                  /* order actually used this step                 */
        bool ok;
        if (m == 1) {
            /* first step: order 1 with explicit-Euler predictor */
            q = 1;
            for (size_t i = 0; i < d; i++) pred[i] = yh[i] + hs * fcur[i];
            ok = nd_newton_theta(P, t1, yh, hs, 1.0, NULL, pred, ynext, tol);
        } else {
            /* order-q needs q history nodes for the corrector and q+1 for the
             * same-order predictor: q <= m-1 (qtarget is dropped on trouble). */
            q = qtarget;
            if (q > m - 1) q = m - 1;
            if (q > BDF_QMAX) q = BDF_QMAX;
            if (q < 1) q = 1;
            bdf_coeffs(t1, th, q, cc);
            for (size_t i = 0; i < d; i++) {
                double acc = 0.0;
                for (int j = 1; j <= q; j++) acc += cc[j] * yh[(size_t)(j - 1) * d + i];
                base[i] = -acc / cc[0];
            }
            lagrange_predict(t1, th, yh, q, d, pred);
            ok = nd_newton_theta(P, t1, base, hs, 1.0 / (cc[0] * hs), NULL, pred, ynext, tol);
        }

        if (!ok) {
            /* Newton diverged: drop order toward backward Euler, halve, retry. */
            no_grow = true;
            if (qtarget > 1) qtarget--;
            h = 0.5 * hs;
            if (fabs(h) < 16.0 * DBL_EPSILON * (fabs(t) + 1.0)) { st = ND_ERR_NONCONV; break; }
            continue;
        }

        /* Milne error estimate = corrector - predictor (same order q). */
        for (size_t i = 0; i < d; i++) base[i] = ynext[i] - pred[i];  /* base: scratch */
        double err = nd_wrms_norm(d, base, yh, ynext, tol);
        double qexp = (double)q + 1.0;

        if (err <= 1.0) {
            /* --- accept --- */
            t = t1;
            if (!nd_rhs_real(P, t, ynext, fnext)) { st = ND_ERR_SAMPLE; break; }
            nd_solution_push(sol, t, ynext, fnext);
            if (o->step_monitor) { Expr* mo = eval_and_free(expr_copy(o->step_monitor)); expr_free(mo); }

            /* Decide the next (order, step).  Between reassessments the step and
             * order are HELD constant for a run of q+1 steps: a uniform mesh makes
             * BDF zero-stable at every order <= 6, which is what lets high order be
             * used safely.  Only at the end of a hold window do we consider a
             * (possibly large) jump or an order change.  History is still t_n.. */
            int    qnew = q;
            double fac  = 1.0;
            if (no_grow) {
                fac = 1.0; hold = q + 1;        /* recovered from trouble: settle  */
            } else if (hold > 0) {
                fac = 1.0; hold--;              /* inside a hold window            */
            } else {
                double fac_same = pow(1.0 / (err > 1e-300 ? err : 1e-300), 1.0 / qexp);
                /* Lowering signal: if the order-(q-1) polynomial predictor implies
                 * a substantially larger step, the high-order differences are not
                 * decaying (advection / near-imaginary spectra where high-order
                 * BDF is unstable) — drop.  The 1.5 bias avoids spurious drops on
                 * smooth, locally-polynomial data (the exponent asymmetry). */
                if (q > 1) {
                    lagrange_predict(t1, th, yh, q - 1, d, fcur);   /* fcur scratch */
                    for (size_t i = 0; i < d; i++) fcur[i] = ynext[i] - fcur[i];
                    double err_dn = nd_wrms_norm(d, fcur, yh, ynext, tol);
                    double fac_dn = pow(1.0 / (err_dn > 1e-300 ? err_dn : 1e-300), 1.0 / (double)q);
                    if (fac_dn > 1.5 * fac_same) qnew = q - 1;
                }
                /* Raising: climb one order per reassessment; the lowering signal
                 * pulls it back where high order does not pay, so the order settles
                 * at the best value for the local spectrum. */
                if (qnew == q && q < BDF_QMAX && q <= m - 1) qnew = q + 1;
                fac = 0.9 * fac_same;
                double cap = bdf_grow_cap(qnew);
                if (fac < 0.2) fac = 0.2;
                if (fac > cap) fac = cap;
                hold = qnew + 1;                /* hold the new (order, step)      */
            }
            /* shift the history ring, record the new front node */
            int keep = (m < CAP) ? m : CAP - 1;
            for (int j = keep; j >= 1; j--) {
                th[j] = th[j - 1];
                memcpy(&yh[(size_t)j * d], &yh[(size_t)(j - 1) * d], sizeof(double) * d);
            }
            th[0] = t;
            memcpy(yh, ynext, sizeof(double) * d);
            if (m < CAP) m++;
            qtarget = qnew;
            no_grow = false;
            h = dir * fabs(hs) * fac;
        } else {
            /* reject: shrink, drop the order, force a reassessment next accept */
            no_grow = true; hold = 0;
            if (qtarget > 1) qtarget--;
            double fac = 0.9 * pow(1.0 / err, 1.0 / qexp);
            if (fac < 0.2) fac = 0.2;
            if (fac > 1.0) fac = 1.0;
            h = dir * fabs(hs) * fac;
        }
    }
done:
    free(yh); free(ynext); free(base); free(pred); free(fcur); free(fnext);
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
    int64_t budget = (o->max_steps > 0) ? o->max_steps : ND_AUTO_MAX_STEPS;
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
    "BDF", ND_IMPLICIT | ND_MULTISTEP, 5, 1, 0, NULL,
    "NDSolve`BDF[eqns, u, {x, xmin, xmax}]\n"
    "\tBackward differentiation formula, the stiff implicit multistep workhorse.\n"
    "\tVariable step size and variable order (1-5) with exact nonuniform-mesh\n"
    "\tcoefficients, predictor-corrector local error control, order ramping and\n"
    "\tNewton-failure step/order recovery; self-started at order 1."
};
