/* Mathilda — NDSolve explicit fixed-step submethods: ExplicitEuler, ExplicitMidpoint. */
#include "ndsolve_common.h"
#include <stdlib.h>

/* Forward Euler: Ynew = Y + h f(t, Y). Order 1. */
static bool euler_step(const NdStepper* S, NdProblem* P, double t, const double* Y,
                       double h, double* Ynew, double* Yerr, double* K) {
    (void)S; (void)Yerr; (void)K;
    size_t d = P->d;
    double* k = malloc(sizeof(double) * d);
    bool ok = nd_rhs_real(P, t, Y, k);
    if (ok) for (size_t i = 0; i < d; i++) Ynew[i] = Y[i] + h * k[i];
    free(k);
    return ok;
}

/* Explicit midpoint (order 2): k1=f(t,Y); k2=f(t+h/2, Y+h/2 k1); Ynew=Y+h k2. */
static bool midpoint_step(const NdStepper* S, NdProblem* P, double t, const double* Y,
                          double h, double* Ynew, double* Yerr, double* K) {
    (void)S; (void)Yerr; (void)K;
    size_t d = P->d;
    double* k1 = malloc(sizeof(double) * d);
    double* k2 = malloc(sizeof(double) * d);
    double* tmp = malloc(sizeof(double) * d);
    bool ok = nd_rhs_real(P, t, Y, k1);
    if (ok) {
        for (size_t i = 0; i < d; i++) tmp[i] = Y[i] + 0.5 * h * k1[i];
        ok = nd_rhs_real(P, t + 0.5 * h, tmp, k2);
    }
    if (ok) for (size_t i = 0; i < d; i++) Ynew[i] = Y[i] + h * k2[i];
    free(k1); free(k2); free(tmp);
    return ok;
}

const NdStepper nd_stepper_explicit_euler = {
    "ExplicitEuler", 0, 1, 0, 0, euler_step,
    "NDSolve`ExplicitEuler[eqns, u, {x, xmin, xmax}]\n"
    "\tForward (explicit) Euler method, order 1.  A fixed-step submethod; the\n"
    "\tdriver adapts the step by step-doubling to meet the accuracy goals."
};

const NdStepper nd_stepper_explicit_midpoint = {
    "ExplicitMidpoint", 0, 2, 0, 0, midpoint_step,
    "NDSolve`ExplicitMidpoint[eqns, u, {x, xmin, xmax}]\n"
    "\tExplicit midpoint (modified Euler) method, order 2."
};
