/* Mathilda — NDSolve explicit Runge–Kutta steppers: classical RK4 and the
 * adaptive Dormand–Prince 5(4) pair (DOPRI5), the Automatic non-stiff default. */
#include "ndsolve_common.h"
#include "ndsolve_tableau.h"
#include <stdlib.h>

/* Classical RK4 (order 4, fixed step). */
static bool rk4_step(const NdStepper* S, NdProblem* P, double t, const double* Y,
                     double h, double* Ynew, double* Yerr, double* K) {
    (void)S; (void)Yerr; (void)K;
    size_t d = P->d;
    double* k1 = malloc(sizeof(double) * d);
    double* k2 = malloc(sizeof(double) * d);
    double* k3 = malloc(sizeof(double) * d);
    double* k4 = malloc(sizeof(double) * d);
    double* tmp = malloc(sizeof(double) * d);
    bool ok = nd_rhs_real(P, t, Y, k1);
    if (ok) { for (size_t i = 0; i < d; i++) tmp[i] = Y[i] + 0.5*h*k1[i];
              ok = nd_rhs_real(P, t + 0.5*h, tmp, k2); }
    if (ok) { for (size_t i = 0; i < d; i++) tmp[i] = Y[i] + 0.5*h*k2[i];
              ok = nd_rhs_real(P, t + 0.5*h, tmp, k3); }
    if (ok) { for (size_t i = 0; i < d; i++) tmp[i] = Y[i] + h*k3[i];
              ok = nd_rhs_real(P, t + h, tmp, k4); }
    if (ok) for (size_t i = 0; i < d; i++)
        Ynew[i] = Y[i] + (h/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    free(k1); free(k2); free(k3); free(k4); free(tmp);
    return ok;
}

/* Dormand–Prince 5(4): 7 stages, propagates the 5th-order solution with a
 * 4th-order embedded error estimate. */
static bool dopri5_step(const NdStepper* S, NdProblem* P, double t, const double* Y,
                        double h, double* Ynew, double* Yerr, double* K) {
    (void)S; (void)K;
    size_t d = P->d;
    double* k = malloc(sizeof(double) * 7 * d);
    double* tmp = malloc(sizeof(double) * d);
    bool ok = true;
    if (P->fsal) memcpy(&k[0], P->fsal_cur, sizeof(double) * d);  /* stage 1 reused (FSAL) */
    else ok = nd_rhs_real(P, t, Y, &k[0]);           /* stage 1 */
    for (int s = 1; s < 7 && ok; s++) {
        for (size_t i = 0; i < d; i++) {
            double acc = 0.0;
            for (int j = 0; j < s; j++) acc += DP_A[s][j] * k[j*d + i];
            tmp[i] = Y[i] + h * acc;
        }
        ok = nd_rhs_real(P, t + DP_C[s] * h, tmp, &k[s*d]);
    }
    if (ok) {
        for (size_t i = 0; i < d; i++) {
            double sol = 0.0, err = 0.0;
            for (int s = 0; s < 7; s++) { sol += DP_B[s] * k[s*d + i];
                                          err += DP_E[s] * k[s*d + i]; }
            Ynew[i] = Y[i] + h * sol;
            if (Yerr) Yerr[i] = h * err;
        }
        /* stage 7 == f(t+h, Ynew) (c7=1, row7=b): hand it to the driver as the
         * FSAL last stage (reused as node slope + next step's first stage). */
        if (P->fsal) memcpy(P->fsal_pending, &k[6*d], sizeof(double) * d);
    }
    free(k); free(tmp);
    return ok;
}

const NdStepper nd_stepper_rk4 = {
    "RK4", 0, 4, 0, 0, rk4_step,
    "NDSolve`RK4[eqns, u, {x, xmin, xmax}]\n"
    "\tClassical fourth-order Runge–Kutta method (fixed step)."
};

const NdStepper nd_stepper_dopri5 = {
    "DOPRI5", ND_ADAPTIVE | ND_FSAL, 5, 4, 7, dopri5_step,
    "NDSolve`DOPRI5[eqns, u, {x, xmin, xmax}]\n"
    "\tDormand–Prince 5(4) adaptive embedded Runge–Kutta pair; the default\n"
    "\tExplicitRungeKutta method used for non-stiff problems."
};
