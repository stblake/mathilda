/* Mathilda — Butcher tableaux for the explicit Runge–Kutta steppers.
 * Coefficients are exact rationals evaluated to double at load. */
#ifndef NDSOLVE_TABLEAU_H
#define NDSOLVE_TABLEAU_H

/* Classical RK4 (order 4). */
static const double RK4_C[4] = { 0.0, 0.5, 0.5, 1.0 };
static const double RK4_B[4] = { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };

/* Dormand–Prince 5(4), 7 stages, FSAL.  A is strictly lower triangular. */
static const double DP_C[7] = { 0.0, 1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0 };
static const double DP_A[7][6] = {
    { 0, 0, 0, 0, 0, 0 },
    { 1.0/5.0, 0, 0, 0, 0, 0 },
    { 3.0/40.0, 9.0/40.0, 0, 0, 0, 0 },
    { 44.0/45.0, -56.0/15.0, 32.0/9.0, 0, 0, 0 },
    { 19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0, 0, 0 },
    { 9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0, 0 },
    { 35.0/384.0, 0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0 }
};
/* 5th-order solution weights (row 7 of A; b7 = 0). */
static const double DP_B[7] = { 35.0/384.0, 0, 500.0/1113.0, 125.0/192.0,
                                -2187.0/6784.0, 11.0/84.0, 0.0 };
/* e = b - b*  (difference to the embedded 4th-order weights). */
static const double DP_E[7] = { 71.0/57600.0, 0, -71.0/16695.0, 71.0/1920.0,
                                -17253.0/339200.0, 22.0/525.0, -1.0/40.0 };

#endif /* NDSOLVE_TABLEAU_H */
