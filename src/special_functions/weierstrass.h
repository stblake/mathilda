/* Mathilda -- the Weierstrass elliptic functions WeierstrassP, WeierstrassPPrime,
 * WeierstrassZeta, WeierstrassSigma [u, {g2, g3}], WeierstrassHalfPeriods,
 * WeierstrassInvariants and InverseWeierstrassP. See weierstrass.c. */
#ifndef MATHILDA_WEIERSTRASS_H
#define MATHILDA_WEIERSTRASS_H

#include <stdbool.h>
#include "expr.h"

/* Function codes of the four lattice functions. */
enum { WP_P = 0, WP_PPRIME = 1, WP_ZETA = 2, WP_SIGMA = 3 };

/* The lattice of invariants (g2, g3) in machine arithmetic: the half-periods
 * in the reference system's normalisation (tau = w3/w1 in the closed
 * fundamental domain, its boundary and sign rules) and the theta constants the
 * function kernels need. kind: 0 generic, 1 g2 = g3 = 0, 2 a double root
 * (g2^3 = 27 g3^2, g2 != 0). */
typedef struct {
    int    kind;
    bool   real;                 /* g2, g3 real */
    double g2r, g2i, g3r, g3i;
    double w1r, w1i, w3r, w3i;   /* the half-periods (kind 0) */
    double tr, ti, trl, til;     /* tau = w3 / w1, double-double */
    double qr, qi, q4r, q4i;     /* nome exp(i pi tau), its root exp(i pi tau / 4) */
    double i2r, i2i, i2rl, i2il; /* 1 / (2 w1), double-double */
    double pAr, pAi;             /* pi^2 theta2^2 theta3^2 / (2 w1)^2 (theta at 0) */
    double pBr, pBi;             /* pi^2 (theta2^4 + theta3^4) / 3 / (2 w1)^2 */
    double kpr, kpi;             /* -2 pi^3 theta1'(0)^2 / (2 w1)^3 */
    double scr, sci;             /* 2 w1 / (pi theta1'(0)), theta1'(0) = th2 th3 th4 */
    double etr, eti, etrl, etil; /* zeta(1/2) of the lattice (1, tau): pi^2 E2(tau) / 6 */
    double er, ei, cr, ci;       /* kind 2: the double root e and c = sqrt(3 e) */
} WpPlan;

/* Build the plan for (g2, g3); false when the machine route cannot certify the
 * lattice (it then answers through arb). */
bool sf_wp_plan(double g2r, double g2i, double g3r, double g3i, WpPlan* P);

/* The machine kernel for function `kind` at u over invariants (g2, g3): the
 * helper behind the scalar builtins, the NDArray / packed map and the Compile[]
 * VM. A value the double route cannot certify (cancellation near a zero, an
 * argument near a pole, a far-out sigma) is answered by arb at 53 bits. Declines
 * only at a pole or outside the double range. */
bool sf_machine_weierstrass(int kind, double g2r, double g2i, double g3r, double g3i,
                            double ur, double ui, double* rr, double* ri);

/* The plan-level evaluator the buffer path uses after building the plan once. */
bool sf_wp_eval_plan(const WpPlan* P, int kind, double ur, double ui, double* rr, double* ri);

/* Normalise a lattice basis the way the reference system does: Gauss-reduce
 * tau = w3/w1 into the fundamental domain, then apply its boundary rules
 * (Re tau = +1/2 rather than -1/2, Re tau >= 0 on |tau| = 1), the rotation
 * rules of the square (sq, g3 = 0) and hexagonal (hex, g2 = 0) lattices, and
 * Re w1 > 0 (Im w1 < 0 on the imaginary axis). M receives the integer matrix
 * with (w1', w3') = (M0 w1 + M1 w3, M2 w1 + M3 w3). False on a degenerate
 * basis. Shared by the machine plan and the arb route. */
bool sf_wp_normalize(double w1r, double w1i, double w3r, double w3i, bool sq, bool hex,
                     double M[4]);

void weierstrass_init(void);

#endif /* MATHILDA_WEIERSTRASS_H */
