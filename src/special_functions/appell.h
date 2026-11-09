/* Mathilda -- the Appell hypergeometric functions of two variables:
 *
 *   AppellF1[a, b1, b2, c, x, y]        = Sum (a)_(m+n) (b1)_m (b2)_n / ((c)_(m+n) m! n!) x^m y^n
 *   AppellF2[a, b1, b2, c1, c2, x, y]   = Sum (a)_(m+n) (b1)_m (b2)_n / ((c1)_m (c2)_n m! n!) x^m y^n
 *   AppellF3[a1, a2, b1, b2, c, x, y]   = Sum (a1)_m (a2)_n (b1)_m (b2)_n / ((c)_(m+n) m! n!) x^m y^n
 *   AppellF4[a, b, c1, c2, x, y]        = Sum (a)_(m+n) (b)_(m+n) / ((c1)_m (c2)_n m! n!) x^m y^n
 *
 * See appell.c for the method. Every numeric surface -- the scalar builtin,
 * the NDArray / packed n-ary map and the Compile[] VM -- goes through
 * sf_machine_appell at machine precision, so the four agree bit for bit; the
 * arbitrary-precision evaluator is appell_acb (appell_arb.c), on the shared
 * precision ladder of flint_num_bridge.c.
 */
#ifndef MATHILDA_APPELL_H
#define MATHILDA_APPELL_H

#include <stdbool.h>
#include "expr.h"

/* kind 1..4 selects F1..F4; the argument count is 6 (F1, F4) or 7 (F2, F3). */
#define APPELL_NARGS(kind) (((kind) == 2 || (kind) == 3) ? 7 : 6)

/* The machine kernel: re[] / im[] hold the arguments in the head's own order.
 * A value the double-precision route cannot certify (cancellation, a region
 * without a convergent expansion, complex parameters beyond the machine 2F1)
 * is answered by arb at 53 bits inside the kernel. False only when there is no
 * finite value (a pole: c a non-positive integer; an unrepresentable result)
 * or no method at all. */
bool sf_machine_appell(int kind, const double* re, const double* im, double* rr, double* ri);

/* ------------------------------------------------------------------ */
/* The expansion plan, shared by the machine and the arb evaluators     */
/* ------------------------------------------------------------------ */

/* Every evaluation is one OUTER series over m whose terms carry an inner Gauss
 * function of the other variable, continued analytically:
 *
 *   S = Sum_m T_m G_m,  T_(m+1)/T_m = (n1 + m)(n2 + m) X / ((d1 + m)(m + 1)),
 *   G_m = 2F1(A0 + sA m, B0 + sB m; C0 + sC m; Y),
 *
 * after a transformation (F1: the six Euler-integral substitutions; F2: the
 * three Pfaff-type ones) that puts a small variable outside. The plan records
 * which transformation, which variable is outer, and the side of the cut the
 * inner variable sits on when it is real and > 1 (inputs on a cut are taken as
 * the limit from below, x - i0, the convention of Hypergeometric2F1). */
typedef struct {
    int    kind;   /* 1..4 */
    int    tr;     /* transformation index (appell.c, AP_TR) */
    int    swap;   /* the second variable is the outer one */
    int    sigY;   /* -1 / +1: the inner variable is on its cut, from below / above */
    double r;      /* modulus of the outer series' ratio (0: it terminates) */
} AppellPlan;

/* Choose the expansion for (kind, args). False when no transformation gives a
 * series with ratio <= rmax (the caller then has the Euler integral, F1). */
bool appell_plan(int kind, const double* re, const double* im, double rmax, AppellPlan* P);

/* The transformation table, read by both evaluators. */
#define AP_MAXP 5          /* parameters before (x, y) */
#define AP_MAXB 2          /* power factors in a prefactor */
enum { AP_V_X, AP_V_Y, AP_V_XM, AP_V_YM, AP_V_YX1, AP_V_XY1, AP_V_YXY, AP_V_XYX,
       AP_V_Y1X, AP_V_X1Y, AP_V_XS, AP_V_YS };
enum { AP_B_1X, AP_B_1Y, AP_B_1XY };
typedef struct {
    int nbase;
    int base[AP_MAXB];                  /* AP_B_*: the base of each power factor */
    signed char expo[AP_MAXB][AP_MAXP]; /* its exponent, a combination of the parameters */
    signed char newp[AP_MAXP][AP_MAXP]; /* the new parameters, combinations of the old */
    int vx, vy;                         /* AP_V_*: the new variables */
    int rule;                           /* validity rule: 0 always, 1 / 2 the sector tests */
} AppellTr;

/* The transformation `tr` of `kind` (NULL past the end). */
const AppellTr* appell_tr(int kind, int tr);

/* The engine's (n1, n2, d1, A0, B0, C0) as indices into the transformed
 * parameters, for `kind` with the second variable outer when `swap`. */
void appell_engine_map(int kind, int swap, int idx[6]);

/* The inner functions' growth rate |G_(m+1)/G_m| and whether their recurrence
 * runs backward (G the minimal solution) at inner variable Y. */
void appell_engine_info(int kind, double yr, double yi, double* growth, int* backward);

/* Number of parameters before (x, y): 4 or 5. */
int appell_np(int kind);

void appell_init(void);

#ifdef USE_FLINT
#include <flint/acb.h>
/* The Appell function at precision `prec`: p[] holds the arguments in the
 * head's order. Sets `res` to an indeterminate ball when no method applies. */
void appell_acb(acb_t res, int kind, acb_srcptr p, slong prec);
#endif

#endif /* MATHILDA_APPELL_H */
