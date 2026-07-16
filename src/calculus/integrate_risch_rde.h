/* integrate_risch_rde.h — Bronstein Chapter 6 Risch differential equation.
 *
 * The Risch DE solver  D[y] + f y = g,  factored out of the transcendental
 * integrator (integrate_risch_transcendental.c) into its own module.  Two
 * layers:
 *
 *   - the BASE field C(x), D = d/dx (rde_base and the classical §6.1–6.5 stack:
 *     WeakNormalizer, RdeNormalDenominator, RdeBoundDegreeBase, SPDE,
 *     PolyRischDENoCancel), an audited, complete decision procedure for that
 *     field;
 *   - the recursive TOWER lift (Gap 1): the same Bronstein boxes over an
 *     arbitrary K_m = C(x)(t_0..t_m), parameterised by the RdeCtx derivation
 *     abstraction, for the primitive (RT_LOG) top, non-cancellation branch
 *     (rde_tower and the *_field helpers).
 *
 * The module reuses a handful of small tower/eval helpers that remain DEFINED in
 * integrate_risch_transcendental.c (declared extern below) so the two files
 * share one implementation; and it exposes the RtTower / RdeCtx / RdeSpde types
 * plus the solver entry points that the integrator and the Risch` builtins call.
 */

#ifndef MATHILDA_INTEGRATE_RISCH_RDE_H
#define MATHILDA_INTEGRATE_RISCH_RDE_H

#include <stdbool.h>
#include <stddef.h>
#include "expr.h"
#include "risch_util.h"    /* shared rt_ eval/predicate/poly helpers */
#include "risch_tower.h"   /* RtKind / RtTower (the differential tower) */

/* Derivation context for the recursive Risch DE: the current top monomial
 * tau = t_m (kind T->kind[m]) with the lower tower vars {x, t_0..t_{m-1}} = k as
 * a transcendental coefficient field.  m < 0 selects the base field C(x). */
typedef struct {
    RtTower* T;
    long     m;      /* current top monomial index; m < 0 => base field C(x)   */
    Expr*    x;
    Expr*    tau;    /* polynomial variable: T->t[m] (m>=0) or x (m<0), borrowed */
} RdeCtx;

/* SPDE reduction result: a solution q of a Dq + b q = c with deg(q) <= n has the
 * form q = alpha*H + beta, where deg(H) <= m and D H + b H = c. */
typedef struct { Expr* b; Expr* c; long m; Expr* alpha; Expr* beta; } RdeSpde;
void rde_spde_free(RdeSpde* s);

/* ---- Solver entry points (defined in integrate_risch_rde.c) --------------- */

/* Base field C(x): solve D y + f y = g for a rational y (Bronstein Ch.6).
 * Returns y (owned), 0 when g == 0, or NULL for the authoritative "no solution". */
Expr* rde_base(Expr* f, Expr* g, Expr* x);

/* Solve q' + i u' q = p for q in C(x) (f = i u'); NULL == no solution. */
Expr* rt_solve_rde(Expr* p, long i, Expr* u, Expr* x);

/* Recursive Risch DE over the tower: solve D_tower[y] + f y = g for y in K_m
 * (C->m).  m < 0 delegates to rde_base; m >= 0 with a primitive (RT_LOG) top
 * runs the literal Bronstein reduction (Gap 1a), certified by the exact tower
 * identity.  NULL == no solution / out of the current increment's scope. */
Expr* rde_tower(Expr* f, Expr* g, RdeCtx* C);

/* Bronstein SPDE box (§6.4) over k(x)[tau]: 1 + *out on success, 0 for no
 * solution of degree <= n. */
int rde_spde_field(Expr* a, Expr* b, Expr* c, RdeCtx* C, long n, RdeSpde* out);

/* PolyRischDENoCancel over k(x)[tau] (§6.5, deg_tau(b) >= 1): q (owned) or NULL. */
Expr* rde_polyrischde_nocancel_field(Expr* b, Expr* c, RdeCtx* C, long n);

/* Small shared eval/predicate/poly helpers (rt_head_is, rt_eval1..3,
 * rt_is_poly, rt_degree, rt_coeff, rt_is_zero, ...) are declared in
 * risch_util.h; rt_tower_deriv in risch_tower.h — both included above. */

/* Field integrator entry (Bronstein mutual recursion): integrate F in K_L over the
 * tower T.  Returns the antiderivative in tower-variable form (owned), or NULL; a
 * non-constant-residue remainder is returned via rem_out when non-NULL.  The RDE's
 * antidifferentiation branch (D h = c) calls this to solve for h in the field. */
Expr* rt_field_integrate(Expr* F, RtTower* T, long L, Expr* x, Expr** rem_out);

#endif /* MATHILDA_INTEGRATE_RISCH_RDE_H */
