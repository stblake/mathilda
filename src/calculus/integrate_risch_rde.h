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

/* ---- Differential transcendental tower (shared with the integrator) ------- */

typedef enum { RT_LOG, RT_EXP, RT_TAN } RtKind; /* monomial kind: Log[u] / E^w / Tan[u] */

/* The per-monomial and member arrays are heap-allocated by rt_tower_build_min,
 * sized to the actual kernel count, so tower DEPTH IS UNBOUNDED (no RT_MAXK cap).
 * rt_tower_free releases them; a not-yet-built tower has all pointers NULL. */
typedef struct {
    size_t n;
    RtKind* kind;
    Expr** kernel;   /* Log[u_i] / Power[E, w_i] / Tan[u_i]      (owned) */
    Expr** arg;      /* u_i (log/tan argument) or w_i (exp exponent) (owned) */
    Expr** t;        /* fresh tower variable t_i                 (owned) */
    Expr** Dcoef;    /* log: u'/u ; exp: w' ; tan: u'           (owned) */
    long*  tsg;      /* RT_TAN special sign sigma: +1 (Tan/Cot), -1 (Tanh/Coth); 0 else */
    Expr* subrules;  /* List of all kernel -> t_i rules          (owned) */
    /* Multiplicatively commensurate non-primitive exp members: a collected
     * kernel E^w whose exponent w = mmult * arg[mprim] is NOT an independent
     * extension — it is (E^arg[mprim])^mmult = t[mprim]^mmult.  Recorded here so
     * rt_subst_kernels can alias it to a power of the primitive's tower var
     * instead of leaving it as a foreign kernel. */
    size_t nm;
    Expr** marg;     /* the member exponent w                    (owned) */
    long*  mprim;    /* tower index of the class primitive              */
    long*  mmult;    /* integer multiplier k (w = k * arg[mprim])       */
} RtTower;

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

/* ---- Shared helpers DEFINED in integrate_risch_transcendental.c ----------- */
/* (small build-and-evaluate / polynomial primitives reused by both files) */

bool  rt_head_is(const Expr* e, const char* name);
Expr* rt_eval1(const char* head, Expr* a);
Expr* rt_eval2(const char* head, Expr* a, Expr* b);
Expr* rt_eval3(const char* head, Expr* a, Expr* b, Expr* c);
bool  rt_is_poly(Expr* e, Expr* x);
long  rt_degree(Expr* e, Expr* x);
Expr* rt_coeff(Expr* e, Expr* x, long k);
bool  rt_is_zero(Expr* e);
Expr* rt_tower_deriv(Expr* e, RtTower* T, Expr* x);

/* Field integrator entry (Bronstein mutual recursion): integrate F in K_L over the
 * tower T.  Returns the antiderivative in tower-variable form (owned), or NULL; a
 * non-constant-residue remainder is returned via rem_out when non-NULL.  The RDE's
 * antidifferentiation branch (D h = c) calls this to solve for h in the field. */
Expr* rt_field_integrate(Expr* F, RtTower* T, long L, Expr* x, Expr** rem_out);

#endif /* MATHILDA_INTEGRATE_RISCH_RDE_H */
