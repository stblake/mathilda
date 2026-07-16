/* risch_tower.h — differential transcendental tower (Risch structure layer).
 *
 * The ordered differential tower `RtTower` over a single integration variable,
 * and the operational Risch structure-theorem machinery that builds and
 * manipulates it: collect the log / exp / tangent kernels of an integrand,
 * decide their commensurability and triangular independence, assemble the
 * tower (rt_tower_build_min), and provide the tower's structural operations
 * (lifecycle, monomial derivation, kernel substitution, generator-minimising
 * pre-normalisations).  Defined in risch_tower.c.
 *
 * This is DISTINCT from risch_structure.c, which implements the abstract
 * structure-theorem DECISION builtins (Risch`RationalSpan / LogReducible /
 * ExpReducible, Bronstein §9.3).  This file is the concrete tower that the
 * recursive integrator (integrate_risch_transcendental.c) and the Risch DE
 * (integrate_risch_rde.c) run over; integrate_risch_rde.h includes this header
 * for the shared `RtTower` type.
 *
 * A handful of small eval / predicate helpers are DEFINED in
 * integrate_risch_transcendental.c and shared back here (declared extern below),
 * so the files share one implementation — the same pattern as the RDE split.
 */

#ifndef MATHILDA_RISCH_TOWER_H
#define MATHILDA_RISCH_TOWER_H

#include <stdbool.h>
#include <stddef.h>
#include "expr.h"
#include "risch_util.h"   /* shared rt_ eval/predicate/poly helpers */

/* ---- Differential transcendental tower ------------------------------------ */

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

/* ---- Tower construction & structural operations (defined in risch_tower.c) - */

/* Build the ordered differential tower of f over C(x): collect every x-dependent
 * kernel, order innermost-first by structural containment, assign tower variables
 * and derivation coefficients, and verify the structure-theorem triangularity.
 * Returns true with T populated (n >= min_n), false otherwise; the caller must
 * still call rt_tower_free either way. */
bool  rt_tower_build_min(Expr* f, Expr* x, RtTower* T, size_t min_n);
void  rt_tower_free(RtTower* T);

/* Monomial derivation of the tower: D_tower[e] = D[e,x] + sum_i Dt_i D[e,t_i]. */
Expr* rt_tower_deriv(Expr* e, RtTower* T, Expr* x);
Expr* rt_dt_i(RtTower* T, size_t i);              /* Dt_i (log/exp/tan), owned */
Expr* rt_build_deriv_rules(RtTower* T, Expr* x);  /* List[Rule[var, Dvar], ...] */

/* Structural (non-evaluating) substitution of tower kernels -> tower variables. */
Expr* rt_subst_kernels(Expr* e, RtTower* T);

/* Generator-minimising / kernel-exposing pre-normalisations. */
Expr* rt_powers_to_exp(Expr* e, Expr* x);   /* b^e -> Power[E, e Log b] (transcendental) */
Expr* rt_expand_logs(Expr* e);              /* Log[a b] -> Log a + Log b, Log[b^p] -> p Log b */
Expr* rt_log_of(Expr* a);                   /* expanded Log[a] (argument borrowed) */
Expr* rt_expand_exp_sums(Expr* e);          /* split merged E^(a+b) isolating nested kernels */
bool  rt_has_explog_kernel(Expr* e);        /* does e contain any Exp/E^/Log kernel? */

/* Kernel collectors (owned, deduplicated). */
void  rt_collect_exp_exponents(Expr* e, Expr* x, Expr*** arr, size_t* n, size_t* cap);
void  rt_collect_logs(Expr* e, Expr* x, Expr*** arr, size_t* n, size_t* cap);
void  rt_collect_tangents(Expr* e, Expr* x, Expr*** args, long** sigs, size_t* n, size_t* cap);

/* Structural helpers. */
bool  rt_contains(Expr* big, Expr* small);                 /* subexpression containment */
Expr* rt_build_monomial(Expr** lv, const long* e, size_t nlv);  /* prod lv[j]^e[j] */
void  rt_decode_mono(long idx, const long* bd, size_t nlv, long* e);  /* flat idx -> exps */

/* Small shared eval/predicate/poly helpers (rt_eval*, rt_is_rat_const,
 * rt_free_of_x, rt_class_primitive, rt_find_log_of_x, ...) are declared in
 * risch_util.h, included above. */

#endif /* MATHILDA_RISCH_TOWER_H */
