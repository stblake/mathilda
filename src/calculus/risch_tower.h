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

/* monomial kind: Log[u] / E^w / Tan[u] / non-elementary primitive SF[u]
 * (Knowles Liouvillian tower — RT_PRIM = a primitive theta with theta' = Dcoef,
 * e.g. Erf/Erfi/Erfc/ExpIntegralEi/LogIntegral; see KNOWLES_DESIGN.md §2.1). */
typedef enum { RT_LOG, RT_EXP, RT_TAN, RT_PRIM } RtKind;

/* The per-monomial and member arrays are heap-allocated by rt_tower_build_min,
 * sized to the actual kernel count, so tower DEPTH IS UNBOUNDED (no RT_MAXK cap).
 * rt_tower_free releases them; a not-yet-built tower has all pointers NULL. */
typedef struct {
    size_t n;
    RtKind* kind;
    Expr** kernel;   /* Log[u_i] / Power[E, w_i] / Tan[u_i]      (owned) */
    Expr** arg;      /* u_i (log/tan argument) or w_i (exp exponent) (owned) */
    Expr** t;        /* fresh tower variable t_i                 (owned) */
    Expr** Dcoef;    /* log: u'/u ; exp: w' ; tan: u' ; prim: theta' (owned) */
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
    /* The member's value in tower variables (owned): t[mprim]^mmult, or, for a
     * member of a class whose primitive is an INTEGER combination of other
     * primitives (E^(2u + u^2) beside E^(2u) and E^(u^2): review F12), the
     * monomial Prod t_s^(k_s).  rt_subst_kernels returns it. */
    Expr** mval;
    /* Dependent logarithm members (campaign step 8): a collected Log[v] whose
     * u'/u is in the Q-span of the lower Log/Exp levels is NOT a level.  It is
     * Log[v] = Sum r_j t_j + Sum s_j w_j + c exactly, c an exact constant
     * (Log of the constant ratio, with the branch multiple of 2 Pi I read off
     * sample points), and rt_subst_kernels returns that value. */
    size_t nlm;
    Expr** lmk;      /* the member kernel Log[v]                (owned) */
    Expr** lmv;      /* its value in tower variables            (owned) */
    /* True when the Log/Exp monomials are NOT algebraically independent (the
     * structure theorem finds Dcoef_i in the Q-span of the earlier ones, e.g.
     * E^(x + x^2) beside E^x and E^(x^2)), or when that could not be decided.
     * Integration over such a tower still verifies its answers, but every
     * "no solution" / non-elementarity theorem presupposes independence, so a
     * certificate read off a dependent tower is demoted to undecided. */
    bool dependent;
} RtTower;

/* ---- Tower construction & structural operations (defined in risch_tower.c) - */

/* The constant e^i as a symbol (rt_dependent_exp_rule): a Gaussian-integer
 * offset c = a + b i of a dependent exponential gives E^a Z^b with Z this
 * symbol.  E and E^I are algebraically independent (Lindemann-Weierstrass:
 * 1 and i are linearly independent over Q), so Z is a parameter. */
#define RT_EI_SYMBOL "Risch`Private`rtEI"
/* The back-substitution of a tower answer (owned List): t_i -> kernel_i for
 * every level, and RT_EI_SYMBOL -> E^I. */
Expr* rt_tower_back_rules(const RtTower* T);

/* Build the ordered differential tower of f over C(x): collect every x-dependent
 * kernel, order innermost-first by structural containment, assign tower variables
 * and derivation coefficients, and verify the structure-theorem triangularity.
 * Returns true with T populated (n >= min_n), false otherwise; the caller must
 * still call rt_tower_free either way. */
bool  rt_tower_build_min(Expr* f, Expr* x, RtTower* T, size_t min_n);
/* Both eliminate a dependent logarithm -- it becomes a member
 * Log[v] = Sum r_j t_j + c with c exact, not a level -- and rebuild, clearing
 * T->dependent where they can.  rt_tower_build_min only when c holds on a wide
 * sample set (both half-planes, negative reals); rt_tower_build_decide also
 * for a constant of the positive reals (the decision's germ). */
bool  rt_tower_build_decide(Expr* f, Expr* x, RtTower* T, size_t min_n);
void  rt_tower_free(RtTower* T);

/* Monomial derivation of the tower: D_tower[e] = D[e,x] + sum_i Dt_i D[e,t_i]. */
Expr* rt_tower_deriv(Expr* e, RtTower* T, Expr* x);
Expr* rt_dt_i(RtTower* T, size_t i);              /* Dt_i (log/exp/tan), owned */
Expr* rt_build_deriv_rules(RtTower* T, Expr* x);  /* List[Rule[var, Dvar], ...] */

/* The Moebius change at a split tangent level m (Tanh, alpha = 1; Tan over
 * k(i), alpha = I; Bronstein, Symbolic Integration I, p. 166):
 * theta = (alpha + t)/(alpha - t) has D theta/theta = 2 alpha Dt/(alpha^2 - t^2)
 * = a in K_{m-1}, so k(t) = k(theta) with theta an exponential monomial.
 * Fills T2 with levels 0..m of T, level m replaced by theta (the lower levels
 * and the member tables borrowed from T; subrules rewritten so the tangent
 * kernel maps to alpha (theta - 1)/(theta + 1)).  *to_theta is the Rule
 * t -> alpha (theta - 1)/(theta + 1), *from_theta theta -> (alpha + t)/(alpha - t)
 * (both owned).  False (T2 untouched) when a is not free of t.  Release T2
 * with rt_tower_mobius_free, never rt_tower_free. */
bool  rt_tower_mobius(RtTower* T, size_t m, Expr* alpha, RtTower* T2,
                      Expr** to_theta, Expr** from_theta);
void  rt_tower_mobius_free(RtTower* T2);

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
void  rt_collect_primitives(Expr* e, Expr* x, Expr*** arr, size_t* n, size_t* cap);
bool  rt_is_primitive_head(const char* h);   /* Erf/Erfi/Erfc/ExpIntegralEi/LogIntegral */

/* Structural helpers. */
bool  rt_contains(Expr* big, Expr* small);                 /* subexpression containment */
Expr* rt_build_monomial(Expr** lv, const long* e, size_t nlv);  /* prod lv[j]^e[j] */
void  rt_decode_mono(long idx, const long* bd, size_t nlv, long* e);  /* flat idx -> exps */

/* Small shared eval/predicate/poly helpers (rt_eval*, rt_is_rat_const,
 * rt_free_of_x, rt_class_primitive, rt_find_log_of_x, ...) are declared in
 * risch_util.h, included above. */

#endif /* MATHILDA_RISCH_TOWER_H */
