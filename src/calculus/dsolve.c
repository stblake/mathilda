/*
 * dsolve.c — DSolve dispatcher (cascade polyalgorithm).
 *
 * Mirrors src/calculus/integrate.c: a Method-option enum selects either the
 * automatic cascade (try each method until one returns a non-NULL result) or a
 * single pinned method (strict, no fallback).  A per-command fail-memo keyed on
 * eval_toplevel_id() suppresses the fixed-point loop's redundant re-entry on a
 * problem the deterministic cascade already declined, and g_dsolve_depth
 * distinguishes the outermost user call from internal recursions.
 *
 * The shared problem substrate (parse / verify / fit / assemble) is in
 * dsolve_common.c; each method is one file src/calculus/dsolve_<method>.c.
 */
#include "dsolve.h"
#include "dsolve_common.h"
#include "../core.h"
#include "dsolve_linsys.h"

#include "../sym_names.h"
#include "../sym_intern.h"
#include "../eval.h"
#include "../symtab.h"
#include "../attr.h"
#include "../expr.h"
#include "../arithmetic.h"
#include "../message.h"         /* nested sub-solves run as a muted probe region */

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ---- method table (extended per milestone) ---- */
typedef enum {
    DS_AUTOMATIC = 0,
    DS_FACTORABLE,
    DS_NTHALGEBRAIC,
    DS_QUADRATURE,
    DS_LINEAR1,
    DS_BERNOULLI,
    DS_HOMOGENEOUS,
    DS_SEPARABLE,
    DS_EXACT,
    DS_CLAIRAUT,
    DS_LAGRANGE,
    DS_UNDETCOEFF,
    DS_CONSTCOEFF,
    DS_EULER,
    DS_EXACTODE,
    DS_SPECIALFORM,
    DS_SFTRANSFORM,
    DS_KOVACIC,
    DS_GENAIRY,
    DS_OPERFACTOR,
    DS_SYMSQUARE,
    DS_REDUCEORDER,
    DS_LOWERREDUCE,
    DS_FOS,
    DS_RICCATI,
    DS_CHINI,
    DS_ABEL,
    DS_ABEL_AIR,
    DS_POLYSHIFT,
    DS_LINCOEFF,
    DS_ALMOSTLINEAR,
    DS_SEPREDUCED,
    DS_LINEARIZABLE,
    DS_LIE,
    DS_AUTONOMOUS,
    DS_SCALING2,
    DS_LIOUVILLE,
    DS_IFACTOR,
    DS_FIRSTINTEGRAL,
    DS_FOPOWERSERIES,
    DS_FROBENIUS,
    DS_SOLVABLEFORY,
    DS_SOLVABLEFORX,
    DS_PARAMQUAD,
    DS_INVALID
} DSolveMethod;

static DSolveMethod ds_method_from_string(const char* s) {
    if (strcmp(s, "Automatic")        == 0) return DS_AUTOMATIC;
    if (strcmp(s, "Factorable")       == 0) return DS_FACTORABLE;
    if (strcmp(s, "NthAlgebraic")     == 0) return DS_NTHALGEBRAIC;
    if (strcmp(s, "Quadrature")       == 0) return DS_QUADRATURE;
    if (strcmp(s, "LinearFirstOrder") == 0) return DS_LINEAR1;
    if (strcmp(s, "Bernoulli")        == 0) return DS_BERNOULLI;
    if (strcmp(s, "Homogeneous")      == 0) return DS_HOMOGENEOUS;
    if (strcmp(s, "Separable")        == 0) return DS_SEPARABLE;
    if (strcmp(s, "Exact")            == 0) return DS_EXACT;
    if (strcmp(s, "Clairaut")         == 0) return DS_CLAIRAUT;
    if (strcmp(s, "Lagrange")         == 0) return DS_LAGRANGE;
    if (strcmp(s, "SolvableForY")     == 0) return DS_SOLVABLEFORY;
    if (strcmp(s, "SolvableForX")     == 0) return DS_SOLVABLEFORX;
    if (strcmp(s, "ParametricQuadrature") == 0) return DS_PARAMQUAD;
    if (strcmp(s, "UndeterminedCoefficients") == 0) return DS_UNDETCOEFF;
    if (strcmp(s, "LinearConstantCoefficients") == 0) return DS_CONSTCOEFF;
    if (strcmp(s, "EulerCauchy")      == 0) return DS_EULER;
    if (strcmp(s, "ExactODE")         == 0) return DS_EXACTODE;
    if (strcmp(s, "SpecialFunctionForm") == 0) return DS_SPECIALFORM;
    if (strcmp(s, "SpecialFunctionTransform") == 0) return DS_SFTRANSFORM;
    if (strcmp(s, "Kovacic")            == 0) return DS_KOVACIC;
    if (strcmp(s, "GeneralizedAiry")    == 0) return DS_GENAIRY;
    if (strcmp(s, "OperatorFactor")     == 0) return DS_OPERFACTOR;
    if (strcmp(s, "SymmetricSquare") == 0) return DS_SYMSQUARE;
    if (strcmp(s, "ReductionOfOrder") == 0) return DS_REDUCEORDER;
    if (strcmp(s, "LowerDerivativeReduction") == 0) return DS_LOWERREDUCE;
    if (strcmp(s, "FirstOrderSubstitution") == 0) return DS_FOS;
    if (strcmp(s, "Riccati")              == 0) return DS_RICCATI;
    if (strcmp(s, "Chini")                == 0) return DS_CHINI;
    if (strcmp(s, "Abel")                 == 0) return DS_ABEL;
    if (strcmp(s, "AbelAIR")              == 0) return DS_ABEL_AIR;
    if (strcmp(s, "PolynomialShiftSubstitution") == 0) return DS_POLYSHIFT;
    if (strcmp(s, "LinearCoefficients")   == 0) return DS_LINCOEFF;
    if (strcmp(s, "AlmostLinear")         == 0) return DS_ALMOSTLINEAR;
    if (strcmp(s, "SeparableReduced")     == 0) return DS_SEPREDUCED;
    if (strcmp(s, "Linearizable")         == 0) return DS_LINEARIZABLE;
    if (strcmp(s, "LieSymmetry")          == 0) return DS_LIE;
    if (strcmp(s, "LieGroup")             == 0) return DS_LIE;
    if (strcmp(s, "AutonomousReduction") == 0) return DS_AUTONOMOUS;
    if (strcmp(s, "ScalingReduction") == 0) return DS_SCALING2;
    if (strcmp(s, "Liouville")            == 0) return DS_LIOUVILLE;
    if (strcmp(s, "ReducibleIntegratingFactor") == 0) return DS_IFACTOR;
    if (strcmp(s, "ReducibleFirstIntegral") == 0) return DS_FIRSTINTEGRAL;
    if (strcmp(s, "FirstOrderPowerSeries") == 0) return DS_FOPOWERSERIES;
    if (strcmp(s, "FrobeniusSeries")     == 0) return DS_FROBENIUS;
    if (strcmp(s, "PowerSeries")         == 0) return DS_FROBENIUS;
    return DS_INVALID;
}

/* method try-functions + registrars (one file each) */
extern Expr** dsolve_factorable_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_nth_algebraic_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_nth_algebraic_implicit_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_nth_algebraic_binomial_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_quadrature_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_linear1_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_bernoulli_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_bernoulli_expanded_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_homogeneous_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_homogeneous_implicit_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_separable_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_separable_implicit_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_exact_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_exact_implicit_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_exact_mu_u_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_exact_mu_u_implicit_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_clairaut_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_lagrange_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_solvablefory_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_solvablefory_implicit_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_homog_g_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_homog_g_implicit_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_solvableforx_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_solvableforx_linear_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_paramquad_try(DSolveProblem* P, size_t* nbranch);
extern void dsolve_paramquad_init(void);
extern Expr** dsolve_undetcoeff_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_constcoeff_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_euler_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_exactode_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_exactode_nl_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_specialform_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_specialform_transform_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_kovacic_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_genairy_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_operfactor_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_symsquare_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_reduce_order_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_lower_reduce_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_fos_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_fos_implicit_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_riccati_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_chini_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_abel_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_abel_air_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_polyshift_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_lincoeff_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_lincoeff_implicit_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_almostlinear_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_sepreduced_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_linearizable_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_lie_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_inverse_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_knownsol_try(DSolveProblem* P, size_t* nbranch);
extern void dsolve_knownsol_init(void);
extern void dsolve_inverse_init(void);
extern void dsolve_homog_g_init(void);
extern Expr** dsolve_autonomous_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_autonomous_implicit_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_autonomous_elliptic_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_liouville_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_lie2_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_ifactor_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_ifactor_first_integral_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_changevar_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_ratsol2_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_nonhomog_vop_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_nonhomog_vop_inert_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_piecewise_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_frobenius_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_first_order_series_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_frobenius_shifted_try(DSolveProblem* P, size_t* nbranch);
extern void dsolve_factorable_init(void);
extern void dsolve_nth_algebraic_init(void);
extern void dsolve_quadrature_init(void);
extern void dsolve_linear1_init(void);
extern void dsolve_bernoulli_init(void);
extern void dsolve_homogeneous_init(void);
extern void dsolve_separable_init(void);
extern void dsolve_exact_init(void);
extern void dsolve_clairaut_init(void);
extern void dsolve_lagrange_init(void);
extern void dsolve_solvefor_init(void);
extern void dsolve_undetcoeff_init(void);
extern void dsolve_constcoeff_init(void);
extern void dsolve_euler_init(void);
extern void dsolve_exactode_init(void);
extern void dsolve_specialform_init(void);
extern void dsolve_kovacic_init(void);
extern void dsolve_genairy_init(void);
extern void dsolve_operfactor_init(void);
extern void dsolve_symsquare_init(void);
extern void dsolve_reduce_order_init(void);
extern void dsolve_lower_reduce_init(void);
extern void dsolve_fos_init(void);
extern void dsolve_riccati_init(void);
extern void dsolve_chini_init(void);
extern void dsolve_abel_init(void);
extern void dsolve_abel_air_init(void);
extern void dsolve_polyshift_init(void);
extern void dsolve_lincoeff_init(void);
extern void dsolve_almostlinear_init(void);
extern void dsolve_sepreduced_init(void);
extern void dsolve_linearizable_init(void);
extern void dsolve_lie_init(void);
extern void dsolve_autonomous_init(void);
extern Expr** dsolve_scaling2_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_scaling2_implicit_try(DSolveProblem* P, size_t* nbranch);
extern void dsolve_scaling2_init(void);
extern Expr** dsolve_liouville_implicit_try(DSolveProblem* P, size_t* nbranch);
extern void dsolve_liouville_init(void);
extern void dsolve_lie2_init(void);
extern void dsolve_ifactor_init(void);
extern void dsolve_changevar_init(void);
extern void dsolve_ratsol2_init(void);
extern void dsolve_nonhomog_vop_init(void);
extern void dsolve_piecewise_init(void);
extern void dsolve_frobenius_init(void);
extern void dsolve_normalform_init(void);
extern Expr** dsolve_pde1_solve(DSolveProblem* P);
extern void dsolve_pde1_init(void);
extern Expr** dsolve_pde2_solve(DSolveProblem* P);
extern void dsolve_pde2_init(void);
extern Expr** dsolve_pdequasi_solve(DSolveProblem* P);   /* PDEQuasilinear (Lagrange) */
extern void dsolve_pdequasi_init(void);
extern Expr** dsolve_pdeclairaut_solve(DSolveProblem* P);/* PDEClairaut (complete integral) */
extern void dsolve_pdeclairaut_init(void);
extern Expr** dsolve_pdecharpit_solve(DSolveProblem* P); /* PDECharpit (nonlinear complete integral) */
extern void dsolve_pdecharpit_init(void);
extern void dsolve_pdesep_init(void);  /* DSolve`SeparationOfVariables (pinned-only) */
extern void dsolve_pdeclassify_init(void);  /* PDEClassify (standalone classifier) */
extern Expr* dsolve_wave_ivp_run(DSolveProblem* P);  /* wave d'Alembert IVP (own runner) */
extern void dsolve_wave_init(void);
extern Expr* dsolve_heat_run(DSolveProblem* P);      /* heat-kernel Cauchy IVP (own runner) */
extern void dsolve_heat_init(void);
extern Expr** dsolve_decouple_solve(DSolveProblem* P);
extern Expr** dsolve_triangular_solve(DSolveProblem* P);
/* dsolve_linsys_solve / dsolve_linsys_varcoeff_solve declared in dsolve_linsys.h */
extern Expr** dsolve_sys_reduce_solve(DSolveProblem* P);
extern Expr** dsolve_linsys_commutative_solve(DSolveProblem* P);
extern Expr** dsolve_autosys_solve(DSolveProblem* P);
extern void dsolve_decouple_init(void);
extern void dsolve_triangular_init(void);
extern void dsolve_sys_reduce_init(void);
extern void dsolve_linsys_init(void);
extern void dsolve_linsys_varcoeff_init(void);
extern void dsolve_linsys_commutative_init(void);
extern void dsolve_autosys_init(void);
extern void dsolve_eigenvalue_init(void);

/* ------------------------------------------------------------------ *
 *  Per-command fail-memo (mirror of integrate.c:720-743)              *
 * ------------------------------------------------------------------ */
int g_dsolve_depth = 0;

#define DSOLVE_FAIL_SLOTS 32
static uint64_t ds_fail_epoch = 0;
static int      ds_fail_count = 0;
static struct { uint64_t heq, hv; int method; } ds_fail_tab[DSOLVE_FAIL_SLOTS];

static void ds_fail_sync_epoch(uint64_t tid) {
    if (tid != ds_fail_epoch) { ds_fail_epoch = tid; ds_fail_count = 0; }
}
static bool ds_fail_seen(uint64_t heq, uint64_t hv, int method) {
    for (int i = 0; i < ds_fail_count; i++)
        if (ds_fail_tab[i].heq == heq && ds_fail_tab[i].hv == hv
            && ds_fail_tab[i].method == method)
            return true;
    return false;
}
static void ds_fail_record(uint64_t heq, uint64_t hv, int method) {
    if (ds_fail_seen(heq, hv, method)) return;
    if (ds_fail_count >= DSOLVE_FAIL_SLOTS) return;
    ds_fail_tab[ds_fail_count].heq = heq;
    ds_fail_tab[ds_fail_count].hv = hv;
    ds_fail_tab[ds_fail_count].method = method;
    ds_fail_count++;
}

/* ------------------------------------------------------------------ *
 *  Automatic-answer contract                                          *
 * ------------------------------------------------------------------ */
/* True when e contains Derivative[...][f] for an unknown f of P. */
static bool ds_has_unknown_derivative(const Expr* e, const DSolveProblem* P) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h && h->type == EXPR_FUNCTION && h->data.function.arg_count == 1
        && h->data.function.head->type == EXPR_FUNCTION) {
        const Expr* dh = h->data.function.head->data.function.head;
        const Expr* f  = h->data.function.args[0];
        if (dh->type == EXPR_SYMBOL && dh->data.symbol.name == SYM_Derivative
            && f->type == EXPR_SYMBOL)
            for (size_t i = 0; i < P->nfun; i++)
                if (f->data.symbol.name == P->fun_names[i]) return true;
    }
    if (ds_has_unknown_derivative(h, P)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ds_has_unknown_derivative(e->data.function.args[i], P)) return true;
    return false;
}

/* An AUTOMATIC DSolve returns only solutions.  A truncated SeriesData or a
 * relation still carrying a derivative of the unknown (a reduction of order, e.g.
 * a first integral F(y, y') == C[1]) is not one; those shapes are produced only by
 * an explicitly pinned Method (PowerSeries / FrobeniusSeries / FirstOrderPowerSeries
 * / ReducibleFirstIntegral).  This backstop can only turn such an answer into
 * "unevaluated" -- it never substitutes a different answer, so it cannot mask a
 * wrong one.  The cascade itself no longer runs those methods; this keeps any
 * method that recurses into one, or a future slot, from breaking the contract. */
static bool ds_automatic_answer_admissible(const Expr* result, const DSolveProblem* P) {
    if (ds_has_head(result, SYM_SeriesData)) return false;
    return !ds_has_unknown_derivative(result, P);
}

/* ------------------------------------------------------------------ *
 *  Dispatcher                                                         *
 * ------------------------------------------------------------------ */
Expr* builtin_dsolve(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 3) return NULL;

    uint64_t heq = expr_hash(res->data.function.args[0]);
    uint64_t hv  = expr_hash(res->data.function.args[1]);
    ds_fail_sync_epoch(eval_toplevel_id());

    DSolveProblem P;
    if (!dsolve_parse(res, &P)) { dsolve_problem_free(&P); return NULL; }

    DSolveMethod method = P.method ? ds_method_from_string(P.method) : DS_AUTOMATIC;
    if (method == DS_INVALID)              { dsolve_problem_free(&P); return NULL; }
    if (ds_fail_seen(heq, hv, (int)method)) { dsolve_problem_free(&P); return NULL; }

    /* A NESTED DSolve is a method's speculative sub-solve (a reduction recursing on
     * the reduced equation, a stage-2 quadrature, ...): its result may be discarded,
     * so its diagnostics -- Solve::nongen from an inversion, Integrate::nonelem from
     * a quadrature that does not close -- are probe noise, not user-facing events.
     * Run it muted and roll back its message firings, so they are neither printed
     * nor caught by an enclosing Check[]. */
    bool nested_probe = g_dsolve_depth > 0;
    unsigned long fired0 = mth_msg_fired_count();
    if (nested_probe) mth_msg_suppress_push();
    g_dsolve_depth++;
    /* Every method attempt is speculative: probes and reductions legitimately
     * form 1/0 while classifying an equation that is not of their form.  Mute the
     * cosmetic Power::infy / Infinity::indet family for the whole cascade (the
     * verifier is the real gate) so a decline never leaks arithmetic warnings —
     * mirroring the per-region mute Frobenius/Kovacic already use. */
    arith_warnings_mute_push();
    Expr* result = NULL;
    if (P.is_pde) {
        result = dsolve_run_pde(&P, dsolve_pde1_solve);
        /* first-order nonlinear PDEs (after the linear const-coeff method, which
         * keeps its cleaner explicit form for the shapes it owns): quasilinear
         * Lagrange (characteristics) BEFORE Clairaut, so an equation that is both
         * semilinear and a degenerate (linear-f) Clairaut gets the general
         * arbitrary-function solution rather than the 2-parameter complete
         * integral.  A genuine (nonlinear-f) Clairaut is nonlinear in the
         * derivatives, which quasilinear declines, so Clairaut still claims it. */
        if (!result) result = dsolve_run_pde_implicit(&P, dsolve_pdequasi_solve);
        if (!result) result = dsolve_run_pde_implicit(&P, dsolve_pdeclairaut_solve);
        /* Charpit (fully nonlinear first order, standard forms): after the
         * quasilinear/Clairaut methods, before the second-order factoring. */
        if (!result) result = dsolve_run_pde_implicit(&P, dsolve_pdecharpit_solve);
        if (!result) result = dsolve_run_pde(&P, dsolve_pde2_solve);
        /* wave d'Alembert IVP: PDE + two initial conditions (neq==3); its own
         * runner (multi-equation verify/assemble, not the single-equation
         * dsolve_run_pde).  Declines anything that is not a 1-D wave IVP. */
        if (!result) result = dsolve_wave_ivp_run(&P);
        /* heat-kernel Cauchy IVP: PDE + one initial condition (neq==2). */
        if (!result) result = dsolve_heat_run(&P);
    } else if (P.nfun > 1) {
        if (!result) result = dsolve_run_system(&P, dsolve_decouple_solve);
        if (!result) result = dsolve_run_system(&P, dsolve_triangular_solve);
        /* Higher-order coupled systems: state-augment to first order, then reuse
         * the constant-coefficient engine.  Declines pure first-order systems
         * (that is LinearFirstOrderSystem below). */
        if (!result) result = dsolve_run_system(&P, dsolve_sys_reduce_solve);
        if (!result) result = dsolve_run_system(&P, dsolve_linsys_solve);
        if (!result) result = dsolve_run_system(&P, dsolve_linsys_solve_trig);
        if (!result) result = dsolve_run_system(&P, dsolve_linsys_varcoeff_solve);
        /* 2x2 commutative variable-coeff class A == a(t)I + b(t)K0 (rotation/
         * decay), via Phi == Exp[∫a] Exp[(∫b)K0]; after the scalar-factor solver. */
        if (!result) result = dsolve_run_system(&P, dsolve_linsys_commutative_solve);
        /* 2-D autonomous (incl. nonlinear) via phase-plane orbit + reconstruction */
        if (!result) result = dsolve_run_system(&P, dsolve_autosys_solve);
    } else
    switch (method) {
        case DS_AUTOMATIC:
            /* Factorable / NthAlgebraic run FIRST (matching SymPy): a product- or
             * power-in-the-derivative form is split before the specialists try to
             * match the whole equation.  Both decline cheaply on the common form
             * (irreducible / linear-in-top-derivative), so the linear specialists
             * below still claim their equations. */
            if (!result) result = dsolve_run(&P, dsolve_factorable_try);
            if (!result) result = dsolve_run(&P, dsolve_nth_algebraic_try);
            /* ...its binomial-separable relations (y'^n == A(x) B(y) with a non-
             * elementary y-integral), stashed by the pass above: the n verified first
             * integrals, before the specialists re-attempt the same separation. */
            if (!result) result = dsolve_run_implicit(&P, dsolve_nth_algebraic_binomial_try);
            if (!result) result = dsolve_run(&P, dsolve_quadrature_try);
            if (!result) result = dsolve_run(&P, dsolve_linear1_try);
            if (!result) result = dsolve_run(&P, dsolve_bernoulli_try);
            if (!result) result = dsolve_run(&P, dsolve_homogeneous_try);
            if (!result) result = dsolve_run(&P, dsolve_separable_try);
            /* Separable first integral (implicit): reached when the separated
             * relation Integrate[1/h,y] == Integrate[g,x] + C[1] is elementary but
             * does not invert for y (Log / Root / transcendental) -- e.g.
             * y'==Cot[t]y/(1+y) -> y+Log[y]==Log[Sin t]+C, or the Root-form cubic
             * separables whose IVP is fitted here as C = G(x0,y0) with no
             * inversion.  Placed right after the explicit Separable, mirroring the
             * Exact / ExactImplicit twin. */
            if (!result) result = dsolve_run_implicit(&P, dsolve_separable_implicit_try);
            /* PolynomialShiftSubstitution: y' == R(x) + g(x)(phi(x)+c y)^p with the
             * substitution u = phi+c y -> separable; the deterministic replacement
             * for abaco2_similar's aborting radical [F(x),G(x)]-symmetry cases.
             * Its detection needs a fractional-power-of-(linear-in-y) atom that the
             * standard methods above never produce, so placing it early is safe and
             * — crucially — reaches these radical equations before Exact's/Lagrange's
             * heavier searches spin on the (parameter-laden) radical.  Implicit,
             * branch-safe output. */
            if (!result) result = dsolve_run_implicit(&P, dsolve_polyshift_try);
            /* Linearizable by u = phi(y) (Log/Exp/Sin/Cos/Tan) -> linear/Bernoulli
             * in u.  Deterministic (fixed substitution table) and gated to
             * transcendental-in-y right-hand sides, so it runs ahead of Exact:
             * it produces a clean explicit y = phi^{-1}(H) for the trig/exp/log-
             * linearizable family, and gets there before Exact's integrating-
             * factor search can spin on a transcendental-in-y equation. */
            if (!result) result = dsolve_run(&P, dsolve_linearizable_try);
            if (!result) result = dsolve_run(&P, dsolve_exact_try);
            /* Implicit first integral for an exact ODE whose potential F is
             * transcendental (F == C[1] does not invert for y): F(x, y[x]) ==
             * C[1], verified by the implicit-function rule.  Placed right after
             * explicit Exact so it claims these (e.g. E^(x y)/Log/Sin y exact
             * equations) before a later method spins on the same residual. */
            if (!result) result = dsolve_run_implicit(&P, dsolve_exact_implicit_try);
            if (!result) result = dsolve_run(&P, dsolve_clairaut_try);
            /* Lagrange/d'Alembert (parametric general solution) — Clairaut, its
             * phi(p)==p special case, runs first. */
            if (!result) result = dsolve_run_parametric(&P, dsolve_lagrange_try);
            /* Undetermined coefficients: a tidy particular for UC forcing of a
             * constant-coefficient linear ODE; runs before constcoeff (which is the
             * general variation-of-parameters fallback for any other forcing). */
            /* Piecewise / step / Heaviside forcing IVP: interval-stitching.  Runs
             * BEFORE UndeterminedCoefficients / LinearConstantCoefficients, whose
             * variation-of-parameters leaves Integrate[UnitStep[...]*Cos[...],t]
             * inert (a non-solution).  Gated inside the try-fn to fire ONLY when the
             * forcing carries Piecewise/UnitStep AND it is an IVP, so smooth-forcing
             * linear ODEs still fall through to undetcoeff; a resolved per-interval
             * sub-problem has smooth forcing, so the recursive DSolve declines this
             * method (plus the re-entry guard) — no infinite recursion. */
            if (!result) result = dsolve_run(&P, dsolve_piecewise_try);
            if (!result) result = dsolve_run(&P, dsolve_undetcoeff_try);
            if (!result) result = dsolve_run(&P, dsolve_constcoeff_try);
            if (!result) result = dsolve_run(&P, dsolve_euler_try);
            /* Exact higher-order linear: total derivative d/dx(M[y]) — reduce
             * order once and recurse.  After Euler (cleaner constants there),
             * before the heavier special-function / Kovacic / series machinery. */
            if (!result) result = dsolve_run(&P, dsolve_exactode_try);
            if (!result) result = dsolve_run(&P, dsolve_specialform_try);
            /* Polynomial fundamental set: a homogeneous 2nd-order linear ODE with
             * rational coefficients whose both solutions are polynomials (e.g.
             * (x^2+1)y''-2xy'+2y==0 -> {x, x^2-1}).  Runs BEFORE Kovacic, which would
             * otherwise return the same set in a complex-radical basis that blocks a
             * downstream VariationOfParameters (2.2.26-2592, 2.2.25-2410). */
            if (!result) result = dsolve_run(&P, dsolve_ratsol2_try);
            if (!result) result = dsolve_run(&P, dsolve_kovacic_try);
            /* ChangeOfVariable: a 2nd-order linear ODE with TRANSCENDENTAL
             * coefficients that the direct rational methods above declined may
             * become rational under t = phi(x) (Cos/Sin/Tan) — transform, recurse,
             * back-substitute (e.g. y''+Cot[x]y'+k(k+1)y==0 -> Legendre). */
            if (!result) result = dsolve_run(&P, dsolve_changevar_try);
            /* SpecialFunctionTransform: the polynomial-coefficient class that reaches a
             * Kummer / Whittaker / Bessel row only after a change of variables (shift to
             * the singular point + exponential factor, degenerate Whittaker, quadratic
             * argument z = (x - x0)^2).  Late: Kovacic / ChangeOfVariable find the
             * simpler Liouvillian answer first when there is one; numerically gated. */
            if (!result) result = dsolve_run(&P, dsolve_specialform_transform_try);
            /* General nonhomogeneous variation-of-parameters BACKSTOP: L[y]==g whose
             * homogeneous L[y]==0 is solvable by a special-function / change-of-
             * variable method that does not itself carry forcing.  Solves the
             * homogeneous (recursively), normalises the fundamental set, adds the VoP
             * particular.  Runs AFTER the closed-form 2nd-order methods so it never
             * preempts their (own-forcing) results — it fires only for the residue
             * they all declined (e.g. Sin[x] y'' + (2Sin-Cos) y' + (Sin-Cos) y ==
             * E^-x, via ChangeOfVariable's transcendental homogeneous set). */
            if (!result) result = dsolve_run(&P, dsolve_nonhomog_vop_try);
            /* Third-order symmetric squares (solution space = products of a
             * second-order basis: Airy^2, Bessel products, ...).  Before
             * operfactor, which would otherwise churn on the non-factorable
             * symmetric-square operator. */
            /* Pure-power potential at order >= 3: after the depression y = w u the
             * equation is u^(n) == A x^m u, whose fundamental set is the generalized
             * Airy family x^j 0F_{n-1}.  A recogniser, so it runs before the
             * symmetric-square / factoring searches it would otherwise burn. */
            if (!result) result = dsolve_run(&P, dsolve_genairy_try);
            if (!result) result = dsolve_run(&P, dsolve_symsquare_try);
            /* Higher-order (>=3) reducible linear operators: factor out a first-order
             * right factor (D-r) and recurse.  After Kovacic (owns order 2), before
             * the reduction/series methods. */
            if (!result) result = dsolve_run(&P, dsolve_operfactor_try);
            if (!result) result = dsolve_run(&P, dsolve_reduce_order_try);
            /* General order>=3 reduction when y and all derivatives below some
             * order m>=1 are absent (p=y^(m)); nonlinear-capable, so it runs
             * after the linear operator/reduce_order methods. */
            if (!result) result = dsolve_run(&P, dsolve_lower_reduce_try);
            if (!result) result = dsolve_run(&P, dsolve_fos_try);
            /* Implicit first integral when the autonomous integrand is non-
             * elementary — chiefly y' == f[a x + b y + c] for an ARBITRARY f
             * (problem 159): the inert-integral relation, verified by the
             * implicit-function rule.  Only fires when explicit fos declined. */
            if (!result) result = dsolve_run_implicit(&P, dsolve_fos_implicit_try);
            /* Riccati after fos: fos owns y'==(a x+b y+c)^2 with the cleaner
             * closed form; genuine Riccati (y'==y^2+x, ...) linearises here. */
            if (!result) result = dsolve_run(&P, dsolve_riccati_try);
            /* Chini / Abel (implicit first integral, reducible-to-autonomous
             * sub-class); n=2 Chini is Riccati, already claimed above. */
            if (!result) result = dsolve_run_implicit(&P, dsolve_chini_try);
            if (!result) result = dsolve_run_implicit(&P, dsolve_abel_try);
            /* Abel-2nd-kind / rational-in-y scaling reduction (u=y/s -> separable):
             * claims the rational-in-y forms the polynomial-RHS specialists above
             * reject (D carries y).  Down payment on the full AIR method (M13). */
            if (!result) result = dsolve_run_implicit(&P, dsolve_abel_air_try);
            /* First-order substitution reductions (after the named specialists):
             * LinearCoefficients (shift/homogeneous or separable), AlmostLinear
             * (u=Integrate[g,y] -> linear), SeparableReduced (w=x^n y -> separable). */
            if (!result) result = dsolve_run(&P, dsolve_lincoeff_try);
            if (!result) result = dsolve_run(&P, dsolve_almostlinear_try);
            if (!result) result = dsolve_run_implicit(&P, dsolve_sepreduced_try);
            if (!result) result = dsolve_run(&P, dsolve_autonomous_try);
            /* ...and its implicit companion EARLY, but only where every reduced
             * quadrature is provably non-elementary (p == c Sqrt[P(y)], P squarefree of
             * degree >= 3: the y'' == polynomial(y) family).  There the inert first
             * integral is the answer, and the integrating-factor / point-symmetry
             * searches below only burned the budget first (6 s). */
            if (!result) result = dsolve_run_implicit(&P, dsolve_autonomous_elliptic_try);
            /* Liouville: y'' + g(y)(y')^2 + h(x)y' == 0 (has both y and x, so
             * missing-y/missing-x reductions above decline).  Two quadratures. */
            if (!result) result = dsolve_run(&P, dsolve_liouville_try);
            /* ScalingReduction (2nd order): a deterministic similarity reduction --
             * homogeneous in y -> u = y'/y (first order); isobaric x -> l x,
             * y -> l^k y -> y = x^k u(Log[x]); translation along (1, -rho) ->
             * v = y + rho x -- the last two AUTONOMOUS in u (then AutonomousReduction).
             * Recursive, top level only.  Ahead of the integrating-factor and
             * point-symmetry SEARCHES (heuristic, and together 7 s on
             * (a x + b y) y'^3 + y'' == 0, which this reduces in 1 s); explicit pass,
             * then its implicit pass (relations mapped back and verified numerically
             * against the second-order equation).  MATHILDA_NO_SCALING2 disables both
             * for A/B. */
            if (!result) result = dsolve_run(&P, dsolve_scaling2_try);
            if (!result) result = dsolve_run_implicit(&P, dsolve_scaling2_implicit_try);
            /* ReducibleIntegratingFactor: nonlinear-2nd-order integrating factor of a
             * restricted form (mu(x,y)/mu(x,y')/mu(y,y')) -> first integral R==C[1] ->
             * first-order cascade.  An ALGEBRAIC search (Cheb-Terrab & Roche 1999),
             * cheaper than the symmetry ansatz below and reaching ODEs with no point
             * symmetry, so it runs first; SecondOrderSymmetry is the heavier backstop. */
            if (!result) result = dsolve_run(&P, dsolve_ifactor_try);
            /* SecondOrderSymmetry: the general nonlinear-2nd-order backstop.  Finds
             * a Lie POINT symmetry (polynomial-ansatz determining system) and
             * reduces the order via canonical coordinates -> first-order ODE ->
             * cascade.  Runs after the 2nd-order specialists and before the series
             * fallback, mirroring how first-order Lie (below) sits before series. */
            if (!result) result = dsolve_run(&P, dsolve_lie2_try);
            /* ReducibleFirstIntegral is PINNED-ONLY (Method -> "ReducibleFirstIntegral"):
             * its answer R(x,y,y') == C[1] still carries y', i.e. it is a reduction of
             * order, not a solution, and Mathematica's DSolve never returns one.  From
             * the automatic cascade it pre-empted the autonomous implicit companion
             * below, which answers the same autonomous family in x and y alone. */
            /* AutonomousReduction implicit companion: an autonomous ODE whose reduced
             * separable quadrature Integrate[1/p,y] is non-elementary is returned as the
             * inert first integral Inactive[Integrate][1/p, y[x]] - x == C[1] (M59) --
             * where the explicit autonomous method above declined (spin-guard). */
            if (!result) result = dsolve_run_implicit(&P, dsolve_autonomous_implicit_try);
            /* Liouville's general integral, implicit: where a Liouville quadrature does
             * not close (or g(y), h(x) are arbitrary) the explicit Liouville above
             * declined; the relation with inert integrals, verified at order two. */
            if (!result) result = dsolve_run_implicit(&P, dsolve_liouville_implicit_try);
            /* implicit first-integral fallback: a homogeneous ODE with no explicit
             * inverse (transcendental log-spiral) is returned as G(x,y[x]) == C[1] */
            if (!result) result = dsolve_run_implicit(&P, dsolve_homogeneous_implicit_try);
            /* deterministic linear-coefficients implicit first integral (log-spiral
             * subset): claimed before the Lie `linear` heuristic that also reaches it */
            if (!result) result = dsolve_run_implicit(&P, dsolve_lincoeff_implicit_try);
            /* SolvableForY: the "dp" differentiation method.  Solve F(x,y,y')==0
             * for y, differentiate, and recurse the cascade on the induced
             * first-order ODE -> parametric solution.  Deterministic but expensive
             * (recurses DSolve), so it runs after every named first-order specialist
             * and the implicit fallbacks, and before the heuristic Lie backstop
             * (matching Maple's late dp ordering).  Generalizes Lagrange (the
             * linear-induced-ODE special case).  The SolvableForX mirror is
             * PINNED-ONLY (not here): its automatic corpus yield is ~0 while a
             * degree-1-in-x form with a cubic-denominator induced ODE produces a
             * transcendental branch whose back-substitution verify is very slow,
             * so it is opt-in like FirstOrderPowerSeries / EigenvalueProblem. */
            /* ParametricQuadrature: a first-order ODE with x or y MISSING that no
             * method above could solve for y' (Sin[y'] + y' == x, y'^3 - a x y' + x^3
             * == 0, E^(y' - y) == y'^2 - 1).  The curve F(u, p) == 0 is parametrised
             * (p = t, t u, t +- u, the present variable u isolated) and the missing
             * coordinate is one quadrature.  Not recursive -- one bounded Solve and
             * Integrate per candidate -- so it costs only the equations it targets.
             * Ahead of SolvableForY and the heuristic backstops, which can spin on
             * the same form (§2.2.57-5635: 13 s); it defers (in the try-fn) to
             * NthAlgebraic's implicit pass when every y'-root has a first integral.
             * MATHILDA_NO_PARAMQUAD disables it for A/B. */
            if (!result) result = dsolve_run_parametric(&P, dsolve_paramquad_try);
            if (!result) result = dsolve_run_parametric(&P, dsolve_solvablefory_try);
            /* ...and its implicit slope elimination: when the induced x(p) ODE solves
             * only implicitly, Phi(p, x) == C, each explicit root p_k(x, y) of the
             * ODE gives the implicit solution Phi(p_k, x) == C (stashed by the pass
             * above, so no second solve; re-verified against the original ODE).
             * MATHILDA_NO_SFY_IMPLICIT disables it for A/B. */
            if (!result) result = dsolve_run_implicit(&P, dsolve_solvablefory_implicit_try);
            /* SolvableForX, automatic only where the induced y(p) ODE is linear or
             * Bernoulli (the x-analogue of Lagrange: y^2 y'^3 - x y' + y == 0) --
             * cheap to solve and to verify; the general method stays pinned. */
            if (!result) result = dsolve_run_parametric(&P, dsolve_solvableforx_linear_try);
            /* Homogeneous class G (isobaric): F(l x, l^k y, l^(k-1) y') == l^w F.
             * y = x^k u, t = Log[x] makes it autonomous in u(t); recursed, mapped back.
             * The weight is solved from samples and checked identically; k == 1 is
             * left to DSolve`Homogeneous.  Explicit pass, then its implicit pass
             * (stashed relations, re-verified against the original ODE).
             * Top-level only; MATHILDA_NO_HOMOG_G disables both for A/B. */
            if (!result) result = dsolve_run(&P, dsolve_homog_g_try);
            if (!result) result = dsolve_run_implicit(&P, dsolve_homog_g_implicit_try);
            /* Lie point-symmetry: the general first-order backstop.  Heuristic
             * (underdetermined determining PDE), so it runs after EVERY
             * deterministic specialist — including the homogeneous implicit
             * fallback, which keeps its cleaner form for its own equations — and
             * immediately before the series fallback. */
            if (!result) result = dsolve_run_implicit(&P, dsolve_lie_try);
            /* NB: DSolve`FirstOrderPowerSeries is pinned-only (not in this automatic
             * cascade): a first-order ODE with no closed form stays unevaluated by
             * default (matching SymPy's opt-in 1st_power_series and Mathematica), so
             * a truncated series is offered only on explicit request. */
            /* Nonlinear exactness (the jet peel) — a general backstop, so it runs
             * after every specialist and every implicit fallback, immediately before
             * the series one.  It claims ANY equation whose left side happens to be a
             * total derivative, which from an earlier slot means preempting methods
             * that answer the same equations directly and better: from ExactODE's own
             * early linear slot it took §2.1.2-1143 (y y''' == y' y'') from 0.08 s to
             * 9.96 s for the IDENTICAL answer, and from just after the first-integral
             * group it preempted M59's autonomous implicit companion on §2.1.2-1168
             * (y^2 y''' == y'^3), which owns that order-3 autonomous family.  Here it
             * can only turn a decline into an answer. */
            if (!result) result = dsolve_run(&P, dsolve_exactode_nl_try);
            /* Inverse orientation: solve the reciprocal ODE dx/dy == 1/F for x(y)
             * (Bernoulli / linear / separable in x(y)) and map the implicit first
             * integral back.  It RECURSES the whole engine and so runs ABSOLUTELY
             * LAST among the first-order deterministic backstops — after the
             * nonlinear-exactness jet peel above — so a case that method already
             * answers (e.g. §2.2.44-4323, a 6 s nonlinear-exact solve) never pays
             * the inverse's recursion latency.  Verify-gated by dsolve_run_implicit,
             * re-entry guarded so the inner solve cannot invert again, and
             * MATHILDA_NO_INVERSE disables it for A/B.  Reaches the x(y)-only classes
             * (§2.2.45-4408 Bernoulli-in-x). */
            /* NthAlgebraic's implicit second pass: a first-order ODE of degree >= 2
             * in y' whose y'-branches each solve only IMPLICITLY (y'^2 == (y-1) y^2 ->
             * x/2 + ArcTan[Sqrt[y-1]] == C[1] per branch).  The explicit pass at the
             * front stashed those relations when it declined, so this slot costs no
             * re-solve; it sits LATE so every specialist (Clairaut, Lagrange,
             * SolvableForY, ...) claims the original equation first, and each relation
             * is re-verified against the original ODE by dsolve_run_implicit.  Ahead of
             * InverseOrientation, which would otherwise recurse the engine for nothing.
             * MATHILDA_NO_NA_IMPLICIT disables it for A/B. */
            if (!result) result = dsolve_run_implicit(&P, dsolve_nth_algebraic_implicit_try);
            /* Bernoulli on Expand[F]: a y^n term hidden inside a product (§2.2.54-5332).
             * Late, because expanding makes every RHS quadratic in y look Bernoulli and
             * the specialists above already close (and close faster) most of those. */
            if (!result) result = dsolve_run(&P, dsolve_bernoulli_expanded_try);
            /* Exact by a mu(u) integrating factor, u in {x+y, x-y, x y, x^2 +- y^2}: late
             * for the same reason -- the relation it finds is often implicit where a
             * specialist above closes explicitly (p' == (x + p)^2). */
            if (!result) result = dsolve_run(&P, dsolve_exact_mu_u_try);
            if (!result) result = dsolve_run_implicit(&P, dsolve_exact_mu_u_implicit_try);
            if (!result) result = dsolve_run_implicit(&P, dsolve_inverse_try);
            /* NB: the Frobenius / power-series methods are PINNED-ONLY (Method ->
             * "FrobeniusSeries" / "PowerSeries"), like FirstOrderPowerSeries above: a
             * truncated series is not a solution, so an ODE with no closed form stays
             * unevaluated by default, as in Mathematica (whose series answers live in
             * AsymptoticDSolveValue, never DSolve). */
            /* M61, absolutely last: a 2nd-order linear ODE with rational coefficients,
             * a singular point and polynomial forcing, whose Bessel/hypergeometric
             * homogeneous set makes the variation-of-parameters integrals non-elementary.
             * The particular is returned with those quadratures held inert
             * (Inactive[Integrate]), which is Mathematica's own shape for the family.
             * Last by design: from here the slot can only turn UNEVAL into an answer, so
             * it cannot cost a case that some earlier method already solves, and its
             * cost lands only on equations that were returning nothing. */
            /* Reduction of order from an elementary KNOWN solution (x^m or E^(lambda x),
             * certified by L[y1] == 0), the second solution / particular held as
             * Inactive[Integrate] when non-elementary and the answer self-verified.
             * Late: it only turns an UNEVAL into an answer (§2.2.60-5933 arbitrary f,
             * §2.2.60-5968 symbolic-rate exponential, §2.2.59-5836 forced).
             * MATHILDA_NO_KNOWNSOL disables it for A/B. */
            if (!result) result = dsolve_run(&P, dsolve_knownsol_try);
            if (!result) result = dsolve_run(&P, dsolve_nonhomog_vop_inert_try);
            break;
        case DS_FACTORABLE:   result = dsolve_run(&P, dsolve_factorable_try);  break;
        case DS_NTHALGEBRAIC: result = dsolve_run(&P, dsolve_nth_algebraic_try); break;
        case DS_QUADRATURE:   result = dsolve_run(&P, dsolve_quadrature_try);  break;
        case DS_LINEAR1:      result = dsolve_run(&P, dsolve_linear1_try);     break;
        case DS_BERNOULLI:    result = dsolve_run(&P, dsolve_bernoulli_try);
                              if (!result) result = dsolve_run(&P, dsolve_bernoulli_expanded_try);
                              break;
        case DS_HOMOGENEOUS:
            result = dsolve_run(&P, dsolve_homogeneous_try);
            if (!result) result = dsolve_run_implicit(&P, dsolve_homogeneous_implicit_try);
            break;
        case DS_SEPARABLE:    result = dsolve_run(&P, dsolve_separable_try);   break;
        case DS_EXACT:        result = dsolve_run(&P, dsolve_exact_try);
                              if (!result) result = dsolve_run(&P, dsolve_exact_mu_u_try);
                              break;
        case DS_CLAIRAUT:     result = dsolve_run(&P, dsolve_clairaut_try);    break;
        case DS_LAGRANGE:     result = dsolve_run_parametric(&P, dsolve_lagrange_try); break;
        case DS_SOLVABLEFORY: result = dsolve_run_parametric(&P, dsolve_solvablefory_try); break;
        case DS_SOLVABLEFORX: result = dsolve_run_parametric(&P, dsolve_solvableforx_try); break;
        case DS_PARAMQUAD:    result = dsolve_run_parametric(&P, dsolve_paramquad_try); break;
        case DS_UNDETCOEFF:   result = dsolve_run(&P, dsolve_undetcoeff_try);  break;
        case DS_CONSTCOEFF:   result = dsolve_run(&P, dsolve_constcoeff_try);  break;
        case DS_EULER:        result = dsolve_run(&P, dsolve_euler_try);       break;
        case DS_EXACTODE:     result = dsolve_run(&P, dsolve_exactode_try);
                              if (!result) result = dsolve_run(&P, dsolve_exactode_nl_try);
                              break;
        case DS_SPECIALFORM:  result = dsolve_run(&P, dsolve_specialform_try); break;
        case DS_SFTRANSFORM:  result = dsolve_run(&P, dsolve_specialform_transform_try); break;
        case DS_KOVACIC:      result = dsolve_run(&P, dsolve_kovacic_try);     break;
        case DS_GENAIRY:      result = dsolve_run(&P, dsolve_genairy_try);     break;
        case DS_OPERFACTOR:   result = dsolve_run(&P, dsolve_operfactor_try);  break;
        case DS_SYMSQUARE:    result = dsolve_run(&P, dsolve_symsquare_try);   break;
        case DS_REDUCEORDER:  result = dsolve_run(&P, dsolve_reduce_order_try); break;
        case DS_LOWERREDUCE:  result = dsolve_run(&P, dsolve_lower_reduce_try); break;
        case DS_FOS:          result = dsolve_run(&P, dsolve_fos_try);          break;
        case DS_RICCATI:      result = dsolve_run(&P, dsolve_riccati_try);      break;
        case DS_CHINI:        result = dsolve_run_implicit(&P, dsolve_chini_try); break;
        case DS_ABEL:         result = dsolve_run_implicit(&P, dsolve_abel_try);  break;
        case DS_ABEL_AIR:     result = dsolve_run_implicit(&P, dsolve_abel_air_try); break;
        case DS_POLYSHIFT:    result = dsolve_run_implicit(&P, dsolve_polyshift_try); break;
        case DS_LINCOEFF:
            result = dsolve_run(&P, dsolve_lincoeff_try);
            if (!result) result = dsolve_run_implicit(&P, dsolve_lincoeff_implicit_try);
            break;
        case DS_ALMOSTLINEAR: result = dsolve_run(&P, dsolve_almostlinear_try);  break;
        case DS_SEPREDUCED:   result = dsolve_run_implicit(&P, dsolve_sepreduced_try); break;
        case DS_LINEARIZABLE: result = dsolve_run(&P, dsolve_linearizable_try);  break;
        case DS_LIE:          result = dsolve_run_implicit(&P, dsolve_lie_try);   break;
        case DS_SCALING2:     result = dsolve_run(&P, dsolve_scaling2_try);
                              if (!result) result = dsolve_run_implicit(&P, dsolve_scaling2_implicit_try);
                              break;
        case DS_AUTONOMOUS:   result = dsolve_run(&P, dsolve_autonomous_try);
                              if (!result) result = dsolve_run_implicit(&P, dsolve_autonomous_implicit_try);
                              break;
        case DS_LIOUVILLE:    result = dsolve_run(&P, dsolve_liouville_try);
                              if (!result) result = dsolve_run_implicit(&P, dsolve_liouville_implicit_try);
                              break;
        case DS_IFACTOR:      result = dsolve_run(&P, dsolve_ifactor_try);      break;
        case DS_FIRSTINTEGRAL: result = dsolve_run_first_integral(&P, dsolve_ifactor_first_integral_try); break;
        case DS_FOPOWERSERIES: result = dsolve_run(&P, dsolve_first_order_series_try); break;
        case DS_FROBENIUS:    result = dsolve_run(&P, dsolve_frobenius_try);
                              /* x=0 irregular/obstructed: expand about a nearby
                               * ordinary point instead */
                              if (!result) result = dsolve_run(&P, dsolve_frobenius_shifted_try);
                              break;
        default: break;
    }
    arith_warnings_mute_pop();
    g_dsolve_depth--;
    if (nested_probe) {
        mth_msg_suppress_pop();
        mth_msg_fired_restore(fired0);
    }

    if (result && method == DS_AUTOMATIC && !ds_automatic_answer_admissible(result, &P)) {
        expr_free(result);
        result = NULL;
    }

    result = ds_collapse_identities(result);

    if (!result && g_dsolve_depth == 0) ds_fail_record(heq, hv, (int)method);
    dsolve_problem_free(&P);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Init                                                               *
 * ------------------------------------------------------------------ */
static Expr* mk_opt(const char* name, Expr* val) {
    return expr_new_function(expr_new_symbol(SYM_Rule),
                             (Expr*[]){ expr_new_symbol(name), val }, 2);
}

void dsolve_init(void) {
    intern_symbol("IncludeSingularSolutions");
    symtab_add_builtin("DSolve", builtin_dsolve);
    symtab_get_def("DSolve")->attributes |= ATTR_PROTECTED;

    Expr* opts = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){
        mk_opt("GeneratedParameters",      expr_new_symbol("C")),
        mk_opt("Assumptions",              expr_new_symbol("True")),
        mk_opt("Method",                   expr_new_symbol("Automatic")),
        mk_opt("IncludeSingularSolutions", expr_new_symbol("False"))
    }, 4);
    symtab_set_options("DSolve", opts);

    /* method registrars */
    dsolve_factorable_init();
    dsolve_nth_algebraic_init();
    dsolve_quadrature_init();
    dsolve_linear1_init();
    dsolve_bernoulli_init();
    dsolve_homogeneous_init();
    dsolve_separable_init();
    dsolve_exact_init();
    dsolve_clairaut_init();
    dsolve_lagrange_init();
    dsolve_solvefor_init();
    dsolve_inverse_init();
    dsolve_knownsol_init();
    dsolve_homog_g_init();
    dsolve_paramquad_init();
    tc_register_unwind_counter(&g_dsolve_depth);   /* restored on a timeout unwind */
    dsolve_undetcoeff_init();
    dsolve_constcoeff_init();
    dsolve_euler_init();
    dsolve_exactode_init();
    dsolve_specialform_init();
    dsolve_kovacic_init();
    dsolve_genairy_init();
    dsolve_operfactor_init();
    dsolve_symsquare_init();
    dsolve_reduce_order_init();
    dsolve_lower_reduce_init();
    dsolve_fos_init();
    dsolve_riccati_init();
    dsolve_chini_init();
    dsolve_abel_init();
    dsolve_abel_air_init();
    dsolve_polyshift_init();
    dsolve_lincoeff_init();
    dsolve_almostlinear_init();
    dsolve_sepreduced_init();
    dsolve_linearizable_init();
    dsolve_lie_init();
    dsolve_autonomous_init();
    dsolve_scaling2_init();
    dsolve_liouville_init();
    dsolve_lie2_init();
    dsolve_ifactor_init();
    dsolve_changevar_init();
    dsolve_ratsol2_init();
    dsolve_nonhomog_vop_init();
    dsolve_piecewise_init();
    dsolve_frobenius_init();
    dsolve_normalform_init();
    dsolve_pde1_init();
    dsolve_pde2_init();
    dsolve_pdequasi_init();
    dsolve_pdeclairaut_init();
    dsolve_pdecharpit_init();
    dsolve_pdesep_init();
    dsolve_pdeclassify_init();
    dsolve_wave_init();
    dsolve_heat_init();
    dsolve_decouple_init();
    dsolve_triangular_init();
    dsolve_sys_reduce_init();
    dsolve_linsys_init();
    dsolve_linsys_varcoeff_init();
    dsolve_linsys_commutative_init();
    dsolve_autosys_init();
    dsolve_eigenvalue_init();
}
