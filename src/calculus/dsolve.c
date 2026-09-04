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
#include "dsolve_linsys.h"

#include "../sym_names.h"
#include "../sym_intern.h"
#include "../eval.h"
#include "../symtab.h"
#include "../attr.h"
#include "../expr.h"
#include "../arithmetic.h"

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
    DS_KOVACIC,
    DS_OPERFACTOR,
    DS_REDUCEORDER,
    DS_FOS,
    DS_RICCATI,
    DS_CHINI,
    DS_ABEL,
    DS_LINCOEFF,
    DS_ALMOSTLINEAR,
    DS_SEPREDUCED,
    DS_LIE,
    DS_AUTONOMOUS,
    DS_LIOUVILLE,
    DS_FOPOWERSERIES,
    DS_FROBENIUS,
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
    if (strcmp(s, "UndeterminedCoefficients") == 0) return DS_UNDETCOEFF;
    if (strcmp(s, "LinearConstantCoefficients") == 0) return DS_CONSTCOEFF;
    if (strcmp(s, "EulerCauchy")      == 0) return DS_EULER;
    if (strcmp(s, "ExactODE")         == 0) return DS_EXACTODE;
    if (strcmp(s, "SpecialFunctionForm") == 0) return DS_SPECIALFORM;
    if (strcmp(s, "Kovacic")            == 0) return DS_KOVACIC;
    if (strcmp(s, "OperatorFactor")     == 0) return DS_OPERFACTOR;
    if (strcmp(s, "ReductionOfOrder") == 0) return DS_REDUCEORDER;
    if (strcmp(s, "FirstOrderSubstitution") == 0) return DS_FOS;
    if (strcmp(s, "Riccati")              == 0) return DS_RICCATI;
    if (strcmp(s, "Chini")                == 0) return DS_CHINI;
    if (strcmp(s, "Abel")                 == 0) return DS_ABEL;
    if (strcmp(s, "LinearCoefficients")   == 0) return DS_LINCOEFF;
    if (strcmp(s, "AlmostLinear")         == 0) return DS_ALMOSTLINEAR;
    if (strcmp(s, "SeparableReduced")     == 0) return DS_SEPREDUCED;
    if (strcmp(s, "LieSymmetry")          == 0) return DS_LIE;
    if (strcmp(s, "LieGroup")             == 0) return DS_LIE;
    if (strcmp(s, "AutonomousReduction") == 0) return DS_AUTONOMOUS;
    if (strcmp(s, "Liouville")            == 0) return DS_LIOUVILLE;
    if (strcmp(s, "FirstOrderPowerSeries") == 0) return DS_FOPOWERSERIES;
    if (strcmp(s, "FrobeniusSeries")     == 0) return DS_FROBENIUS;
    if (strcmp(s, "PowerSeries")         == 0) return DS_FROBENIUS;
    return DS_INVALID;
}

/* method try-functions + registrars (one file each) */
extern Expr** dsolve_factorable_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_nth_algebraic_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_quadrature_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_linear1_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_bernoulli_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_homogeneous_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_homogeneous_implicit_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_separable_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_exact_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_clairaut_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_lagrange_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_undetcoeff_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_constcoeff_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_euler_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_exactode_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_specialform_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_kovacic_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_operfactor_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_reduce_order_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_fos_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_riccati_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_chini_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_abel_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_lincoeff_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_lincoeff_implicit_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_almostlinear_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_sepreduced_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_lie_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_autonomous_try(DSolveProblem* P, size_t* nbranch);
extern Expr** dsolve_liouville_try(DSolveProblem* P, size_t* nbranch);
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
extern void dsolve_undetcoeff_init(void);
extern void dsolve_constcoeff_init(void);
extern void dsolve_euler_init(void);
extern void dsolve_exactode_init(void);
extern void dsolve_specialform_init(void);
extern void dsolve_kovacic_init(void);
extern void dsolve_operfactor_init(void);
extern void dsolve_reduce_order_init(void);
extern void dsolve_fos_init(void);
extern void dsolve_riccati_init(void);
extern void dsolve_chini_init(void);
extern void dsolve_abel_init(void);
extern void dsolve_lincoeff_init(void);
extern void dsolve_almostlinear_init(void);
extern void dsolve_sepreduced_init(void);
extern void dsolve_lie_init(void);
extern void dsolve_autonomous_init(void);
extern void dsolve_liouville_init(void);
extern void dsolve_frobenius_init(void);
extern void dsolve_normalform_init(void);
extern Expr** dsolve_pde1_solve(DSolveProblem* P);
extern void dsolve_pde1_init(void);
extern Expr** dsolve_pde2_solve(DSolveProblem* P);
extern void dsolve_pde2_init(void);
extern void dsolve_pdesep_init(void);  /* DSolve`SeparationOfVariables (pinned-only) */
extern void dsolve_pdeclassify_init(void);  /* PDEClassify (standalone classifier) */
extern Expr** dsolve_decouple_solve(DSolveProblem* P);
extern Expr** dsolve_triangular_solve(DSolveProblem* P);
/* dsolve_linsys_solve / dsolve_linsys_varcoeff_solve declared in dsolve_linsys.h */
extern void dsolve_decouple_init(void);
extern void dsolve_triangular_init(void);
extern void dsolve_linsys_init(void);
extern void dsolve_linsys_varcoeff_init(void);
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
        if (!result) result = dsolve_run_pde(&P, dsolve_pde2_solve);
    } else if (P.nfun > 1) {
        if (!result) result = dsolve_run_system(&P, dsolve_decouple_solve);
        if (!result) result = dsolve_run_system(&P, dsolve_triangular_solve);
        if (!result) result = dsolve_run_system(&P, dsolve_linsys_solve);
        if (!result) result = dsolve_run_system(&P, dsolve_linsys_varcoeff_solve);
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
            if (!result) result = dsolve_run(&P, dsolve_quadrature_try);
            if (!result) result = dsolve_run(&P, dsolve_linear1_try);
            if (!result) result = dsolve_run(&P, dsolve_bernoulli_try);
            if (!result) result = dsolve_run(&P, dsolve_homogeneous_try);
            if (!result) result = dsolve_run(&P, dsolve_separable_try);
            if (!result) result = dsolve_run(&P, dsolve_exact_try);
            if (!result) result = dsolve_run(&P, dsolve_clairaut_try);
            /* Lagrange/d'Alembert (parametric general solution) — Clairaut, its
             * phi(p)==p special case, runs first. */
            if (!result) result = dsolve_run_parametric(&P, dsolve_lagrange_try);
            /* Undetermined coefficients: a tidy particular for UC forcing of a
             * constant-coefficient linear ODE; runs before constcoeff (which is the
             * general variation-of-parameters fallback for any other forcing). */
            if (!result) result = dsolve_run(&P, dsolve_undetcoeff_try);
            if (!result) result = dsolve_run(&P, dsolve_constcoeff_try);
            if (!result) result = dsolve_run(&P, dsolve_euler_try);
            /* Exact higher-order linear: total derivative d/dx(M[y]) — reduce
             * order once and recurse.  After Euler (cleaner constants there),
             * before the heavier special-function / Kovacic / series machinery. */
            if (!result) result = dsolve_run(&P, dsolve_exactode_try);
            if (!result) result = dsolve_run(&P, dsolve_specialform_try);
            if (!result) result = dsolve_run(&P, dsolve_kovacic_try);
            /* Higher-order (>=3) reducible linear operators: factor out a first-order
             * right factor (D-r) and recurse.  After Kovacic (owns order 2), before
             * the reduction/series methods. */
            if (!result) result = dsolve_run(&P, dsolve_operfactor_try);
            if (!result) result = dsolve_run(&P, dsolve_reduce_order_try);
            if (!result) result = dsolve_run(&P, dsolve_fos_try);
            /* Riccati after fos: fos owns y'==(a x+b y+c)^2 with the cleaner
             * closed form; genuine Riccati (y'==y^2+x, ...) linearises here. */
            if (!result) result = dsolve_run(&P, dsolve_riccati_try);
            /* Chini / Abel (implicit first integral, reducible-to-autonomous
             * sub-class); n=2 Chini is Riccati, already claimed above. */
            if (!result) result = dsolve_run_implicit(&P, dsolve_chini_try);
            if (!result) result = dsolve_run_implicit(&P, dsolve_abel_try);
            /* First-order substitution reductions (after the named specialists):
             * LinearCoefficients (shift/homogeneous or separable), AlmostLinear
             * (u=Integrate[g,y] -> linear), SeparableReduced (w=x^n y -> separable). */
            if (!result) result = dsolve_run(&P, dsolve_lincoeff_try);
            if (!result) result = dsolve_run(&P, dsolve_almostlinear_try);
            if (!result) result = dsolve_run_implicit(&P, dsolve_sepreduced_try);
            if (!result) result = dsolve_run(&P, dsolve_autonomous_try);
            /* Liouville: y'' + g(y)(y')^2 + h(x)y' == 0 (has both y and x, so
             * missing-y/missing-x reductions above decline).  Two quadratures. */
            if (!result) result = dsolve_run(&P, dsolve_liouville_try);
            /* implicit first-integral fallback: a homogeneous ODE with no explicit
             * inverse (transcendental log-spiral) is returned as G(x,y[x]) == C[1] */
            if (!result) result = dsolve_run_implicit(&P, dsolve_homogeneous_implicit_try);
            /* deterministic linear-coefficients implicit first integral (log-spiral
             * subset): claimed before the Lie `linear` heuristic that also reaches it */
            if (!result) result = dsolve_run_implicit(&P, dsolve_lincoeff_implicit_try);
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
            /* series fallback: always-available, so it runs last */
            if (!result) result = dsolve_run(&P, dsolve_frobenius_try);
            /* very last resort: if x=0 was an irregular/obstructed singular point,
             * expand the series about a nearby ordinary point so a linear ODE with
             * any ordinary point never returns unevaluated (auto-dispatch only). */
            if (!result) result = dsolve_run(&P, dsolve_frobenius_shifted_try);
            break;
        case DS_FACTORABLE:   result = dsolve_run(&P, dsolve_factorable_try);  break;
        case DS_NTHALGEBRAIC: result = dsolve_run(&P, dsolve_nth_algebraic_try); break;
        case DS_QUADRATURE:   result = dsolve_run(&P, dsolve_quadrature_try);  break;
        case DS_LINEAR1:      result = dsolve_run(&P, dsolve_linear1_try);     break;
        case DS_BERNOULLI:    result = dsolve_run(&P, dsolve_bernoulli_try);   break;
        case DS_HOMOGENEOUS:
            result = dsolve_run(&P, dsolve_homogeneous_try);
            if (!result) result = dsolve_run_implicit(&P, dsolve_homogeneous_implicit_try);
            break;
        case DS_SEPARABLE:    result = dsolve_run(&P, dsolve_separable_try);   break;
        case DS_EXACT:        result = dsolve_run(&P, dsolve_exact_try);       break;
        case DS_CLAIRAUT:     result = dsolve_run(&P, dsolve_clairaut_try);    break;
        case DS_LAGRANGE:     result = dsolve_run_parametric(&P, dsolve_lagrange_try); break;
        case DS_UNDETCOEFF:   result = dsolve_run(&P, dsolve_undetcoeff_try);  break;
        case DS_CONSTCOEFF:   result = dsolve_run(&P, dsolve_constcoeff_try);  break;
        case DS_EULER:        result = dsolve_run(&P, dsolve_euler_try);       break;
        case DS_EXACTODE:     result = dsolve_run(&P, dsolve_exactode_try);    break;
        case DS_SPECIALFORM:  result = dsolve_run(&P, dsolve_specialform_try); break;
        case DS_KOVACIC:      result = dsolve_run(&P, dsolve_kovacic_try);     break;
        case DS_OPERFACTOR:   result = dsolve_run(&P, dsolve_operfactor_try);  break;
        case DS_REDUCEORDER:  result = dsolve_run(&P, dsolve_reduce_order_try); break;
        case DS_FOS:          result = dsolve_run(&P, dsolve_fos_try);          break;
        case DS_RICCATI:      result = dsolve_run(&P, dsolve_riccati_try);      break;
        case DS_CHINI:        result = dsolve_run_implicit(&P, dsolve_chini_try); break;
        case DS_ABEL:         result = dsolve_run_implicit(&P, dsolve_abel_try);  break;
        case DS_LINCOEFF:
            result = dsolve_run(&P, dsolve_lincoeff_try);
            if (!result) result = dsolve_run_implicit(&P, dsolve_lincoeff_implicit_try);
            break;
        case DS_ALMOSTLINEAR: result = dsolve_run(&P, dsolve_almostlinear_try);  break;
        case DS_SEPREDUCED:   result = dsolve_run_implicit(&P, dsolve_sepreduced_try); break;
        case DS_LIE:          result = dsolve_run_implicit(&P, dsolve_lie_try);   break;
        case DS_AUTONOMOUS:   result = dsolve_run(&P, dsolve_autonomous_try);   break;
        case DS_LIOUVILLE:    result = dsolve_run(&P, dsolve_liouville_try);    break;
        case DS_FOPOWERSERIES: result = dsolve_run(&P, dsolve_first_order_series_try); break;
        case DS_FROBENIUS:    result = dsolve_run(&P, dsolve_frobenius_try);    break;
        default: break;
    }
    arith_warnings_mute_pop();
    g_dsolve_depth--;

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
    dsolve_undetcoeff_init();
    dsolve_constcoeff_init();
    dsolve_euler_init();
    dsolve_exactode_init();
    dsolve_specialform_init();
    dsolve_kovacic_init();
    dsolve_operfactor_init();
    dsolve_reduce_order_init();
    dsolve_fos_init();
    dsolve_riccati_init();
    dsolve_chini_init();
    dsolve_abel_init();
    dsolve_lincoeff_init();
    dsolve_almostlinear_init();
    dsolve_sepreduced_init();
    dsolve_lie_init();
    dsolve_autonomous_init();
    dsolve_liouville_init();
    dsolve_frobenius_init();
    dsolve_normalform_init();
    dsolve_pde1_init();
    dsolve_pde2_init();
    dsolve_pdesep_init();
    dsolve_pdeclassify_init();
    dsolve_decouple_init();
    dsolve_triangular_init();
    dsolve_linsys_init();
    dsolve_linsys_varcoeff_init();
    dsolve_eigenvalue_init();
}
