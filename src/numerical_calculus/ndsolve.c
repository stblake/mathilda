/* Mathilda — NDSolve front-end: compile ODE initial-value problems into an
 * NdProblem (first-order reduced system), dispatch to a stepper, and assemble
 * the InterpolatingFunction result.  Numerical kernels live in the ndsolve_*
 * modules; see ndsolve_common.h. */
#include "ndsolve.h"
#include "ndsolve_common.h"
#include "fem.h"
#include "ndsolve_compile.h"
#include "ndsolve_events.h"
#include "../sym_names.h"
#include "../sym_intern.h"
#include "../eval.h"
#include "../symtab.h"
#include "../loadmodule.h"
#include "../attr.h"
#include "../common.h"
#include "../numeric.h"
#include "../nc_accuracy.h"   /* shared AccuracyGoal/PrecisionGoal handling */
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* Interned option-name symbols not already in sym_names (set in ndsolve_init). */
static const char* OPT_MaxSteps;
static const char* OPT_MaxStepSize;
static const char* OPT_MaxStepFraction;
static const char* OPT_StartingStepSize;
static const char* OPT_DependentVariables;

/* ------------------------------------------------------------------ *
 *  Small expr helpers (nd_call1/nd_call2/nd_eval_to_double/            *
 *  nd_replace_all are shared, defined in ndsolve_common.c)             *
 * ------------------------------------------------------------------ */

/* Build the literal  u[x]  (order 0) or  Derivative[m][u][x]  (m>=1). */
static Expr* nd_make_funcapp(const char* fname, int order, const char* xvar) {
    if (order == 0) return nd_call1(fname, expr_new_symbol(xvar));
    Expr* mord[1] = { expr_new_integer(order) };
    Expr* d1 = expr_new_function(expr_new_symbol(SYM_Derivative), mord, 1);
    Expr* uu[1] = { expr_new_symbol(fname) };
    Expr* d2 = expr_new_function(d1, uu, 1);
    Expr* xa[1] = { expr_new_symbol(xvar) };
    return expr_new_function(d2, xa, 1);
}

/* Match  u[arg]  (order 0) or  Derivative[m][u][arg].  Fills interned fname,
 * order and the (borrowed) argument. */
static bool nd_match_funcapp(const Expr* e, const char** fname, int* order, const Expr** arg) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 1) return false;
    const Expr* head = e->data.function.head;
    const Expr* a = e->data.function.args[0];
    if (head->type == EXPR_SYMBOL) { *fname = head->data.symbol.name; *order = 0; *arg = a; return true; }
    if (head->type == EXPR_FUNCTION && head->data.function.arg_count == 1
        && head->data.function.args[0]->type == EXPR_SYMBOL) {
        const Expr* d1 = head->data.function.head;   /* Derivative[m] */
        if (d1->type == EXPR_FUNCTION && d1->data.function.arg_count == 1
            && d1->data.function.head->type == EXPR_SYMBOL
            && d1->data.function.head->data.symbol.name == SYM_Derivative
            && d1->data.function.args[0]->type == EXPR_INTEGER) {
            *fname = head->data.function.args[0]->data.symbol.name;
            *order = (int)d1->data.function.args[0]->data.integer;
            *arg = a;
            return true;
        }
    }
    return false;
}

/* Recursively scan `e`, updating maxorder[k] for each dependent function and
 * flagging whether any funcapp is at the independent variable (=> an ODE). */
static void nd_scan(const Expr* e, const char** funcs, size_t nfun, const char* xvar,
                    int* maxorder, bool* has_x) {
    if (!e) return;
    const char* fn; int ord; const Expr* arg;
    if (nd_match_funcapp(e, &fn, &ord, &arg)) {
        for (size_t k = 0; k < nfun; k++)
            if (fn == funcs[k]) {
                if (ord > maxorder[k]) maxorder[k] = ord;
                if (arg->type == EXPR_SYMBOL && arg->data.symbol.name == xvar) *has_x = true;
            }
    }
    if (e->type == EXPR_FUNCTION) {
        nd_scan(e->data.function.head, funcs, nfun, xvar, maxorder, has_x);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            nd_scan(e->data.function.args[i], funcs, nfun, xvar, maxorder, has_x);
    }
}

/* ------------------------------------------------------------------ *
 *  Option parsing                                                     *
 * ------------------------------------------------------------------ */
static bool nd_is_option(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return false;
    if (h->data.symbol.name != SYM_Rule && h->data.symbol.name != SYM_RuleDelayed) return false;
    return e->data.function.arg_count == 2 && e->data.function.args[0]->type == EXPR_SYMBOL;
}

/* Parse a goal option value.  Delegates to the shared parser (nc_accuracy.h):
 * Automatic -> NC_GOAL_AUTO (-1), Infinity -> NC_GOAL_OFF (+Inf ==
 * HUGE_VAL), MachinePrecision -> ~15.95 digits (the newly-accepted setting),
 * positive number -> itself.  On any unrecognised value *out is left untouched,
 * keeping the caller's default.  The sole consumer, nd_resolve_tol, already
 * treats `>= HUGE_VAL` as Infinity (criterion disabled) and `< 0` as Automatic
 * (-> WorkingPrecision/2), so the +Inf sentinel disables exactly as before. */
static void nd_parse_goal(Expr* v, double* out) {
    nc_parse_goal(v, out);
}

static void nd_apply_option(const Expr* opt, NdOpts* o) {
    const Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    if (lhs->data.symbol.name == SYM_Method) {
        /* "Name" or {"Name", subopts...} or Symbol */
        Expr* mv = eval_and_free(expr_copy(rhs));
        if (mv && mv->type == EXPR_STRING) o->method = intern_symbol(mv->data.string);
        else if (mv && mv->type == EXPR_SYMBOL) o->method = mv->data.symbol.name;
        else if (head_is(mv, SYM_List) && mv->data.function.arg_count >= 1
                 && mv->data.function.args[0]->type == EXPR_STRING)
            o->method = intern_symbol(mv->data.function.args[0]->data.string);
        expr_free(mv);
    } else if (lhs->data.symbol.name == SYM_AccuracyGoal) {
        Expr* v = eval_and_free(expr_copy(rhs)); nd_parse_goal(v, &o->acc_goal); expr_free(v);
    } else if (lhs->data.symbol.name == SYM_PrecisionGoal) {
        Expr* v = eval_and_free(expr_copy(rhs)); nd_parse_goal(v, &o->prec_goal); expr_free(v);
    } else if (lhs->data.symbol.name == SYM_WorkingPrecision) {
        Expr* v = eval_and_free(expr_copy(rhs));
        if (v && v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_MachinePrecision) {
            o->spec = numeric_machine_spec(); o->wp_bits = 53;
        } else {
            double dig;
            if (nd_to_double(v, &dig) && dig > NUMERIC_MACHINE_PRECISION_DIGITS) {
#ifdef USE_MPFR
                o->wp_bits = numeric_digits_to_bits(dig);
                o->spec.mode = NUMERIC_MODE_MPFR; o->spec.bits = o->wp_bits;
#else
                o->wp_bits = 53;
#endif
            } else { o->spec = numeric_machine_spec(); o->wp_bits = 53; }
        }
        expr_free(v);
    } else if (lhs->data.symbol.name == OPT_MaxSteps) {
        Expr* v = eval_and_free(expr_copy(rhs)); double d;
        if (v && v->type == EXPR_INTEGER) o->max_steps = v->data.integer;
        else if (nd_to_double(v, &d) && d > 0) o->max_steps = (int64_t)d;
        expr_free(v);
    } else if (lhs->data.symbol.name == OPT_MaxStepSize) {
        Expr* v = eval_and_free(expr_copy(rhs)); double d;
        if (nd_to_double(v, &d) && d > 0) o->max_step_size = d;
        expr_free(v);
    } else if (lhs->data.symbol.name == OPT_MaxStepFraction) {
        Expr* v = eval_and_free(expr_copy(rhs)); double d;
        if (nd_to_double(v, &d) && d > 0) o->max_step_fraction = d;
        expr_free(v);
    } else if (lhs->data.symbol.name == OPT_StartingStepSize) {
        Expr* v = eval_and_free(expr_copy(rhs)); double d;
        if (nd_to_double(v, &d) && d > 0) o->starting_step = d;
        expr_free(v);
    } else if (lhs->data.symbol.name == SYM_InterpolationOrder) {
        Expr* v = eval_and_free(expr_copy(rhs));
        if (v && v->type == EXPR_INTEGER) o->interp_order = (int)v->data.integer;
        expr_free(v);
    } else if (lhs->data.symbol.name == SYM_StepMonitor) {
        o->step_monitor = rhs;
    } else if (lhs->data.symbol.name == SYM_EvaluationMonitor) {
        o->eval_monitor = rhs;
    } else if (lhs->data.symbol.name == SYM_NormFunction) {
        o->norm_function = rhs;
    }
    (void)OPT_DependentVariables;
}

/* ------------------------------------------------------------------ *
 *  Problem cleanup                                                    *
 * ------------------------------------------------------------------ */
static void nd_problem_free(NdProblem* P) {
    nd_operator_free(P->op);
    nd_compiled_free(P->compiled);
    nd_ac_prec_free(P);
    if (P->f) { for (size_t i = 0; i < P->d; i++) expr_free(P->f[i]); free(P->f); }
    if (P->ysym) { for (size_t i = 0; i < P->d; i++) expr_free(P->ysym[i]); free(P->ysym); }
    if (P->jac) {
        for (size_t i = 0; i < P->d; i++) {
            for (size_t j = 0; j < P->d; j++) expr_free(P->jac[i][j]);
            free(P->jac[i]);
        }
        free(P->jac);
    }
    free(P->Y0);
    expr_free(P->t0_expr); expr_free(P->tmin_expr); expr_free(P->tmax_expr);
    if (P->Y0_expr) { for (size_t i = 0; i < P->y0_expr_n; i++) expr_free(P->Y0_expr[i]); free(P->Y0_expr); }
    free(P->bind_y);
    free(P->fun_names);
    free(P->fun_state0);
    nd_events_free(P->events);
}

/* ------------------------------------------------------------------ *
 *  Core solver                                                        *
 * ------------------------------------------------------------------ */
static void nd_warn(const char* tag, const char* msg) {
    mth_message("NDSolve", tag, "%s", msg);
}

/* ------------------------------------------------------------------ *
 *  Boundary-value problems                                            *
 * ------------------------------------------------------------------ */

/* Every dependent-function application in e at a numeric point: appends the
 * point to pts (deduplicated, |p - q| <= 1e-12 (1 + |p|) is one point). */
static void nd_collect_points(const Expr* e, const char** fnames, size_t nfun,
                              NumericSpec spec, double* pts, size_t* npts, size_t cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const char* fn; int ord; const Expr* arg;
    if (nd_match_funcapp(e, &fn, &ord, &arg)) {
        for (size_t k = 0; k < nfun; k++) if (fn == fnames[k]) {
            double p;
            if (nd_eval_to_double((Expr*)arg, spec, &p)) {
                bool seen = false;
                for (size_t i = 0; i < *npts; i++)
                    if (fabs(pts[i] - p) <= 1e-12 * (1.0 + fabs(p))) seen = true;
                if (!seen && *npts < cap) pts[(*npts)++] = p;
            }
            return;
        }
    }
    nd_collect_points(e->data.function.head, fnames, nfun, spec, pts, npts, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        nd_collect_points(e->data.function.args[i], fnames, nfun, spec, pts, npts, cap);
}

/* Is the condition the plain initial value  f^(k)[p] == v  (k below f's
 * order, v numeric, possibly complex)? */
static bool nd_plain_ic(const Expr* c, const char** fnames, size_t nfun,
                        const int* maxorder, NumericSpec spec) {
    if (!head_is((Expr*)c, SYM_Equal) || c->data.function.arg_count != 2) return false;
    for (int side = 0; side < 2; side++) {
        const Expr* fa = c->data.function.args[side];
        Expr* val = c->data.function.args[1 - side];
        const char* fn; int ord; const Expr* arg;
        if (!nd_match_funcapp(fa, &fn, &ord, &arg)) continue;
        for (size_t k = 0; k < nfun; k++) if (fn == fnames[k]) {
            double p, re, im;
            if (ord >= maxorder[k] || !nd_eval_to_double((Expr*)arg, spec, &p)) return false;
            return nd_eval_to_double(val, spec, &re) || nd_split_reim(val, spec, &re, &im);
        }
    }
    return false;
}

/* The conditions of an ODE system pose an initial-value problem when every
 * one is a plain initial value, all at ONE point, and no state entry is given
 * twice.  Anything else -- conditions at two or more points (u[0] == 0,
 * u'[1] == 0), a condition that couples entries (u[0] + u'[0] == 1,
 * u[0] == u[1]) -- is a boundary-value problem.  Before this test such a
 * system was read as an IVP whose start point was the LAST condition's point,
 * with every condition's value taken as an initial value there: a silent
 * wrong answer (-u'' == E^(2x), u[0] == 0, u'[1] == 0 gave u(1) = 0, not
 * 2.0973). */
static bool nd_is_bvp(Expr** eqitems, size_t neq, const bool* is_ode,
                      const char** fnames, size_t nfun, const int* maxorder,
                      NumericSpec spec) {
    double pts[8]; size_t npts = 0;
    size_t nplain = 0, ncond = 0;
    for (size_t e = 0; e < neq; e++) {
        if (is_ode[e]) continue;
        ncond++;
        if (nd_plain_ic(eqitems[e], fnames, nfun, maxorder, spec)) nplain++;
        nd_collect_points(eqitems[e], fnames, nfun, spec, pts, &npts, 8);
    }
    if (ncond == 0) return false;
    if (npts > 1 || nplain < ncond) return true;
    /* one point, all plain: a repeated entry (u[0] == 1, u[0] == 2) */
    for (size_t e = 0; e < neq; e++) {
        if (is_ode[e]) continue;
        for (size_t f = e + 1; f < neq; f++) {
            if (is_ode[f]) continue;
            const Expr* ce = eqitems[e]; const Expr* cf = eqitems[f];
            for (int se = 0; se < 2; se++) for (int sf = 0; sf < 2; sf++) {
                const char* n1; const char* n2; int o1, o2; const Expr* a1; const Expr* a2;
                if (nd_match_funcapp(ce->data.function.args[se], &n1, &o1, &a1)
                    && nd_match_funcapp(cf->data.function.args[sf], &n2, &o2, &a2)
                    && n1 == n2 && o1 == o2) {
                    bool known = false;
                    for (size_t k = 0; k < nfun; k++) if (n1 == fnames[k]) known = true;
                    if (known) return true;
                }
            }
        }
    }
    return false;
}

static unsigned bvp_load_gen = 0;      /* symtab_generation() of the attempt */
static bool bvp_load_succeeded = false;

/* Solve the boundary-value problem by shooting (src/internal/NDSolveBVP.m).
 * Returns the NDSolve result, or NULL (unevaluated) after a message. */
static Expr* nd_solve_bvp(Expr** A, size_t argc, Expr** eqitems, size_t neq,
                          const bool* is_ode, const char** fnames, size_t nfun,
                          const int* maxorder) {
    if (bvp_load_gen != symtab_generation()) {
        bvp_load_gen = symtab_generation();
        mth_msg_suppress_push();
        bvp_load_succeeded = mathilda_load_module("NDSolveBVP.m") != 0;
        mth_msg_suppress_pop();
    }
    if (!bvp_load_succeeded) {
        nd_warn("bvnf", "cannot locate src/internal/NDSolveBVP.m for the boundary-value solver");
        return NULL;
    }
    size_t nodes = 0;
    for (size_t e = 0; e < neq; e++) if (is_ode[e]) nodes++;
    Expr** odes = malloc(sizeof(Expr*) * (nodes ? nodes : 1));
    Expr** conds = malloc(sizeof(Expr*) * (neq - nodes ? neq - nodes : 1));
    size_t io = 0, ic = 0;
    for (size_t e = 0; e < neq; e++) {
        if (is_ode[e]) odes[io++] = expr_copy(eqitems[e]);
        else conds[ic++] = expr_copy(eqitems[e]);
    }
    Expr** heads = malloc(sizeof(Expr*) * nfun);
    Expr** ords = malloc(sizeof(Expr*) * nfun);
    for (size_t k = 0; k < nfun; k++) {
        heads[k] = expr_new_symbol(fnames[k]);
        ords[k] = expr_new_integer(maxorder[k]);
    }
    size_t nargs = 6 + (argc - 3);
    Expr** args = malloc(sizeof(Expr*) * nargs);
    args[0] = expr_new_function(expr_new_symbol(SYM_List), odes, io);
    args[1] = expr_new_function(expr_new_symbol(SYM_List), conds, ic);
    args[2] = expr_copy(A[1]);
    args[3] = expr_copy(A[2]);
    args[4] = expr_new_function(expr_new_symbol(SYM_List), heads, nfun);
    args[5] = expr_new_function(expr_new_symbol(SYM_List), ords, nfun);
    for (size_t i = 3; i < argc; i++) args[6 + i - 3] = expr_copy(A[i]);
    free(odes); free(conds); free(heads); free(ords);
    Expr* call = expr_new_function(expr_new_symbol(intern_symbol("NDSolveBVP`Shoot")),
                                   args, nargs);
    free(args);
    Expr* r = eval_and_free(call);
    Expr* out = NULL;
    const char* tag = NULL;
    if (r && head_is(r, SYM_List) && r->data.function.arg_count >= 1
        && r->data.function.args[0]->type == EXPR_STRING)
        tag = r->data.function.args[0]->data.string;
    if (tag && (strcmp(tag, "ok") == 0 || strcmp(tag, "berr") == 0)
        && r->data.function.arg_count >= 2 && head_is(r->data.function.args[1], SYM_List)) {
        if (strcmp(tag, "berr") == 0)
            nd_warn("berr", "the boundary-value residual did not reach the precision goal; "
                            "the returned solution may be inaccurate");
        out = r->data.function.args[1];
        r->data.function.args[1] = NULL;
    } else if (tag && strcmp(tag, "bvdim") == 0) {
        nd_warn("bvdim", "the number of boundary conditions does not match the order of the system");
    } else {
        nd_warn("bvfail", "the shooting method for the boundary-value problem failed");
    }
    expr_free(r);
    return out;
}

static void nd_collect_unknowns(const Expr* e, const char*** names, size_t* n, size_t* cap);

static Expr* ndsolve_core(Expr* res, const char* forced_method) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 3) return NULL;
    Expr** A = res->data.function.args;
    size_t argc = res->data.function.arg_count;

    /* ---- an elliptic PDE on a region, {x, y} \[Element] reg: finite elements ---- */
    if (fem_is_region_spec(A[2])) return fem_ndsolve(res);

    /* ---- options (trailing) ---- */
    NdOpts o; nd_opts_default(&o);
    /* Default goals.  PrecisionGoal is Automatic (-> WorkingPrecision/2).  At
     * machine precision the default AccuracyGoal is MachinePrecision (a
     * ~15.95-digit absolute-error floor); at a higher WorkingPrecision it is
     * Automatic too (-> WorkingPrecision/2).  A fixed 1e-16 absolute floor
     * there would cap every small-valued solution at machine accuracy
     * whatever the working precision: a decaying solution (Airy Ai integrated
     * forward to x = 8, ~5e-8) kept only ~8 correct digits at
     * WorkingPrecision 50, and an unstable direction amplified the floor's
     * error.  NaN marks "not given", so an explicit AccuracyGoal -> ... (any
     * value, MachinePrecision included) still overrides either default. */
    o.acc_goal = NAN;
    size_t pos_end = argc;
    while (pos_end > 3 && nd_is_option(A[pos_end - 1])) pos_end--;
    /* allow exactly 3 positional args; options after */
    for (size_t i = 3; i < argc; i++) if (nd_is_option(A[i])) nd_apply_option(A[i], &o);
    if (isnan(o.acc_goal))
        o.acc_goal = o.wp_bits > 53 ? -1.0 : NUMERIC_MACHINE_PRECISION_DIGITS;
    /* "MethodOfLines" is a PDE controller, not a time-integration stepper — it
     * must not clobber a Method-> time-integration choice given alongside it. */
    if (forced_method && strcmp(forced_method, "MethodOfLines") != 0)
        o.method = intern_symbol(forced_method);

    /* ---- PDE?  Two or more independent-variable ranges {t,..},{x,..}[,{y,..}]
     * means a partial differential equation: hand off to the Method-of-Lines
     * front-end, which discretizes space into the ODE system solved below. ---- */
    {
        size_t nranges = 0;
        for (size_t i = 2; i < pos_end; i++) {
            const Expr* r = A[i];
            if (head_is((Expr*)r, SYM_List) && r->data.function.arg_count == 3
                && r->data.function.args[0]->type == EXPR_SYMBOL) nranges++;
        }
        if (nranges >= 2) return nd_mol_solve(res, &o, forced_method);
    }

    /* ---- range {x, xmin, xmax} ---- */
    Expr* range = A[2];
    if (!head_is(range, SYM_List) || range->data.function.arg_count != 3) return NULL;
    if (range->data.function.args[0]->type != EXPR_SYMBOL) return NULL;
    const char* xvar = range->data.function.args[0]->data.symbol.name;
    double tmin, tmax;
    if (!nd_eval_to_double(range->data.function.args[1], o.spec, &tmin)) return NULL;
    if (!nd_eval_to_double(range->data.function.args[2], o.spec, &tmax)) return NULL;
    /* {x, a, b} with a > b is the interval [b, a] (integrated backward from a
     * initial condition at its top); it was left unevaluated */
    bool range_rev = tmin > tmax;
    if (range_rev) { double sw = tmin; tmin = tmax; tmax = sw; }

    /* ---- dependent functions ---- */
    Expr* funcs = A[1];
    Expr** funcitems; size_t nfun; bool own_items = false;
    if (head_is(funcs, SYM_List)) {
        nfun = funcs->data.function.arg_count; funcitems = funcs->data.function.args;
    } else { nfun = 1; funcitems = malloc(sizeof(Expr*)); funcitems[0] = funcs; own_items = true; }
    if (nfun == 0) { if (own_items) free(funcitems); return NULL; }
    const char** fnames = malloc(sizeof(char*) * nfun);
    bool applied = false;
    for (size_t k = 0; k < nfun; k++) {
        Expr* it = funcitems[k];
        if (it->type == EXPR_SYMBOL) fnames[k] = it->data.symbol.name;
        else if (it->type == EXPR_FUNCTION && it->data.function.head->type == EXPR_SYMBOL
                 && it->data.function.arg_count == 1) {
            fnames[k] = it->data.function.head->data.symbol.name; applied = true;
        } else { free(fnames); if (own_items) free(funcitems); return NULL; }
    }
    if (own_items) free(funcitems);

    /* ---- normalize equations to a flat list ---- */
    Expr* eqns = A[0];
    Expr** eqall; size_t nall;
    if (head_is(eqns, SYM_List)) { nall = eqns->data.function.arg_count; eqall = eqns->data.function.args; }
    else { nall = 1; eqall = &eqns; }
    /* WhenEvent[...] items are events, not equations: set them aside (borrowed). */
    Expr** eqitems = malloc(sizeof(Expr*) * (nall ? nall : 1));
    Expr** evitems = malloc(sizeof(Expr*) * (nall ? nall : 1));
    size_t neq = 0, nev = 0;
    for (size_t e = 0; e < nall; e++) {
        if (nd_is_when_event(eqall[e])) evitems[nev++] = eqall[e];
        else eqitems[neq++] = eqall[e];
    }

    /* ---- unknowns the equations differentiate but the call does not list
     * (NDSolve[{x' == y, y' == -x, ...}, x, ...]): solved for with the rest
     * and not returned, as the reference does.  They used to be read as
     * unknown symbols in the right-hand side (NDSolve::nrnum) or their
     * conditions as boundary conditions. ---- */
    size_t nreq = nfun;
    {
        const char** extra = NULL; size_t ne = 0, ce = 0;
        for (size_t e = 0; e < neq; e++) nd_collect_unknowns(eqitems[e], &extra, &ne, &ce);
        for (size_t i = 0; i < ne; i++) {
            bool have = false;
            for (size_t k = 0; k < nfun && !have; k++) have = fnames[k] == extra[i];
            if (!have) { fnames = realloc(fnames, sizeof(char*) * (nfun + 1)); fnames[nfun++] = extra[i]; }
        }
        free(extra);
    }

    /* ---- classify + orders ---- */
    int* maxorder = calloc(nfun, sizeof(int));
    bool* is_ode = calloc(neq, sizeof(bool));
    for (size_t e = 0; e < neq; e++) {
        bool hx = false;
        nd_scan(eqitems[e], fnames, nfun, xvar, maxorder, &hx);
        is_ode[e] = hx;
    }
    size_t d = 0;
    size_t* base = malloc(sizeof(size_t) * nfun);
    bool order_ok = true;
    for (size_t k = 0; k < nfun; k++) {
        base[k] = d;
        if (maxorder[k] < 1) order_ok = false;
        d += (size_t)(maxorder[k] < 1 ? 1 : maxorder[k]);
    }
    if (!order_ok || d == 0) {
        free(maxorder); free(is_ode); free(base); free(fnames); free(eqitems); free(evitems);
        nd_warn("nodeqs", "could not identify a well-posed ODE system");
        return NULL;
    }

    /* ---- boundary-value problem?  Hand off to the shooting solver. ---- */
    if (nd_is_bvp(eqitems, neq, is_ode, fnames, nfun, maxorder, o.spec)) {
        Expr* r = NULL;
        if (nev > 0)
            nd_warn("bvev", "WhenEvent is not supported for boundary-value problems");
        else
            r = nd_solve_bvp(A, argc, eqitems, neq, is_ode, fnames, nfun, maxorder);
        free(maxorder); free(is_ode); free(base); free(fnames); free(eqitems); free(evitems);
        return r;
    }

    /* ---- reduced-state symbols + substitution pairs ---- */
    Expr** ysym = malloc(sizeof(Expr*) * d);
    const char** ysym_name = malloc(sizeof(char*) * d);
    Expr** lits = malloc(sizeof(Expr*) * d);   /* substitution literals  */
    Expr** subs = malloc(sizeof(Expr*) * d);   /* -> reduced-state syms   */
    for (size_t k = 0; k < nfun; k++) {
        int nk = maxorder[k];
        for (int m = 0; m < nk; m++) {
            size_t gi = base[k] + (size_t)m;
            char buf[64];
            snprintf(buf, sizeof(buf), "NDSolve`y%zu", gi);
            ysym_name[gi] = intern_symbol(buf);
            ysym[gi] = expr_new_symbol(ysym_name[gi]);
            lits[gi] = nd_make_funcapp(fnames[k], m, xvar);
            subs[gi] = expr_new_symbol(ysym_name[gi]);
        }
    }

    /* ---- Block-localize x and the reduced-state symbols ---- */
    NdProblem P; memset(&P, 0, sizeof(P));
    P.d = d; P.spec = o.spec; P.tvar = xvar;
    P.eval_monitor = o.eval_monitor;
    P.nfun = nreq; P.fun_applied = applied;   /* only the requested ones are returned */
    P.fun_names = fnames;
    P.fun_state0 = malloc(sizeof(size_t) * nfun);
    for (size_t k = 0; k < nfun; k++) P.fun_state0[k] = base[k];
    P.ysym = ysym;
    P.bind_y = malloc(sizeof(NdBind) * d);
    nd_bind_snapshot(&P.bind_t, xvar);
    for (size_t i = 0; i < d; i++) nd_bind_snapshot(&P.bind_y[i], ysym_name[i]);

    /* ---- build reduced RHS f[] ---- */
    P.f = calloc(d, sizeof(Expr*));
    /* chain equations: f[base_k+m] = y_{base_k+m+1} for m < n_k-1 */
    for (size_t k = 0; k < nfun; k++)
        for (int m = 0; m + 1 < maxorder[k]; m++)
            P.f[base[k] + (size_t)m] = expr_copy(ysym[base[k] + (size_t)m + 1]);

    /* solve each ODE for its top derivative -> f[base_k + n_k - 1] */
    bool build_ok = true;
    for (size_t e = 0; e < neq && build_ok; e++) {
        if (!is_ode[e]) continue;
        /* leading function = one whose derivative order is highest in this eqn */
        int* eo = calloc(nfun, sizeof(int)); bool hx = false;
        nd_scan(eqitems[e], fnames, nfun, xvar, eo, &hx);
        int besto = 0; size_t bestk = 0; bool found = false;
        for (size_t k = 0; k < nfun; k++) if (eo[k] > besto) { besto = eo[k]; bestk = k; found = true; }
        free(eo);
        if (!found || besto != maxorder[bestk]) { continue; }   /* not the top eqn */

        if (!head_is(eqitems[e], SYM_Equal) || eqitems[e]->data.function.arg_count != 2) {
            build_ok = false; break;
        }
        /* residual R = lhs - rhs */
        Expr* R = nd_call2(SYM_Subtract, expr_copy(eqitems[e]->data.function.args[0]),
                                         expr_copy(eqitems[e]->data.function.args[1]));
        Expr* topLit = nd_make_funcapp(fnames[bestk], besto, xvar);
        Expr* Psym = expr_new_symbol("NDSolve`P");
        /* Rp = R with topLit -> P */
        Expr* Rp = nd_replace_all(R, &topLit, &Psym, 1);
        /* a = D[Rp, P] ; b = Rp|P=0 ; g = -b/a */
        Expr* dargs[2] = { expr_copy(Rp), expr_copy(Psym) };
        Expr* aE = eval_and_free(expr_new_function(expr_new_symbol(SYM_D), dargs, 2));
        Expr* zero = expr_new_integer(0);
        Expr* bE = nd_replace_all(expr_copy(Rp), &Psym, &zero, 1);
        expr_free(zero);
        Expr* inv[2] = { expr_copy(aE), expr_new_integer(-1) };
        Expr* invE = expr_new_function(expr_new_symbol(SYM_Power), inv, 2);
        Expr* g3[3] = { expr_new_integer(-1), expr_copy(bE), invE };
        Expr* g = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), g3, 3));
        /* substitute state forms */
        Expr* gsub = nd_replace_all(g, lits, subs, d);
        expr_free(P.f[base[bestk] + (size_t)maxorder[bestk] - 1]);
        P.f[base[bestk] + (size_t)maxorder[bestk] - 1] = gsub;
        /* R was consumed by nd_replace_all; do not free it again. */
        expr_free(topLit); expr_free(Psym); expr_free(Rp);
        expr_free(aE); expr_free(bE);
    }
    for (size_t i = 0; i < d; i++) if (!P.f[i]) build_ok = false;

    /* ---- initial conditions -> Y0 (Y0im for complex ICs), t0 ---- */
    P.Y0 = calloc(d, sizeof(double));
    P.Y0_expr = calloc(d, sizeof(Expr*));
    P.y0_expr_n = d;
    P.tmin_expr = expr_copy(range->data.function.args[range_rev ? 2 : 1]);
    P.tmax_expr = expr_copy(range->data.function.args[range_rev ? 1 : 2]);
    double* Y0im = calloc(d, sizeof(double));
    bool* have_ic = calloc(d, sizeof(bool));
    double t0 = tmin; bool t0_set = false;
    for (size_t e = 0; e < neq; e++) {
        if (is_ode[e]) continue;
        if (!head_is(eqitems[e], SYM_Equal) || eqitems[e]->data.function.arg_count != 2) continue;
        Expr* s0 = eqitems[e]->data.function.args[0];
        Expr* s1 = eqitems[e]->data.function.args[1];
        const char* fn; int ord; const Expr* arg;
        Expr* fa = NULL; Expr* val = NULL;
        if (nd_match_funcapp(s0, &fn, &ord, &arg)) { fa = s0; val = s1; }
        else if (nd_match_funcapp(s1, &fn, &ord, &arg)) { fa = s1; val = s0; }
        if (!fa) continue;
        (void)fa;
        size_t k; bool fk = false;
        for (k = 0; k < nfun; k++) if (fn == fnames[k]) { fk = true; break; }
        if (!fk || ord >= maxorder[k]) continue;
        double pt;
        if (nd_eval_to_double((Expr*)arg, o.spec, &pt)) {
            t0 = pt; t0_set = true;
            expr_free(P.t0_expr); P.t0_expr = expr_copy((Expr*)arg);
        }
        double re, im = 0.0;
        if (nd_eval_to_double(val, o.spec, &re)) { /* real IC */ }
        else if (!nd_split_reim(val, o.spec, &re, &im)) continue;   /* complex IC */
        size_t gi = base[k] + (size_t)ord;
        P.Y0[gi] = re; Y0im[gi] = im; have_ic[gi] = true;
        expr_free(P.Y0_expr[gi]); P.Y0_expr[gi] = expr_copy(val);
    }
    (void)t0_set;
    bool ic_ok = true;
    for (size_t i = 0; i < d; i++) if (!have_ic[i]) ic_ok = false;
    free(have_ic);
    P.t0 = t0; P.tmin = tmin; P.tmax = tmax;

    /* ---- complex-valued ODE: realify (split Re/Im, double the dimension) ----
     * Detected when any reduced RHS carries I or any IC is non-real.  The real
     * machinery (sampler, steppers, Hermite output) is then reused; the output
     * recombines the paired components into a complex value. */
    bool cplx = false;
    for (size_t i = 0; i < d; i++) if (Y0im[i] != 0.0) cplx = true;
    for (size_t i = 0; i < d && !cplx; i++) if (nd_has_imaginary(P.f[i])) cplx = true;
    if (cplx && build_ok && ic_ok) {
#ifdef USE_MPFR
        if (numeric_spec_is_mpfr(o.spec)) {
            nd_warn("cmplx", "complex-valued ODEs are integrated at machine precision");
            o.spec = numeric_machine_spec(); P.spec = o.spec;
        }
#endif
        double* Y0re = malloc(d * sizeof(double));
        memcpy(Y0re, P.Y0, d * sizeof(double));
        if (!nd_realify(&P, d, Y0re, Y0im)) build_ok = false;
        free(Y0re);
    }
    free(Y0im);

    /* ---- events: WhenEvent[event, action] in terms of the reduced state ---- */
    bool ev_ok = true;
    if (nev > 0 && build_ok && ic_ok) {
        if (cplx) {
            nd_warn("cev", "WhenEvent is not supported for complex-valued ODEs");
            ev_ok = false;
        } else {
            /* u^(m)[x] -> reduced symbol (m below the order), and the top
             * derivative u^(n)[x] -> its right-hand side */
            size_t nsub = d + nfun;
            Expr** el = malloc(sizeof(Expr*) * nsub);
            Expr** es = malloc(sizeof(Expr*) * nsub);
            for (size_t i = 0; i < d; i++) { el[i] = lits[i]; es[i] = subs[i]; }
            for (size_t k = 0; k < nfun; k++) {
                el[d + k] = nd_make_funcapp(fnames[k], maxorder[k], xvar);
                es[d + k] = P.f[base[k] + (size_t)maxorder[k] - 1];
            }
            ev_ok = nd_events_build(&P, evitems, nev, el, es, nsub, fnames, base, maxorder, nfun);
            for (size_t k = 0; k < nfun; k++) expr_free(el[d + k]);
            free(el); free(es);
        }
    }

    Expr* result = NULL;
    if (build_ok && ic_ok && ev_ok) {
        const NdStepper* S = nd_lookup_stepper(o.method);
        if (!S) { nd_warn("method", "unknown Method; using Automatic"); S = nd_default_stepper(); }
#ifdef USE_MPFR
        if (numeric_spec_is_mpfr(o.spec)) {
            /* Arbitrary-precision path: explicit (DOPRI5/RK4) or, for stiff
             * implicit/multistep methods, the MPFR variable-order BDF. */
            result = nd_solve_mpfr(&P, &o, S);
            nd_bind_restore(&P.bind_t);
            for (size_t i = 0; i < P.d; i++) nd_bind_restore(&P.bind_y[i]);
        } else
#endif
        {
        NdSolution sol; nd_solution_init(&sol, P.d);
        NdStatus st = nd_integrate(&P, S, &o, &sol);
        /* restore bindings before building the result (uses x, u symbolically) */
        nd_bind_restore(&P.bind_t);
        for (size_t i = 0; i < P.d; i++) nd_bind_restore(&P.bind_y[i]);
        if (st == ND_ERR_MAXSTEPS) {
            /* Ran out of the step budget before reaching the endpoint: the
             * returned solution is truncated, so its error over the requested
             * interval is unbounded.  Report the shortfall against the documented
             * combined tolerance at unit solution scale (NDSolve's error control
             * is scale-relative) and return the partial (best) solution. */
            double wp = o.wp_bits > 0 ? numeric_bits_to_digits(o.wp_bits)
                                      : NUMERIC_MACHINE_PRECISION_DIGITS;
            nc_warn_goal("NDSolve", HUGE_VAL,
                         nc_combined_tol(o.acc_goal, o.prec_goal, 1.0, wp));
        }
        else if (st == ND_ERR_STEPSIZE) nd_warn("ndsz", "step size effectively zero; singularity or stiffness suspected");
        else if (st == ND_ERR_NONCONV) nd_warn("ndcf", "corrector failed to converge");
        else if (st == ND_ERR_SAMPLE) nd_warn("nrnum", "right-hand side did not evaluate to a number");
        result = cplx ? nd_build_result_complex(&P, &o, &sol)
                      : nd_build_result(&P, &o, &sol);
        nd_solution_free(&sol);
        }
    } else {
        nd_bind_restore(&P.bind_t);
        for (size_t i = 0; i < P.d; i++) nd_bind_restore(&P.bind_y[i]);
        if (!ic_ok) nd_warn("underdet", "insufficient initial conditions to determine the solution");
    }

    /* cleanup */
    for (size_t i = 0; i < d; i++) { expr_free(lits[i]); expr_free(subs[i]); }
    free(lits); free(subs); free(ysym_name);
    free(maxorder); free(is_ode); free(base); free(eqitems); free(evitems);
    nd_problem_free(&P);
    return result;   /* NULL leaves NDSolve[...] unevaluated */
}

/* ------------------------------------------------------------------ *
 *  Builtins                                                           *
 * ------------------------------------------------------------------ */
Expr* builtin_ndsolve(Expr* res) {
    return ndsolve_core(res, NULL);
}

/* Every function the equations differentiate, Derivative[...][g][...] with g a
 * symbol, in order of first appearance (WhenEvent items are not equations). */
static void nd_collect_unknowns(const Expr* e, const char*** names, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION || nd_is_when_event(e)) return;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_FUNCTION && h->data.function.arg_count == 1
        && h->data.function.args[0]->type == EXPR_SYMBOL
        && head_is((Expr*)h->data.function.head, SYM_Derivative)) {
        const char* g = h->data.function.args[0]->data.symbol.name;
        bool seen = false;
        for (size_t i = 0; i < *n && !seen; i++) seen = (*names)[i] == g;
        if (!seen) {
            if (*n == *cap) { *cap = *cap ? 2 * *cap : 4; *names = realloc(*names, sizeof(char*) * *cap); }
            (*names)[(*n)++] = g;
        }
    }
    nd_collect_unknowns(h, names, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        nd_collect_unknowns(e->data.function.args[i], names, n, cap);
}

/* NDSolveValue[eqns, expr, {x, a, b}, ...]: solve as NDSolve does, for every
 * function the equations differentiate, and return expr with the solution
 * substituted (y -> InterpolatingFunction, so y gives the function, y[1] a
 * number, {y[1], y'[1]} a list).  Left unevaluated when NDSolve is. */
static Expr* builtin_ndsolvevalue(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 3) return NULL;
    Expr** A = res->data.function.args;
    size_t argc = res->data.function.arg_count;
    const char** names = NULL; size_t n = 0, cap = 0;
    nd_collect_unknowns(A[0], &names, &n, &cap);
    if (n == 0) {
        free(names);
        mth_message("NDSolveValue", "dvnf", "the equations contain no derivative of an unknown function");
        return NULL;
    }
    Expr** fs = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) fs[i] = expr_new_symbol(names[i]);
    free(names);
    Expr** args = malloc(sizeof(Expr*) * argc);
    args[0] = expr_copy(A[0]);
    args[1] = expr_new_function(expr_new_symbol(SYM_List), fs, n);
    free(fs);
    for (size_t i = 2; i < argc; i++) args[i] = expr_copy(A[i]);
    Expr* call = expr_new_function(expr_new_symbol(intern_symbol("NDSolve")), args, argc);
    free(args);
    Expr* sol = eval_and_free(call);
    Expr* out = NULL;
    if (head_is(sol, SYM_List) && sol->data.function.arg_count >= 1
        && head_is(sol->data.function.args[0], SYM_List)) {
        Expr* rules = expr_copy(sol->data.function.args[0]);
        out = eval_and_free(nd_call2(SYM_ReplaceAll, expr_copy(A[1]), rules));
    }
    expr_free(sol);
    return out;
}

/* One shared wrapper for every NDSolve`Method builtin: read the method from the
 * call's own head (the text after the backtick) and force it. */
static Expr* builtin_ndsolve_method(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* full = res->data.function.head->data.symbol.name;
    const char* tick = strrchr(full, '`');
    const char* method = tick ? tick + 1 : full;
    return ndsolve_core(res, method);
}

void ndsolve_init(void) {
    OPT_MaxSteps         = intern_symbol("MaxSteps");
    OPT_MaxStepSize      = intern_symbol("MaxStepSize");
    OPT_MaxStepFraction  = intern_symbol("MaxStepFraction");
    OPT_StartingStepSize = intern_symbol("StartingStepSize");
    OPT_DependentVariables = intern_symbol("DependentVariables");

    symtab_add_builtin("NDSolve", builtin_ndsolve);
    symtab_get_def("NDSolve")->attributes |= ATTR_PROTECTED;

    /* Options[NDSolve] — defaults mirror nd_opts_default() (see
     * ndsolve_common.c): every value shown here is the resolved behavior when
     * the option is omitted. `Automatic` corresponds to the sentinel defaults
     * (NULL method, -1 goals/steps/order, 0 step sizes) the solver interprets
     * internally; MaxStepFraction is a concrete 1/10. */
    {
        struct { const char* lhs; Expr* rhs; } defs[] = {
            { SYM_Method,             expr_new_symbol(SYM_Automatic) },
            { SYM_WorkingPrecision,   expr_new_symbol(SYM_MachinePrecision) },
            { SYM_AccuracyGoal,       expr_new_symbol(SYM_Automatic) },
            { SYM_PrecisionGoal,      expr_new_symbol(SYM_Automatic) },
            { OPT_MaxSteps,           expr_new_symbol(SYM_Automatic) },
            { OPT_MaxStepSize,        expr_new_symbol(SYM_Automatic) },
            { OPT_MaxStepFraction,    expr_new_function(expr_new_symbol(SYM_Rational),
                                        (Expr*[]){ expr_new_integer(1), expr_new_integer(10) }, 2) },
            { OPT_StartingStepSize,   expr_new_symbol(SYM_Automatic) },
            { SYM_InterpolationOrder, expr_new_symbol(SYM_Automatic) },
            { SYM_StepMonitor,        expr_new_symbol(SYM_None) },
            { SYM_EvaluationMonitor,  expr_new_symbol(SYM_None) },
            { SYM_NormFunction,       expr_new_symbol(SYM_Automatic) },
            { OPT_DependentVariables, expr_new_symbol(SYM_Automatic) },
        };
        size_t nopt = sizeof(defs) / sizeof(defs[0]);
        Expr** rules = malloc(nopt * sizeof(Expr*));
        for (size_t i = 0; i < nopt; i++)
            rules[i] = expr_new_function(expr_new_symbol(SYM_Rule),
                (Expr*[]){ expr_new_symbol(defs[i].lhs), defs[i].rhs }, 2);
        Expr* opts = expr_new_function(expr_new_symbol(SYM_List), rules, nopt);
        free(rules);
        symtab_set_options("NDSolve", opts);  /* takes ownership */
    }

    symtab_set_docstring("NDSolve",
        "NDSolve[eqns, u, {x, xmin, xmax}]\n"
        "\tsolves the ordinary differential equations eqns numerically for the\n"
        "\tfunction u over xmin <= x <= xmax, returning {{u -> "
        "InterpolatingFunction[...]}}.\n"
        "NDSolve[eqns, {u1, u2, ...}, {x, xmin, xmax}] solves a system.\n"
        "NDSolve[eqns, u[x], {x, xmin, xmax}] gives u[x] -> InterpolatingFunction[...][x].\n"
        "\tHigher-order equations (u''[x] == ...) are reduced to first order.\n"
        "\tConditions at two or more points (u[0] == 0, u'[1] == 0), or coupling\n"
        "\tvalues (u[0] + u'[0] == 1, u[0] == u[1]), pose a boundary-value\n"
        "\tproblem, solved by shooting with Newton's method.\n"
        "\tWhenEvent[event, action] items among eqns are events (see WhenEvent).\n"
        "NDSolve[eqns, u, {t, tmin, tmax}, {x, xmin, xmax}] solves a partial\n"
        "\tdifferential equation over a rectangular region by the method of lines,\n"
        "\tgiving a 2-D InterpolatingFunction applied as u[t, x].\n"
        "\tThe arguments are evaluated (equations may be held in a symbol).\n"
        "\tAt WorkingPrecision above machine precision the Automatic method is\n"
        "\ta Taylor-series integrator when the right-hand side is elementary.\n"
        "\tOptions: Method, WorkingPrecision, AccuracyGoal, PrecisionGoal,\n"
        "\tMaxSteps, MaxStepSize, MaxStepFraction, StartingStepSize,\n"
        "\tInterpolationOrder, StepMonitor, EvaluationMonitor.");

    symtab_add_builtin("NDSolveValue", builtin_ndsolvevalue);
    symtab_get_def("NDSolveValue")->attributes |= ATTR_PROTECTED;
    {
        Expr* nopts = symtab_get_options("NDSolve");
        if (nopts) symtab_set_options("NDSolveValue", expr_copy(nopts));
    }
    symtab_set_docstring("NDSolveValue",
        "NDSolveValue[eqns, expr, {x, xmin, xmax}]\n"
        "\tsolves the differential equations eqns numerically, as NDSolve does,\n"
        "\tfor every function they differentiate, and gives expr with the\n"
        "\tsolution substituted: NDSolveValue[eqns, y, ...] is the\n"
        "\tInterpolatingFunction, NDSolveValue[eqns, y[1], ...] a number,\n"
        "\tNDSolveValue[eqns, {y[1], y'[1]}, ...] a list.  PDE ranges\n"
        "\t{t, ...}, {x, ...} and every NDSolve option are accepted.");

    /* WhenEvent[event, action]: an inert, held event specification read by
     * NDSolve (ndsolve_events.c). */
    symtab_get_def("WhenEvent")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("WhenEvent",
        "WhenEvent[event, action]\n"
        "\tgiven among the equations of NDSolve, runs action whenever event\n"
        "\toccurs during the integration.  event is lhs == rhs (a zero crossing\n"
        "\tin either direction), an inequality (fires when it becomes True), any\n"
        "\tother condition (fires when it becomes True), or a list of events.\n"
        "\tThe event time is located to the integrator's accuracy.  action is\n"
        "\t\"StopIntegration\" (the solution's domain ends at the event),\n"
        "\t\"RemoveEvent\", a state reset u'[t] -> expr or a list of resets\n"
        "\t(integration restarts from the event with the new state), any\n"
        "\texpression (evaluated at the event with t and the state bound, e.g.\n"
        "\tSow[t]), or a list of actions.  WhenEvent holds its arguments.");

    /* Register one builtin per stepper as NDSolve`Name. */
    char buf[64];
    for (size_t i = 0; i < nd_stepper_count(); i++) {
        const NdStepper* S = nd_stepper_at(i);
        snprintf(buf, sizeof(buf), "NDSolve`%s", S->name);
        symtab_add_builtin(buf, builtin_ndsolve_method);
        symtab_get_def(buf)->attributes |= ATTR_PROTECTED;
        symtab_set_docstring(buf, S->doc);
    }
    /* Wolfram-facing method-name aliases (resolved by nd_lookup_stepper). */
    static const char* const aliases[] = {
        "ExplicitRungeKutta", "RungeKutta", "StiffnessSwitching"
    };
    for (size_t i = 0; i < sizeof(aliases)/sizeof(aliases[0]); i++) {
        snprintf(buf, sizeof(buf), "NDSolve`%s", aliases[i]);
        symtab_add_builtin(buf, builtin_ndsolve_method);
        symtab_get_def(buf)->attributes |= ATTR_PROTECTED;
        symtab_set_docstring(buf,
            "NDSolve`<method>[eqns, u, {x, xmin, xmax}] — solve forcing this "
            "time-integration method (see NDSolve).");
    }

    /* PDE controller: the Method-of-Lines spatial discretization. */
    symtab_add_builtin("NDSolve`MethodOfLines", builtin_ndsolve_method);
    symtab_get_def("NDSolve`MethodOfLines")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("NDSolve`MethodOfLines",
        "NDSolve`MethodOfLines[eqns, u, {t, tmin, tmax}, {x, xmin, xmax}]\n"
        "\tsolves a partial differential equation by the method of lines: the\n"
        "\tspatial operator is discretized on a uniform grid (central finite\n"
        "\tdifferences) reducing the PDE to an ODE system integrated in time.\n"
        "\tThe grid resolution is set via Method -> {\"MethodOfLines\",\n"
        "\t\"SpatialDiscretization\" -> {\"TensorProductGrid\", \"MinPoints\" -> n}};\n"
        "\twithout MinPoints a 1-D grid is refined until it resolves the initial\n"
        "\tdata, up to \"MaxPoints\" (1537 at machine precision).");
}
