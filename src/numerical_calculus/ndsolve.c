/* Mathilda — NDSolve front-end: compile ODE initial-value problems into an
 * NdProblem (first-order reduced system), dispatch to a stepper, and assemble
 * the InterpolatingFunction result.  Numerical kernels live in the ndsolve_*
 * modules; see ndsolve_common.h. */
#include "ndsolve.h"
#include "ndsolve_common.h"
#include "ndsolve_compile.h"
#include "../sym_names.h"
#include "../sym_intern.h"
#include "../eval.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include "../numeric.h"
#include "../nc_accuracy.h"   /* shared AccuracyGoal/PrecisionGoal handling */
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
                o->spec.preserve_inexact = false;
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
    free(P->bind_y);
    free(P->fun_names);
    free(P->fun_state0);
}

/* ------------------------------------------------------------------ *
 *  Core solver                                                        *
 * ------------------------------------------------------------------ */
static void nd_warn(const char* tag, const char* msg) {
    fprintf(stderr, "NDSolve::%s: %s\n", tag, msg);
}

static Expr* ndsolve_core(Expr* res, const char* forced_method) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 3) return NULL;
    Expr** A = res->data.function.args;
    size_t argc = res->data.function.arg_count;

    /* ---- options (trailing) ---- */
    NdOpts o; nd_opts_default(&o);
    /* Default AccuracyGoal -> MachinePrecision (a ~15.95-digit absolute-error
     * floor); PrecisionGoal stays Automatic (-> WorkingPrecision/2).  Set before
     * option parsing so an explicit AccuracyGoal-> still overrides it. */
    o.acc_goal = NUMERIC_MACHINE_PRECISION_DIGITS;
    size_t pos_end = argc;
    while (pos_end > 3 && nd_is_option(A[pos_end - 1])) pos_end--;
    /* allow exactly 3 positional args; options after */
    for (size_t i = 3; i < argc; i++) if (nd_is_option(A[i])) nd_apply_option(A[i], &o);
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
    Expr** eqitems; size_t neq;
    if (head_is(eqns, SYM_List)) { neq = eqns->data.function.arg_count; eqitems = eqns->data.function.args; }
    else { neq = 1; eqitems = &eqns; }

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
        free(maxorder); free(is_ode); free(base); free(fnames);
        nd_warn("nodeqs", "could not identify a well-posed ODE system");
        return NULL;
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
    P.nfun = nfun; P.fun_applied = applied;
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
        if (nd_eval_to_double((Expr*)arg, o.spec, &pt)) { t0 = pt; t0_set = true; }
        double re, im = 0.0;
        if (nd_eval_to_double(val, o.spec, &re)) { /* real IC */ }
        else if (!nd_split_reim(val, o.spec, &re, &im)) continue;   /* complex IC */
        size_t gi = base[k] + (size_t)ord;
        P.Y0[gi] = re; Y0im[gi] = im; have_ic[gi] = true;
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

    Expr* result = NULL;
    if (build_ok && ic_ok) {
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
    free(maxorder); free(is_ode); free(base);
    nd_problem_free(&P);
    return result;   /* NULL leaves NDSolve[...] unevaluated */
}

/* ------------------------------------------------------------------ *
 *  Builtins                                                           *
 * ------------------------------------------------------------------ */
Expr* builtin_ndsolve(Expr* res) {
    return ndsolve_core(res, NULL);
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
    symtab_get_def("NDSolve")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;

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
        "NDSolve[eqns, u, {t, tmin, tmax}, {x, xmin, xmax}] solves a partial\n"
        "\tdifferential equation over a rectangular region by the method of lines,\n"
        "\tgiving a 2-D InterpolatingFunction applied as u[t, x].\n"
        "\tOptions: Method, WorkingPrecision, AccuracyGoal, PrecisionGoal,\n"
        "\tMaxSteps, MaxStepSize, MaxStepFraction, StartingStepSize,\n"
        "\tInterpolationOrder, StepMonitor, EvaluationMonitor.");

    /* Register one builtin per stepper as NDSolve`Name. */
    char buf[64];
    for (size_t i = 0; i < nd_stepper_count(); i++) {
        const NdStepper* S = nd_stepper_at(i);
        snprintf(buf, sizeof(buf), "NDSolve`%s", S->name);
        symtab_add_builtin(buf, builtin_ndsolve_method);
        symtab_get_def(buf)->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
        symtab_set_docstring(buf, S->doc);
    }
    /* Wolfram-facing method-name aliases (resolved by nd_lookup_stepper). */
    static const char* const aliases[] = {
        "ExplicitRungeKutta", "RungeKutta", "StiffnessSwitching"
    };
    for (size_t i = 0; i < sizeof(aliases)/sizeof(aliases[0]); i++) {
        snprintf(buf, sizeof(buf), "NDSolve`%s", aliases[i]);
        symtab_add_builtin(buf, builtin_ndsolve_method);
        symtab_get_def(buf)->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
        symtab_set_docstring(buf,
            "NDSolve`<method>[eqns, u, {x, xmin, xmax}] — solve forcing this "
            "time-integration method (see NDSolve).");
    }

    /* PDE controller: the Method-of-Lines spatial discretization. */
    symtab_add_builtin("NDSolve`MethodOfLines", builtin_ndsolve_method);
    symtab_get_def("NDSolve`MethodOfLines")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("NDSolve`MethodOfLines",
        "NDSolve`MethodOfLines[eqns, u, {t, tmin, tmax}, {x, xmin, xmax}]\n"
        "\tsolves a partial differential equation by the method of lines: the\n"
        "\tspatial operator is discretized on a uniform grid (central finite\n"
        "\tdifferences) reducing the PDE to an ODE system integrated in time.\n"
        "\tThe grid resolution is set via Method -> {\"MethodOfLines\",\n"
        "\t\"SpatialDiscretization\" -> {\"TensorProductGrid\", \"MinPoints\" -> n}}.");
}
