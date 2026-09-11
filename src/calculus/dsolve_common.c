/*
 * dsolve_common.c — shared substrate for the DSolve method cascade.
 *
 * See dsolve_common.h.  Contains: the problem parser (dsolve_parse), the
 * per-branch verify + condition-fit + result-assembly pipeline (dsolve_run),
 * and the small expression helpers the methods build on.  The design mirrors
 * NDSolve's held-equation parsing (nd_match_funcapp/nd_scan) and Solve's
 * verify-by-back-substitution policy, but keeps everything symbolic.
 */
#include "dsolve_common.h"
#include "dsolve.h"          /* g_dsolve_depth */

#include "../sym_names.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../eval.h"
#include "../common.h"
#include "../internal.h"
#include "../zero_test.h"
#include "../ndarray.h"
#include "../parse.h"           /* parse_expression: verify-residual normalisation rule */
#include "integrate.h"          /* g_integrate_quiet: silence speculative nonelem */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>               /* isnan / isfinite / NAN: numeric verify probe */

/* ds_has_undefined_function is public (declared in dsolve_common.h): kovacic's
 * forcing closure needs it to accept an arbitrary-forcing VoP integral. */

/* ------------------------------------------------------------------ *
 *  Small expression helpers                                           *
 * ------------------------------------------------------------------ */
Expr* ds_call1(const char* head, Expr* a) {
    return expr_new_function(expr_new_symbol(head), (Expr*[]){ a }, 1);
}
Expr* ds_call2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(head), (Expr*[]){ a, b }, 2);
}
Expr* ds_d(Expr* e, Expr* v)         { return eval_and_free(ds_call2(SYM_D, e, v)); }
Expr* ds_integrate(Expr* e, Expr* v) { return eval_and_free(ds_call2(SYM_Integrate, e, v)); }
Expr* ds_solve(Expr* eq, Expr* v)    { return eval_and_free(ds_call2(SYM_Solve, eq, v)); }

Expr* ds_subst(Expr* body, Expr* from, Expr* to) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule), (Expr*[]){ from, to }, 2);
    return eval_and_free(internal_replace_all((Expr*[]){ body, rule }, 2));
}

Expr* ds_make_funcapp(const char* fname, int order, const char* xvar) {
    if (order == 0) return ds_call1(fname, expr_new_symbol(xvar));
    Expr* d1 = expr_new_function(expr_new_symbol(SYM_Derivative),
                                 (Expr*[]){ expr_new_integer(order) }, 1);
    Expr* d2 = expr_new_function(d1, (Expr*[]){ expr_new_symbol(fname) }, 1);
    return expr_new_function(d2, (Expr*[]){ expr_new_symbol(xvar) }, 1);
}

Expr* ds_const(int k) {
    return expr_new_function(expr_new_symbol("C"), (Expr*[]){ expr_new_integer(k) }, 1);
}

bool ds_contains(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == name;
    if (e->type == EXPR_FUNCTION) {
        if (ds_contains(e->data.function.head, name)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (ds_contains(e->data.function.args[i], name)) return true;
    }
    return false;
}

bool ds_is_zero(const Expr* e)    { return zero_test_decide(e) == ZERO_TEST_TRUE; }
bool ds_is_nonzero(const Expr* e) { return zero_test_decide(e) == ZERO_TEST_FALSE; }
bool ds_has_head(const Expr* e, const char* head) { return ds_contains(e, head); }

bool ds_free_of(const Expr* e, const char* var) {
    if (!ds_contains(e, var)) return true;                 /* cheap syntactic pass */
    Expr* d = ds_d(expr_copy((Expr*)e), expr_new_symbol(var));
    bool zero = ds_is_zero(d);
    expr_free(d);
    return zero;
}

Expr* ds_simplify(Expr* e) { return eval_and_free(ds_call1("Simplify", e)); }

/* Materialize a packed list to a plain List (recursively over its direct
 * elements) so components can be walked via data.function.args.  Consumes e. */
Expr* ds_delist(Expr* e) {
    if (e && is_packed_list(e)) { Expr* p = ndarray_to_nested_list(e); expr_free(e); return p; }
    if (!e || e->type != EXPR_FUNCTION
        || e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol.name != SYM_List) return e;
    size_t n = e->data.function.arg_count; bool any = false;
    for (size_t i = 0; i < n; i++) if (is_packed_list(e->data.function.args[i])) any = true;
    if (!any) return e;
    Expr** args = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) args[i] = ds_delist(expr_copy(e->data.function.args[i]));
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), args, n); free(args); expr_free(e);
    return r;
}

/* ------------------------------------------------------------------ *
 *  Derivative matching / normalization (mirrors NDSolve, kept local)  *
 * ------------------------------------------------------------------ */
/* Match  y[arg] (order 0) or Derivative[m][y][arg]; fills interned fname,
 * order and the (borrowed) argument. */
static bool ds_match_funcapp(const Expr* e, const char** fname, int* order, const Expr** arg) {
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

/* Update maxorder[k] for each dependent function; flag has_x if any funcapp is
 * at the independent variable (=> an ODE term, not a boundary condition). */
static void ds_scan(const Expr* e, const char** funcs, size_t nfun, const char* xvar,
                    int* maxorder, bool* has_x) {
    if (!e) return;
    const char* fn; int ord; const Expr* arg;
    if (ds_match_funcapp(e, &fn, &ord, &arg)) {
        for (size_t k = 0; k < nfun; k++)
            if (fn == funcs[k]) {
                if (ord > maxorder[k]) maxorder[k] = ord;
                if (arg->type == EXPR_SYMBOL && arg->data.symbol.name == xvar) *has_x = true;
            }
    }
    if (e->type == EXPR_FUNCTION) {
        ds_scan(e->data.function.head, funcs, nfun, xvar, maxorder, has_x);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            ds_scan(e->data.function.args[i], funcs, nfun, xvar, maxorder, has_x);
    }
}

/* Return a fresh copy of `e` with every D[...] subtree evaluated into its
 * canonical Derivative[...] form (so methods pattern-match one shape). */
static Expr* ds_normalize_derivs(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_D)
        return eval_and_free(expr_copy((Expr*)e));
    Expr* nh = ds_normalize_derivs(h);
    size_t n = e->data.function.arg_count;
    Expr** args = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) args[i] = ds_normalize_derivs(e->data.function.args[i]);
    Expr* r = expr_new_function(nh, args, n);
    free(args);
    return r;
}

/* ------------------------------------------------------------------ *
 *  Generated-constant helpers (local copies of reduce_int's statics)  *
 * ------------------------------------------------------------------ */
static void ds_collect_consts(const Expr* e, Expr*** list, size_t* n) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && strcmp(h->data.symbol.name, "C") == 0
        && e->data.function.arg_count == 1) {
        for (size_t i = 0; i < *n; i++) if (expr_eq((*list)[i], e)) return;
        *list = realloc(*list, (*n + 1) * sizeof(Expr*));
        (*list)[(*n)++] = expr_copy((Expr*)e);
        return;
    }
    ds_collect_consts(h, list, n);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        ds_collect_consts(e->data.function.args[i], list, n);
}

/* Rewrite every generated constant C[k] to head[k]; fresh tree (no-op copy for
 * head == "C"). */
Expr* ds_rename_param(const Expr* e, const char* head) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && strcmp(h->data.symbol.name, "C") == 0
        && e->data.function.arg_count == 1)
        return expr_new_function(expr_new_symbol(head),
            (Expr*[]){ ds_rename_param(e->data.function.args[0], head) }, 1);
    Expr* nh = ds_rename_param(h, head);
    size_t n = e->data.function.arg_count;
    Expr** args = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) args[i] = ds_rename_param(e->data.function.args[i], head);
    Expr* r = expr_new_function(nh, args, n);
    free(args);
    return r;
}

/* ------------------------------------------------------------------ *
 *  Option parsing                                                     *
 * ------------------------------------------------------------------ */
static bool ds_is_option(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return false;
    if (h->data.symbol.name != SYM_Rule && h->data.symbol.name != SYM_RuleDelayed) return false;
    return e->data.function.arg_count == 2 && e->data.function.args[0]->type == EXPR_SYMBOL;
}

static void ds_apply_option(const Expr* opt, DSolveProblem* P) {
    const Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    const char* nm = lhs->data.symbol.name;
    if (nm == SYM_GeneratedParameters) {
        Expr* v = eval_and_free(expr_copy(rhs));
        if (v && v->type == EXPR_SYMBOL) P->param_head = v->data.symbol.name;
        expr_free(v);
    } else if (nm == SYM_Assumptions) {
        if (P->assumptions) expr_free(P->assumptions);
        P->assumptions = eval_and_free(expr_copy(rhs));
    } else if (nm == SYM_Method) {
        Expr* v = eval_and_free(expr_copy(rhs));
        if (v && v->type == EXPR_STRING) P->method = intern_symbol(v->data.string);
        else if (v && v->type == EXPR_SYMBOL) P->method = v->data.symbol.name;
        else if (head_is(v, SYM_List) && v->data.function.arg_count >= 1
                 && v->data.function.args[0]->type == EXPR_STRING)
            P->method = intern_symbol(v->data.function.args[0]->data.string);
        expr_free(v);
    } else if (nm == intern_symbol("IncludeSingularSolutions")) {
        Expr* v = eval_and_free(expr_copy(rhs));
        if (v && v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True)
            P->include_singular = true;
        expr_free(v);
    }
}

/* ------------------------------------------------------------------ *
 *  Condition detection                                                *
 * ------------------------------------------------------------------ */
/* A normalized item is a boundary/initial condition iff it is Equal[a,b] with
 * one side a funcapp of a dependent function whose argument is NOT the
 * independent variable (i.e. a fixed point). */
static bool ds_is_condition(const Expr* item, const char** funcs, size_t nfun, const char* xvar,
                            size_t* fi, int* order, const Expr** point, const Expr** value) {
    if (!head_is((Expr*)item, SYM_Equal) || item->data.function.arg_count != 2) return false;
    const Expr* s0 = item->data.function.args[0];
    const Expr* s1 = item->data.function.args[1];
    const char* fn; int ord; const Expr* arg;
    const Expr* fa = NULL; const Expr* val = NULL;
    if (ds_match_funcapp(s0, &fn, &ord, &arg)) { fa = s0; val = s1; }
    else if (ds_match_funcapp(s1, &fn, &ord, &arg)) { fa = s1; val = s0; }
    if (!fa) return false;
    if (arg->type == EXPR_SYMBOL && arg->data.symbol.name == xvar) return false; /* an ODE term */
    size_t k; bool ok = false;
    for (k = 0; k < nfun; k++) if (fn == funcs[k]) { ok = true; break; }
    if (!ok) return false;
    *fi = k; *order = ord; *point = arg; *value = val;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Parser                                                             *
 * ------------------------------------------------------------------ */
bool dsolve_parse(Expr* res, DSolveProblem* P) {
    memset(P, 0, sizeof(*P));
    P->param_head = "C";
    if (!res || res->type != EXPR_FUNCTION) return false;
    Expr** A = res->data.function.args;
    size_t argc = res->data.function.arg_count;
    if (argc < 3) return false;

    /* options */
    for (size_t i = 3; i < argc; i++) if (ds_is_option(A[i])) ds_apply_option(A[i], P);

    /* ---- independent variables (A[2]) ---- */
    Expr* vars = A[2];
    if (vars->type == EXPR_SYMBOL) {
        P->nind = 1;
        P->ind_names = malloc(sizeof(char*));
        P->ind_names[0] = vars->data.symbol.name;
    } else if (head_is(vars, SYM_List)) {
        size_t nv = vars->data.function.arg_count;
        Expr** vi = vars->data.function.args;
        if (nv == 3 && vi[0]->type == EXPR_SYMBOL) {
            /* {x, xmin, xmax} : single ODE variable with a range */
            P->nind = 1;
            P->ind_names = malloc(sizeof(char*));
            P->ind_names[0] = vi[0]->data.symbol.name;
            P->xmin = expr_copy(vi[1]);
            P->xmax = expr_copy(vi[2]);
        } else {
            /* list of independent variables (PDE when >= 2) */
            P->nind = nv;
            P->ind_names = malloc(nv * sizeof(char*));
            for (size_t i = 0; i < nv; i++) {
                Expr* it = vi[i];
                if (it->type == EXPR_SYMBOL) P->ind_names[i] = it->data.symbol.name;
                else if (head_is(it, SYM_List) && it->data.function.arg_count >= 1
                         && it->data.function.args[0]->type == EXPR_SYMBOL)
                    P->ind_names[i] = it->data.function.args[0]->data.symbol.name;
                else { free(P->ind_names); P->ind_names = NULL; return false; }
            }
            if (nv >= 2) P->is_pde = true;
        }
    } else {
        return false;
    }
    if (P->nind == 0) return false;
    const char* xvar = P->ind_names[0];

    /* ---- dependent functions (A[1]) ---- */
    Expr* funcs = A[1];
    Expr** fi; size_t nfun; bool own_fi = false;
    if (head_is(funcs, SYM_List)) { nfun = funcs->data.function.arg_count; fi = funcs->data.function.args; }
    else { nfun = 1; fi = malloc(sizeof(Expr*)); fi[0] = funcs; own_fi = true; }
    if (nfun == 0) { if (own_fi) free(fi); return false; }
    P->nfun = nfun;
    P->fun_names = malloc(nfun * sizeof(char*));
    for (size_t k = 0; k < nfun; k++) {
        Expr* it = fi[k];
        if (it->type == EXPR_SYMBOL) P->fun_names[k] = it->data.symbol.name;
        else if (it->type == EXPR_FUNCTION && it->data.function.head->type == EXPR_SYMBOL
                 && it->data.function.arg_count == 1) {
            P->fun_names[k] = it->data.function.head->data.symbol.name;
            P->applied = true;
        } else { if (own_fi) free(fi); free(P->fun_names); P->fun_names = NULL; return false; }
    }
    if (own_fi) free(fi);

    /* ---- equations / conditions (A[0]) ---- */
    Expr* eqns = A[0];
    Expr** items; size_t nitems; bool own_items = false;
    if (head_is(eqns, SYM_List)) { nitems = eqns->data.function.arg_count; items = eqns->data.function.args; }
    else { nitems = 1; items = malloc(sizeof(Expr*)); items[0] = eqns; own_items = true; }

    /* Normalize each item (D -> Derivative). */
    Expr** norm = malloc(nitems * sizeof(Expr*));
    for (size_t i = 0; i < nitems; i++) norm[i] = ds_normalize_derivs(items[i]);
    if (own_items) free(items);

    P->max_order = calloc(nfun ? nfun : 1, sizeof(int));
    P->eq_residuals = malloc((nitems ? nitems : 1) * sizeof(Expr*));
    P->conds = malloc((nitems ? nitems : 1) * sizeof(DSolveCond));
    P->neq = 0; P->ncond = 0;

    for (size_t i = 0; i < nitems; i++) {
        size_t cfi; int cord; const Expr* cpt; const Expr* cval;
        if (ds_is_condition(norm[i], P->fun_names, nfun, xvar, &cfi, &cord, &cpt, &cval)) {
            P->conds[P->ncond].fi = cfi;
            P->conds[P->ncond].order = cord;
            P->conds[P->ncond].point = expr_copy((Expr*)cpt);
            P->conds[P->ncond].value = expr_copy((Expr*)cval);
            P->ncond++;
        } else if (head_is(norm[i], SYM_Equal) && norm[i]->data.function.arg_count == 2) {
            bool hx = false;
            ds_scan(norm[i], P->fun_names, nfun, xvar, P->max_order, &hx);
            Expr* lhs = norm[i]->data.function.args[0];
            Expr* rhs = norm[i]->data.function.args[1];
            P->eq_residuals[P->neq++] =
                eval_and_free(ds_call2(SYM_Subtract, expr_copy(lhs), expr_copy(rhs)));
        }
        /* other item heads are ignored (left symbolic by declining later) */
        expr_free(norm[i]);
    }
    free(norm);

    if (P->neq == 0) return false;
    return true;
}

void dsolve_problem_free(DSolveProblem* P) {
    if (!P) return;
    free(P->fun_names);
    free(P->max_order);
    free(P->ind_names);
    if (P->xmin) expr_free(P->xmin);
    if (P->xmax) expr_free(P->xmax);
    for (size_t i = 0; i < P->neq; i++) expr_free(P->eq_residuals[i]);
    free(P->eq_residuals);
    for (size_t i = 0; i < P->ncond; i++) { expr_free(P->conds[i].point); expr_free(P->conds[i].value); }
    free(P->conds);
    if (P->assumptions) expr_free(P->assumptions);
    memset(P, 0, sizeof(*P));
}

/* ------------------------------------------------------------------ *
 *  Verify / fit / assemble                                            *
 * ------------------------------------------------------------------ */
/* Substitute the candidate body into each equation residual and require none to
 * be decidably non-zero.  Rather than substitute y -> Function[{x}, body] and
 * let the evaluator reduce Derivative[k][y][x], we replace each derivative term
 * Derivative[k][y][x] -> D[body, {x,k}] and y[x] -> body directly.  This is
 * equivalent for elementary bodies but essential for a SeriesData body, whose
 * pure-function derivative Derivative[k][Function[{x}, SeriesData]][x] the
 * evaluator does not reduce (it returns 0) — matching the PDE-verify workaround. */
/* True if e carries an unevaluated DEFINITE Integrate (Integrate[_, _List]) or a
 * distributional / piecewise head (DiracDelta / HeavisideTheta / UnitStep /
 * Piecewise): a Green's-function / impulse / step residual that zero_test cannot
 * decide and whose numeric probe cannot run (differentiating a Piecewise body
 * yields a boundary term the probe reads as an undefined function, so it bails;
 * and zero_test evaluates a piecewise condition on the wrong branch, so it can
 * return a SPURIOUS nonzero -- it rejected the correct step-forced answer of
 * y''+4y==Sin[t]-UnitStep[t-2Pi]Sin[t], §2.2.15-1497).  Such a branch is accepted
 * on its construction plus the producing method's own numeric verify (the
 * per-interval + per-IC pw_num_ok of DSolve`PiecewiseForcing), never driven
 * through zero_test. */
static bool ds_residual_is_distributional(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* hd = e->data.function.head;
    if (hd->type == EXPR_SYMBOL) {
        const char* hn = hd->data.symbol.name;
        if (hn == SYM_HeavisideTheta || hn == SYM_UnitStep || hn == SYM_Piecewise
            || strcmp(hn, "DiracDelta") == 0) return true;
        if (hn == SYM_Integrate && e->data.function.arg_count >= 2 &&
            e->data.function.args[1]->type == EXPR_FUNCTION &&
            e->data.function.args[1]->data.function.head->type == EXPR_SYMBOL &&
            e->data.function.args[1]->data.function.head->data.symbol.name == SYM_List)
            return true;
    }
    if (ds_residual_is_distributional(hd)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ds_residual_is_distributional(e->data.function.args[i])) return true;
    return false;
}

/* Collect argument-position free symbols (a symbolic parameter such as k, a, b
 * that a numeric probe must instantiate).  Recurses into arguments only, never a
 * function head, so Sin/Cos/Log/Exp/C/Integrate stay untouched. */
static void ds_collect_arg_syms(const Expr* e, const char** out, int* n, int cap) {
    if (!e || *n >= cap) return;
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < *n; i++) if (out[i] == e->data.symbol.name) return;
        out[(*n)++] = e->data.symbol.name;
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        ds_collect_arg_syms(e->data.function.args[i], out, n, cap);
}

/* NUMERIC back-substitution probe: is the (already body-substituted) residual
 * `sub` numerically ZERO at a spread of clean real sample points?  Instantiates
 * the generated constants C[1..8] and any free parameter symbols to distinct
 * generic reals, then samples |sub| off the origin.  Returns true only when
 * CONFIDENT (>= 3 real samples numericize, >= 3 are tiny, none is a real
 * clearly-nonzero value); a complex/non-finite sample (a Log/ArcTan/Sqrt branch
 * cut crossed by the sample point) is skipped, never counted as nonzero, so a
 * correct answer is never rejected here.
 *
 * This is the same numeric verify sf_num_ok / l2_num_ok use, lifted into the
 * generic scalar verify purely as a fast KEEP short-circuit: zero_test_decide's
 * precision ladder can climb for many seconds on a residual that IS zero but
 * carries Log branch cuts (the variation-of-parameters answer of y''+y==Tan[x],
 * y''+y==2 Sec[x/2] -- >8 s, so the harness times them out).  When this cannot
 * confirm zero it returns false and the branch falls through to the unchanged
 * symbolic zero_test, so the REJECT path (a genuinely nonzero residual) is
 * unaffected. */
/* Instantiate the generated constants C[1..8] and every argument-position free
 * parameter symbol of `R` to distinct generic reals (skipping the independent
 * variable and the named numeric constants), so `R` can be sampled numerically in
 * `xv`.  Consumes `R`, returns the substituted expression (owned).  Shared by the
 * numeric KEEP short-circuit (ds_residual_numeric_zero) and the numeric REJECT
 * filter (ds_branch_num_ok) so both instantiate identically. */
static Expr* ds_subst_generics(Expr* R, const char* xv) {
    for (int k = 1; k <= 8; k++)
        R = ds_subst(R, ds_const(k), expr_new_real(0.31 + 0.17 * (double)k));
    const char* skip[] = { xv, intern_symbol("E"), intern_symbol("Pi"),
        intern_symbol("I"), intern_symbol("EulerGamma"), intern_symbol("Degree"),
        intern_symbol("GoldenRatio"), intern_symbol("Catalan"),
        intern_symbol("Infinity") };
    const int nskip = (int)(sizeof(skip)/sizeof(skip[0]));
    const char* syms[64]; int ns = 0;
    ds_collect_arg_syms(R, syms, &ns, 64);
    int pi = 0;
    for (int i = 0; i < ns; i++) {
        bool sk = false;
        for (int j = 0; j < nskip; j++) if (syms[i] == skip[j]) { sk = true; break; }
        if (sk) continue;
        R = ds_subst(R, expr_new_symbol(syms[i]), expr_new_real(0.37 + 0.11 * (double)pi));
        pi++;
    }
    return R;
}

static bool ds_residual_numeric_zero(const Expr* sub0, const char* xv) {
    Expr* R = ds_subst_generics(expr_copy((Expr*)sub0), xv);
    /* An arbitrary-function or inert-Integrate residual (Bessel operator with an
     * undefined forcing g[x]) never numericizes -- decline the numeric route so
     * the symbolic keep-on-undecidable policy handles it. */
    if (ds_has_undefined_function(R) || ds_has_head(R, SYM_Integrate)) {
        expr_free(R); return false;
    }
    static const double xs[] = { 0.31, 0.53, 0.74, 1.13, 1.47, 1.92, 2.31, 2.68 };
    int nsmall = 0, nbig = 0, ngood = 0;
    for (int i = 0; i < 8; i++) {
        Expr* at = ds_subst(expr_copy(R), expr_new_symbol(xv), expr_new_real(xs[i]));
        Expr* mg = eval_and_free(ds_call1("Abs", eval_and_free(ds_call1("N", expr_copy(at)))));
        Expr* ig = eval_and_free(ds_call1("Abs",
                       eval_and_free(ds_call1("N", ds_call1("Im", at)))));
        double mag = (mg && mg->type == EXPR_REAL) ? mg->data.real
                   : (mg && mg->type == EXPR_INTEGER) ? (double)mg->data.integer : NAN;
        double imag = (ig && ig->type == EXPR_REAL) ? ig->data.real
                    : (ig && ig->type == EXPR_INTEGER) ? (double)ig->data.integer : NAN;
        expr_free(mg); expr_free(ig);
        if (isnan(mag) || !isfinite(mag)) continue;           /* couldn't numericize */
        if (mag < 1e-8) { nsmall++; ngood++; }
        else if (mag > 1e-4 && isfinite(imag) && imag < 1e-6) { nbig++; ngood++; }
        /* else: a complex sample (branch cut) or borderline -- skip, don't reject */
    }
    expr_free(R);
    return ngood >= 3 && nsmall >= 3 && nbig == 0;
}

/* Numeric REJECT filter (see dsolve_common.h): false only when the body's residual
 * against the original ODE is robustly, finitely NONZERO at a majority of clean real
 * sample points with NOTHING looking like zero.  Uses the complex-aware modulus
 * Abs[N[residual]] (matching the corpus prelude), so a spurious complex principal-root
 * branch — the wrong +/- of a cube/square root — is rejected, while dsolve_verify_body
 * (symbolic, keeps undecidable branch-cut residuals) would let it through.  Conservative:
 * a residual that is ~0, non-numericizable, distributional, or mixed (some samples ~0,
 * a partial-domain-valid branch such as 1+x^3 for x>0) is KEPT. */
bool ds_has_radical_power(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2
        && e->data.function.args[1]->type != EXPR_INTEGER)
        return true;                               /* Rational/Real/symbolic exponent */
    if (ds_has_radical_power(h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ds_has_radical_power(e->data.function.args[i])) return true;
    return false;
}

bool ds_branch_num_ok(const DSolveProblem* P, const Expr* body) {
    if (P->nfun != 1) return true;                 /* systems verified separately */
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    int maxord = P->max_order[0];
    for (size_t e = 0; e < P->neq; e++) {
        Expr* sub = expr_copy(P->eq_residuals[e]);
        for (int k = maxord; k >= 1; k--) {
            Expr* dk = expr_copy((Expr*)body);
            for (int i = 0; i < k; i++) dk = ds_d(dk, expr_new_symbol(xvar));
            sub = ds_subst(sub, ds_make_funcapp(yname, k, xvar), dk);
        }
        sub = ds_subst(sub, ds_make_funcapp(yname, 0, xvar), expr_copy((Expr*)body));
        /* Distributional / Green's-function residual: cannot numerically judge -- keep. */
        if (ds_residual_is_distributional(sub)) { expr_free(sub); continue; }
        Expr* R = ds_subst_generics(sub, xvar);    /* consumes sub */
        if (ds_has_undefined_function(R) || ds_has_head(R, SYM_Integrate)) {
            expr_free(R); continue;                /* non-numericizable -- keep */
        }
        static const double xs[] = { 0.31, 0.53, 0.74, 1.13, 1.47, 1.92, 2.31, 2.68 };
        int nsmall = 0, nbig = 0, ngood = 0;
        for (int i = 0; i < 8; i++) {
            Expr* at = ds_subst(expr_copy(R), expr_new_symbol(xvar), expr_new_real(xs[i]));
            /* Abs[N[.]] is the complex modulus, so a complex branch value counts as
             * nonzero here (unlike the KEEP short-circuit, which skips complex). */
            Expr* mg = eval_and_free(ds_call1("Abs", eval_and_free(ds_call1("N", at))));
            double mag = (mg && mg->type == EXPR_REAL)    ? mg->data.real
                       : (mg && mg->type == EXPR_INTEGER) ? (double)mg->data.integer : NAN;
            expr_free(mg);
            if (isnan(mag) || !isfinite(mag)) continue;    /* couldn't numericize */
            ngood++;
            if (mag < 1e-8) nsmall++;
            else if (mag > 1e-4) nbig++;
        }
        expr_free(R);
        /* Reject only when CONFIDENT: robustly nonzero at every point it numericized
         * and no sample looks like zero (a mixed profile is a partial-domain branch). */
        if (ngood >= 4 && nbig >= 4 && nsmall == 0) return false;
    }
    return true;
}

/* Would the corpus harness (DSolve_test_status/dsolve_corpus_prelude.m) score this
 * FITTED first-order scalar body as verifying?  Samples the residual at the prelude's
 * own grid (11/10 + k*5/13, k=0..5) and applies its majority rule (OK iff the residual
 * is small at >= Ceiling[n/2] of the numericizing points).  Returns false only for the
 * prelude's "BAD".  Restricted to first-order scalar ODEs, so it never touches the
 * large-eigenvalue catastrophic-cancellation cases (2nd-order / systems) the prelude
 * rescues at 200-digit precision.  Because it reuses the prelude's grid+majority, a
 * branch the harness would PASS (>= half the points small) is NEVER dropped -- it only
 * removes would-be-FAIL explicit branches: a wrong fitted root, or a closed form valid
 * only on a sub-interval the fixed grid overshoots (a sqrt/cube-root IVP whose principal
 * branch flips past a pole -- §2.2.17-1636), so the cascade can fall through to a later
 * method's verifiable (often implicit) form. */
static bool ds_branch_corpus_verifiable(const DSolveProblem* P, const Expr* body) {
    if (P->nfun != 1 || P->max_order[0] != 1) return true;   /* first-order scalar only */
    if (!ds_has_radical_power(body)) return true;            /* only a radical can hide a
                                                             * spurious principal-root /
                                                             * pole-crossing branch */
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    for (size_t e = 0; e < P->neq; e++) {
        Expr* sub = expr_copy(P->eq_residuals[e]);
        Expr* dk = ds_d(expr_copy((Expr*)body), expr_new_symbol(xvar));
        sub = ds_subst(sub, ds_make_funcapp(yname, 1, xvar), dk);
        sub = ds_subst(sub, ds_make_funcapp(yname, 0, xvar), expr_copy((Expr*)body));
        if (ds_residual_is_distributional(sub)) { expr_free(sub); continue; }
        Expr* R = ds_subst_generics(sub, xvar);              /* consumes sub */
        if (ds_has_undefined_function(R) || ds_has_head(R, SYM_Integrate)) {
            expr_free(R); continue;                          /* non-numericizable -> keep */
        }
        int nsmall = 0, nnum = 0;
        for (int k = 0; k < 6; k++) {
            double xk = 1.1 + (double)k * (5.0 / 13.0);      /* the prelude's sweep grid */
            Expr* at = ds_subst(expr_copy(R), expr_new_symbol(xvar), expr_new_real(xk));
            Expr* mg = eval_and_free(ds_call1("Abs", eval_and_free(ds_call1("N", at))));
            double mag = (mg && mg->type == EXPR_REAL)    ? mg->data.real
                       : (mg && mg->type == EXPR_INTEGER) ? (double)mg->data.integer : NAN;
            expr_free(mg);
            if (isnan(mag) || !isfinite(mag)) continue;      /* prelude skips these too */
            nnum++;
            if (mag < 1e-6) nsmall++;                        /* $dsTol */
        }
        expr_free(R);
        if (nnum >= 2 && nsmall < (nnum + 1) / 2) return false;   /* prelude verdict "BAD" */
    }
    return true;
}

static bool dsolve_verify_body(const DSolveProblem* P, const Expr* body) {
    if (P->nfun != 1) return true;    /* systems: verified separately */
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];
    int maxord = P->max_order[0];
    for (size_t e = 0; e < P->neq; e++) {
        Expr* sub = expr_copy(P->eq_residuals[e]);
        for (int k = maxord; k >= 1; k--) {
            Expr* dk = expr_copy((Expr*)body);
            for (int i = 0; i < k; i++) dk = ds_d(dk, expr_new_symbol(xvar));
            sub = ds_subst(sub, ds_make_funcapp(yname, k, xvar), dk);
        }
        sub = ds_subst(sub, ds_make_funcapp(yname, 0, xvar), expr_copy((Expr*)body));
        /* Distributional / Green's-function residual (definite Integrate,
         * DiracDelta, HeavisideTheta): zero_test cannot decide it and its ladder
         * could spin, so keep the branch on construction rather than reject. */
        if (ds_residual_is_distributional(sub)) { expr_free(sub); continue; }
        /* A Gaussian x Erf residual — the integrating-factor solution of an exact
         * ODE, e.g. y''+x y'+y==0 whose closed form carries Erf[-I x/Sqrt[2]] —
         * used to defeat zero_test's numeric precision ladder (the E^(-x^2/2) of
         * the solution and the E^(x^2/2) from differentiating Erf sat in separate
         * summands, never combined to E^0, and numericalised as tiny*huge, so the
         * ladder climbed to 1000 bits on every sample and effectively hung).
         * zero_test_decide now performs the exponential-combining ExpandAll
         * normalisation itself (POSSIBLE_ZEROQ_IMPROVEMENTS.md #1), so no
         * Erf-gated pre-pass or FALSE-path re-check is needed here. */
        /* Fast KEEP when the residual is NUMERICALLY zero: zero_test's precision
         * ladder can spin for >8 s on a residual that is zero but carries Log/
         * ArcTan branch cuts (the variation-of-parameters answer of y''+y==Tan[x],
         * §2.2.14-1337/1341), so a numerically-confirmed zero keeps the branch
         * without the symbolic test.  Only confirms zero; a nonzero / undecidable
         * residual falls through to zero_test unchanged (reject path preserved). */
        if (ds_residual_numeric_zero(sub, xvar)) { expr_free(sub); continue; }
        ZeroTestResult zt = zero_test_decide(sub);
        if (zt == ZERO_TEST_FALSE) { expr_free(sub); return false; }
        expr_free(sub);
    }
    return true;
}

/* How the constant-fit turned out, so dsolve_run can decide per branch WITH
 * knowledge of the sibling branches (a decision no single branch can make alone). */
enum { FIT_OK = 0,     /* Solve fixed >=1 constant (fully fit, or under-determined). */
       FIT_EMPTY = 1,  /* scalar Solve returned {} -- the condition is unsatisfiable
                        * on this branch (a wrong +/- sign / Root index) IF a sibling
                        * fits, else a genuine basis singularity to keep. */
       FIT_UNDEF = 2,  /* the fit substituted Undefined/$Failed -- never a solution. */
       FIT_UNDECIDED = 3 };/* Solve bubbled back unevaluated with a condition present --
                        * the branch still carries its generated constant, so this
                        * method did NOT solve the IVP (distinct from an under-determined
                        * fit, where Solve SUCCEEDS and the free constant is genuine). */

/* Fit generated constants to the initial/boundary conditions; returns a fresh
 * body (the general body copied when there is nothing to fit).  Sets *no_solution
 * (when non-NULL) true only when Solve PROVES the conditions inconsistent (the LIST
 * form returns {} -- an over-determined BVP with no solution).  Sets *fit_state to
 * FIT_OK / FIT_EMPTY / FIT_UNDEF (see the enum) so the caller can drop a branch that
 * an initial condition cannot be met on while keeping its siblings.  An undecided
 * scalar fit (Solve stays unevaluated) keeps the general solution as FIT_OK. */
static Expr* dsolve_fit_constants(const DSolveProblem* P, const Expr* body,
                                  bool* no_solution, int* fit_state) {
    if (no_solution) *no_solution = false;
    if (fit_state) *fit_state = FIT_OK;
    if (P->ncond == 0) return expr_copy((Expr*)body);
    Expr** params = NULL; size_t npar = 0;
    ds_collect_consts(body, &params, &npar);
    if (npar == 0) return expr_copy((Expr*)body);
    const char* xvar = P->ind_names[0];

    Expr** eqs = malloc(P->ncond * sizeof(Expr*));
    size_t neq = 0;
    /* A SeriesData body does not reduce under x->point (SeriesData[x0,x0,{a...}]
     * stays a SeriesData, so the fit equation never exposes a[0]=C[1], a[1]=C[2],
     * ...).  Take its Normal (the truncated polynomial) for the fit equations,
     * which evaluates and differentiates at the IC point normally -- the
     * Frobenius/ordinary-point series IVP (§2.2.14 1381/1384/1385).  The FINAL
     * body keeps its SeriesData form (line ~616 substitutes the fitted constants
     * into the original body); Normal is the identity on an elementary body. */
    const Expr* fitbody = body;
    Expr* normbody = NULL;
    if (ds_has_head(body, SYM_SeriesData)) {
        normbody = eval_and_free(ds_call1("Normal", expr_copy((Expr*)body)));
        fitbody = normbody;
    }
    for (size_t c = 0; c < P->ncond; c++) {
        if (P->conds[c].fi != 0) continue;               /* single-function M0 */
        Expr* bexpr = expr_copy((Expr*)fitbody);
        for (int d = 0; d < P->conds[c].order; d++)
            bexpr = ds_d(bexpr, expr_new_symbol(xvar));
        bexpr = ds_subst(bexpr, expr_new_symbol(xvar), expr_copy(P->conds[c].point));
        eqs[neq++] = expr_new_function(expr_new_symbol(SYM_Equal),
                        (Expr*[]){ bexpr, eval_and_free(expr_copy(P->conds[c].value)) }, 2);
    }
    if (normbody) expr_free(normbody);
    /* A single condition fitting a single constant is solved in Solve's SCALAR
     * form (Solve[eq, C[1]]), never the list form (Solve[{eq}, {C[1]}]): only the
     * scalar form applies inverse-function inversion, so a constant sitting inside
     * a transcendental -- Sqrt[C[1]]==1, Log[3+C[1]]==0, Tan[..C[1]..]==k -- fits
     * (the common first-order IVP), where the list form bubbles back unevaluated
     * and leaked the general solution's C[1].  Multi-condition (2nd-order IVP) or
     * under-determined fits keep the list form, which handles the linear system. */
    bool scalar_fit = (neq == 1 && npar == 1);
    Expr* solres;
    if (scalar_fit) {
        solres = ds_solve(eqs[0], params[0]);
        free(eqs); free(params);
    } else {
        Expr* eqlist = expr_new_function(expr_new_symbol(SYM_List), eqs, neq);
        free(eqs);
        Expr* varlist = expr_new_function(expr_new_symbol(SYM_List), params, npar);
        free(params);
        solres = ds_solve(eqlist, varlist);
    }
    Expr* fitted = NULL;
    bool inconsistent = false;
    bool scalar_empty = false;
    if (solres && head_is(solres, SYM_List)) {
        if (solres->data.function.arg_count == 0) {
            /* Empty Solve result.  From the multi-condition LIST form this is a
             * genuine inconsistency -> no solution (BVP soundness, M11).  From the
             * single-condition SCALAR form it is EITHER an unsatisfiable inverse
             * branch (a wrong +/- sign / Root index whose sibling fits -- drop) OR a
             * basis SINGULARITY (a lone branch fitted at an infinite point -- keep).
             * dsolve_fit_constants cannot tell which alone, so it reports FIT_EMPTY
             * and dsolve_run decides using the sibling branches. */
            if (!scalar_fit) inconsistent = true;
            else scalar_empty = true;
        } else {
            /* Multiple solution branches can return (a +/- root, several fitted-constant
             * values): choose the FIRST whose fitted body numerically satisfies the ODE,
             * so a spurious root that meets the initial condition but NOT the equation is
             * skipped -- the Bernoulli sqrt-y IVP y'-2y==2Sqrt[y], y(0)=1 solves
             * (C[1]e^x-1)^2 and Solve[(C[1]-1)^2==1,C[1]] returns BOTH C[1]->0 (=> y=1,
             * which fails the ODE) and C[1]->2 (correct); the old args[0]-only pick could
             * take the wrong one.  Falls back to args[0] when none is confidently verified
             * (numeric check inconclusive), preserving prior single-branch behavior. */
            size_t nbr = solres->data.function.arg_count;
            Expr* branch = solres->data.function.args[0];   /* List[Rule[C[k],val],...] */
            for (size_t bi = 0; bi < nbr && nbr > 1; bi++) {
                Expr* cb = solres->data.function.args[bi];
                if (!head_is(cb, SYM_List)) continue;
                Expr* cand = eval_and_free(internal_replace_all(
                    (Expr*[]){ expr_copy((Expr*)body), expr_copy(cb) }, 2));
                bool ok = ds_branch_num_ok(P, cand);
                expr_free(cand);
                if (ok) { branch = cb; break; }
            }
            if (head_is(branch, SYM_List)) {
                fitted = eval_and_free(internal_replace_all(
                    (Expr*[]){ expr_copy((Expr*)body), expr_copy(branch) }, 2));
                /* A multivalued inverse (Solve[Tan[..C..]==k, C]) fits as a
                 * ConditionalExpression over an integer family, reintroducing the
                 * generated constant.  For a fully-determined single-constant fit
                 * every branch satisfies the condition, so take the PRINCIPAL one:
                 * strip ConditionalExpression, then collapse the residual family
                 * index C[_] -> 0.  Guarded to the ConditionalExpression case so a
                 * clean fit is untouched and an unfitted constant is never zeroed. */
                if (scalar_fit && fitted &&
                    ds_has_head(fitted, SYM_ConditionalExpression)) {
                    Expr* strip = parse_expression("ConditionalExpression[e_, _] :> e");
                    fitted = eval_and_free(internal_replace_all(
                        (Expr*[]){ fitted, strip }, 2));
                    Expr* zero = parse_expression("C[_] -> 0");
                    fitted = eval_and_free(internal_replace_all(
                        (Expr*[]){ fitted, zero }, 2));
                }
            }
        }
    }
    if (solres) expr_free(solres);
    (void)neq;
    if (inconsistent) {            /* over-determined BVP: no solution */
        if (no_solution) *no_solution = true;
        return expr_copy((Expr*)body);
    }
    bool applied = (fitted != NULL);                /* Solve produced a fit */
    if (!fitted) fitted = expr_copy((Expr*)body);   /* undecided / singularity: keep general */

    /* Classify the fit for dsolve_run's per-branch keep/drop decision.  An
     * Undefined/$Failed value (the fit had no consistent constant on this branch --
     * 2.2.12-1147's wrong answer) is never a solution.  A scalar empty-Solve is an
     * unsatisfiable branch (drop if a sibling fits) or a singularity (keep if not).
     * A condition present but Solve bubbled back unevaluated (no fit applied) means
     * the branch still carries its generated constant -- the IVP is UNSOLVED by this
     * method (FIT_UNDECIDED); a later cascade method may fit it (e.g. Homogeneous's
     * transcendental log-form does not invert for the constant, but Exact's
     * polynomial first integral does).  This is DISTINCT from an under-determined
     * fit, where Solve SUCCEEDS and the leftover free constant is genuine
     * (y''+y==0, y[0]==0, y[Pi]==0 -> C[2] Sin[x], FIT_OK). */
    if (fit_state) {
        if (ds_contains(fitted, SYM_Undefined) || ds_contains(fitted, intern_symbol("$Failed")))
            *fit_state = FIT_UNDEF;
        else if (scalar_empty)
            *fit_state = FIT_EMPTY;
        else if (!applied && P->ncond > 0)
            *fit_state = FIT_UNDECIDED;
        else
            *fit_state = FIT_OK;
    }
    return fitted;
}

/* Wrap finished branch bodies into {{u -> Function[{x}, b]}} / {{u[x] -> b}}. */
static Expr* dsolve_assemble(const DSolveProblem* P, Expr** bodies, size_t nb) {
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];
    Expr** branches = malloc(nb * sizeof(Expr*));
    for (size_t b = 0; b < nb; b++) {
        Expr* body = ds_rename_param(bodies[b], P->param_head);
        Expr* lhs; Expr* rhs;
        if (P->applied) {
            lhs = expr_new_function(expr_new_symbol(yname),
                                    (Expr*[]){ expr_new_symbol(xvar) }, 1);
            rhs = body;
        } else {
            lhs = expr_new_symbol(yname);
            Expr* plist = expr_new_function(expr_new_symbol(SYM_List),
                                            (Expr*[]){ expr_new_symbol(xvar) }, 1);
            rhs = expr_new_function(expr_new_symbol(SYM_Function),
                                    (Expr*[]){ plist, body }, 2);
        }
        Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule), (Expr*[]){ lhs, rhs }, 2);
        branches[b] = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ rule }, 1);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), branches, nb);
    free(branches);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Root analysis + variation of parameters (shared by the linear ODE  *
 *  methods: constant-coefficient and Euler-Cauchy)                    *
 * ------------------------------------------------------------------ */
static int poly_root_multiplicity(const Expr* poly, const char* var, const Expr* root, int cap) {
    int m = 0;
    Expr* d = expr_copy((Expr*)poly);
    while (m <= cap) {
        Expr* at = ds_subst(expr_copy(d), expr_new_symbol(var), expr_copy((Expr*)root));
        bool z = ds_is_zero(at);
        expr_free(at);
        if (!z) break;
        m++;
        d = ds_d(d, expr_new_symbol(var));
    }
    expr_free(d);
    return m;
}

bool dsolve_analyze_roots(const Expr* poly, const char* var, int degree, DSolveRoots* out) {
    memset(out, 0, sizeof(*out));
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                   (Expr*[]){ expr_copy((Expr*)poly), expr_new_integer(0) }, 2);
    Expr* sol = ds_solve(eq, expr_new_symbol(var));
    size_t nd = 0;
    Expr** roots = dsolve_extract_solutions(sol, var, &nd);
    if (sol) expr_free(sol);
    if (!roots) return false;

    Expr** dist = malloc(nd * sizeof(Expr*));
    size_t ndist = 0;
    for (size_t i = 0; i < nd; i++) {
        bool dup = false;
        for (size_t j = 0; j < ndist && !dup; j++) {
            Expr* d = eval_and_free(ds_call2(SYM_Subtract, expr_copy(roots[i]), expr_copy(dist[j])));
            if (ds_is_zero(d)) dup = true;
            expr_free(d);
        }
        if (dup) expr_free(roots[i]); else dist[ndist++] = roots[i];
    }
    free(roots);

    out->roots = dist; out->ndist = ndist;
    out->mult = malloc(ndist * sizeof(int));
    out->im = malloc(ndist * sizeof(Expr*));
    out->isreal = malloc(ndist * sizeof(bool));
    out->total = 0;
    for (size_t i = 0; i < ndist; i++) {
        out->mult[i] = poly_root_multiplicity(poly, var, dist[i], degree);
        out->total += out->mult[i];
        out->im[i] = eval_and_free(ds_call1("Im", expr_copy(dist[i])));
        out->isreal[i] = ds_is_zero(out->im[i]);
    }
    return true;
}

bool dsolve_linear_coeffs(DSolveProblem* P, Expr*** coeffs, Expr** forcing, int* order) {
    if (P->nfun != 1 || P->neq != 1) return false;
    int n = P->max_order[0];
    if (n < 1) return false;
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];

    Expr* Ralg = expr_copy(P->eq_residuals[0]);
    const char** Dn = malloc((size_t)(n + 1) * sizeof(char*));
    for (int k = 0; k <= n; k++) {
        char buf[32]; snprintf(buf, sizeof(buf), "DSolve`D%d", k);
        Dn[k] = intern_symbol(buf);
        Ralg = ds_subst(Ralg, ds_make_funcapp(yname, k, xvar), expr_new_symbol(Dn[k]));
    }
    /* c_k = R_{D_k}; linear iff each is free of every D_j */
    Expr** c = malloc((size_t)(n + 1) * sizeof(Expr*));
    bool ok = true;
    for (int k = 0; k <= n; k++) {
        c[k] = ds_d(expr_copy(Ralg), expr_new_symbol(Dn[k]));
        if (ok) for (int j = 0; j <= n && ok; j++) if (!ds_free_of(c[k], Dn[j])) ok = false;
    }
    Expr* R0 = NULL; Expr* g = NULL;
    if (ok) {
        R0 = expr_copy(Ralg);
        for (int k = 0; k <= n; k++) R0 = ds_subst(R0, expr_new_symbol(Dn[k]), expr_new_integer(0));
        g = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(R0)));
        /* verify R == Σ c_k D_k + R0 */
        Expr** t = malloc((size_t)(n + 1) * sizeof(Expr*));
        for (int k = 0; k <= n; k++)
            t[k] = ds_call2(SYM_Times, expr_copy(c[k]), expr_new_symbol(Dn[k]));
        Expr* lin = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), t, (size_t)(n + 1)));
        free(t);
        Expr* recon = eval_and_free(ds_call2(SYM_Plus, lin, expr_copy(R0)));
        Expr* chk = eval_and_free(ds_call2(SYM_Subtract, expr_copy(Ralg), recon));
        ok = ds_is_zero(chk);
        expr_free(chk);
    }
    expr_free(Ralg); free(Dn);
    if (R0) expr_free(R0);
    if (!ok) { for (int k = 0; k <= n; k++) expr_free(c[k]); free(c); if (g) expr_free(g); return false; }
    *coeffs = c; *forcing = g; *order = n;
    return true;
}

/* True iff e is a rational function of xvar: after Together, both its Numerator
 * and Denominator are polynomials in xvar.  Gates the polynomial/rational
 * normaliser below away from transcendental-in-x coefficients (E^(f(x)), Sin[x],
 * …), which FLINT's rational canonicaliser cannot process — see the note in
 * dsolve_linear_normalize.  e is borrowed. */
bool ds_is_rational_in(const Expr* e, const char* xvar) {
    Expr* tg  = eval_and_free(ds_call1("Together", expr_copy((Expr*)e)));
    Expr* num = eval_and_free(ds_call1("Numerator", expr_copy(tg)));
    Expr* den = eval_and_free(ds_call1("Denominator", tg));            /* consumes tg */
    Expr* pn  = eval_and_free(ds_call2("PolynomialQ", num, expr_new_symbol(xvar)));
    Expr* pd  = eval_and_free(ds_call2("PolynomialQ", den, expr_new_symbol(xvar)));
    bool ok = pn && pn->type == EXPR_SYMBOL && pn->data.symbol.name == SYM_True
           && pd && pd->type == EXPR_SYMBOL && pd->data.symbol.name == SYM_True;
    if (pn) expr_free(pn);
    if (pd) expr_free(pd);
    return ok;
}

/* Normalize an extracted linear-ODE coefficient vector for the constant-
 * coefficient / Euler detectors: (1) clear denominators by multiplying c[] and g
 * through by the product of their denominators, so the coefficients become
 * polynomials — an equation solved for the top derivative with a rational RHS
 * (y''' == (24x+24y)/x^3, giving c_0 = -24/x^3) becomes the polynomial-coefficient
 * Euler form x^3 y''' - 24 y == 24 x; then (2) divide through by the polynomial
 * GCD of the coefficients, so a common factor is removed — x(y'''+2y''-y'-2y)==1
 * becomes the CONSTANT-coefficient y'''+2y''-y'-2y == 1/x.  Value-preserving; the
 * forcing g is scaled the same way (it may acquire the 1/x).  Applied only by the
 * const-coeff / Euler methods (NOT the shared extractor) so the OperatorFactor /
 * SymmetricSquare factorization searches keep the coefficient form they expect.
 * c[] (length n+1) and *g are updated in place. */
void dsolve_linear_normalize(Expr** c, Expr** g, int n, const char* xvar) {
    /* Only the user's (possibly non-canonical) direct input needs normalizing;
     * equations generated by internal recursions (OperatorFactor peeling a
     * quotient, Riccati's linearization, …) are already canonical, and paying the
     * Denominator/PolynomialGCD/Cancel cost on every deep recursion is what pushed
     * the OperatorFactor stress corpus over its time budget.  g_dsolve_depth == 1
     * at the outermost call, >= 2 inside a recursion. */
    if (g_dsolve_depth > 1) return;
    /* Gate: normalize's Denominator / PolynomialGCD / Cancel machinery is a
     * polynomial/rational-function transform in x.  A coefficient that is
     * transcendental in x (e.g. a_0 = a(lam E^(lam x) - a E^(2 lam x)), or a
     * denominator (1 + E^(x^2/2))^2) is NOT a rational function of x; feeding it
     * to PolynomialGCD/Cancel makes the FLINT rational-canonicaliser emit a
     * non-finite content and then recurse on GCD(-Infinity, 1) until the stack
     * overflows (a SIGSEGV on Kamke exp-coefficient equations 876/879, and the
     * (1+E^(x^2/2))^2 case 298).  Such equations are never constant-coefficient
     * or Euler anyway, so skipping normalize loses nothing and lets them fall
     * through to Kovacic / Frobenius.  Guard by requiring every coefficient and
     * the forcing to be rational in x (numerator and denominator polynomial). */
    for (int k = 0; k <= n; k++)
        if (!ds_is_rational_in(c[k], xvar)) return;
    if (!ds_is_rational_in(*g, xvar)) return;
    /* (1) clear denominators: mult = Π Denominator[c_k] · Denominator[g] */
    Expr* mult = expr_new_integer(1);
    for (int k = 0; k <= n; k++)
        mult = eval_and_free(ds_call2(SYM_Times, mult,
                   eval_and_free(ds_call1("Denominator", expr_copy(c[k])))));
    mult = eval_and_free(ds_call2(SYM_Times, mult,
               eval_and_free(ds_call1("Denominator", expr_copy(*g)))));
    if (!ds_free_of(mult, xvar)) {
        for (int k = 0; k <= n; k++)
            c[k] = eval_and_free(ds_call1("Cancel",
                       ds_call2(SYM_Times, c[k], expr_copy(mult))));
        *g = eval_and_free(ds_call1("Cancel",
                 ds_call2(SYM_Times, *g, expr_copy(mult))));
    }
    expr_free(mult);
    /* (2) divide by the polynomial GCD of the nonzero coefficients */
    Expr* d = NULL;
    for (int k = 0; k <= n; k++) {
        if (ds_is_zero(c[k])) continue;
        d = d ? eval_and_free(ds_call2("PolynomialGCD", d, expr_copy(c[k])))
              : expr_copy(c[k]);
    }
    if (d && !ds_free_of(d, xvar) && !ds_is_zero(d)) {
        for (int k = 0; k <= n; k++)
            c[k] = eval_and_free(ds_call1("Cancel",
                       ds_call2(SYM_Times, c[k],
                           ds_call2(SYM_Power, expr_copy(d), expr_new_integer(-1)))));
        *g = eval_and_free(ds_call1("Cancel",
                 ds_call2(SYM_Times, *g,
                     ds_call2(SYM_Power, expr_copy(d), expr_new_integer(-1)))));
    }
    if (d) expr_free(d);
}

static bool second_order_PQ_impl(DSolveProblem* P, Expr** Pc, Expr** Qc,
                                 bool require_homog) {
    if (P->nfun != 1 || P->neq != 1) return false;
    if (P->max_order[0] != 2) return false;
    Expr** c; Expr* g; int n;
    if (!dsolve_linear_coeffs(P, &c, &g, &n)) return false;
    bool homog = ds_is_zero(g);
    expr_free(g);
    if (n != 2 || (require_homog && !homog)) {
        for (int k = 0; k <= n; k++) expr_free(c[k]);
        free(c);
        return false;
    }
    /* normalized P = c1/c2, Q = c0/c2 (from the homogeneous part; any forcing is
     * the caller's responsibility) */
    *Pc = ds_simplify(ds_call2(SYM_Times, expr_copy(c[1]),
              expr_new_function(expr_new_symbol(SYM_Power),
                  (Expr*[]){ expr_copy(c[2]), expr_new_integer(-1) }, 2)));
    *Qc = ds_simplify(ds_call2(SYM_Times, expr_copy(c[0]),
              expr_new_function(expr_new_symbol(SYM_Power),
                  (Expr*[]){ expr_copy(c[2]), expr_new_integer(-1) }, 2)));
    for (int k = 0; k <= 2; k++) expr_free(c[k]);
    free(c);
    return true;
}

bool dsolve_second_order_PQ(DSolveProblem* P, Expr** Pc, Expr** Qc) {
    return second_order_PQ_impl(P, Pc, Qc, true);
}

/* As dsolve_second_order_PQ, but accepts an INHOMOGENEOUS equation too, returning
 * P, Q of the homogeneous part.  The caller (Kovacic) then adds the particular
 * solution by variation of parameters over the fundamental set it recovers. */
bool dsolve_second_order_PQ_forced(DSolveProblem* P, Expr** Pc, Expr** Qc) {
    return second_order_PQ_impl(P, Pc, Qc, false);
}

Expr* dsolve_normal_form(const Expr* Pc, const Expr* Qc, const char* xvar,
                         Expr** recovery_out) {
    /* r = P^2/4 + P'/2 - Q = (P^2 + 2 P' - 4 Q) / 4 */
    Expr* Psq = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy((Expr*)Pc), expr_new_integer(2) }, 2));
    Expr* dP  = ds_d(expr_copy((Expr*)Pc), expr_new_symbol(xvar));
    Expr* num = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){
                    Psq,
                    ds_call2(SYM_Times, expr_new_integer(2), dP),
                    ds_call2(SYM_Times, expr_new_integer(-4), expr_copy((Expr*)Qc))
                }, 3));
    Expr* r = ds_simplify(ds_call2(SYM_Times, num,
                  expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_new_integer(4), expr_new_integer(-1) }, 2)));

    if (recovery_out) {
        *recovery_out = NULL;
        /* recovery w = Exp[-Integrate[P/2, x]]; guard D[Integrate[P/2]] == P/2. */
        Expr* half = eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)Pc),
                         expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ expr_new_integer(2), expr_new_integer(-1) }, 2)));
        g_integrate_quiet++;   /* speculative recovery ∫P/2; non-elementary => decline */
        Expr* integ = ds_integrate(expr_copy(half), expr_new_symbol(xvar));
        g_integrate_quiet--;
        if (!ds_has_head(integ, SYM_Integrate)) {
            Expr* back = ds_d(expr_copy(integ), expr_new_symbol(xvar));
            Expr* diff = eval_and_free(ds_call2(SYM_Subtract, back, expr_copy(half)));
            if (ds_is_zero(diff))
                *recovery_out = eval_and_free(ds_call1("Exp",
                    eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(integ)))));
            expr_free(diff);
        }
        expr_free(half); expr_free(integ);
    }
    return r;
}

void dsolve_roots_free(DSolveRoots* r) {
    if (!r->roots) return;
    for (size_t i = 0; i < r->ndist; i++) { expr_free(r->roots[i]); expr_free(r->im[i]); }
    free(r->roots); free(r->mult); free(r->im); free(r->isreal);
    memset(r, 0, sizeof(*r));
}

/* n x n matrix List of dv[k][j], with column `repl` (if >= 0) set to the vector
 * (0, ..., 0, gn). */
static Expr* vp_matrix(Expr*** dv, size_t n, long repl, const Expr* gn) {
    Expr** rows = malloc(n * sizeof(Expr*));
    for (size_t k = 0; k < n; k++) {
        Expr** cols = malloc(n * sizeof(Expr*));
        for (size_t j = 0; j < n; j++) {
            if (repl >= 0 && j == (size_t)repl)
                cols[j] = (k == n - 1) ? expr_copy((Expr*)gn) : expr_new_integer(0);
            else
                cols[j] = expr_copy(dv[k][j]);
        }
        rows[k] = expr_new_function(expr_new_symbol(SYM_List), cols, n);
        free(cols);
    }
    Expr* m = expr_new_function(expr_new_symbol(SYM_List), rows, n);
    free(rows);
    return m;
}

/* True if e carries an arbitrary/undefined function of the variable -- a head
 * symbol with no builtin and no DownValues (e.g. the forcing f in y'' = f(x)),
 * or an inert Derivative.  Indefinite integration of such a term never closes
 * and can hang, so variation of parameters routes it to the definite
 * convolution instead.  (Mirrors lie_has_undefined_function in dsolve_lie.c.) */
bool ds_has_undefined_function(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        if (h->data.symbol.name == SYM_Derivative) return true;
        SymbolDef* d = symtab_lookup(h->data.symbol.name);
        /* No SymbolDef at all, or one with neither a builtin nor DownValues, is
         * an arbitrary/undefined function head (e.g. the forcing f in f[x]). */
        if (!d || (!d->builtin_func && !d->down_values)) return true;
    } else if (ds_has_undefined_function(h)) {
        return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ds_has_undefined_function(e->data.function.args[i])) return true;
    return false;
}

/* 1/e (fresh Power[e, -1]). */
static Expr* vp_recip(Expr* e) {
    return expr_new_function(expr_new_symbol(SYM_Power),
               (Expr*[]){ e, expr_new_integer(-1) }, 2);
}

/* Definite-convolution (Green's-function) fallback for variation of parameters,
 * used when the indefinite integral does not close in elementary form -- the
 * arbitrary forcing g = f(x) and impulse g = DiracDelta(x - a) cases.  Returns
 *   x_p(x) = Integrate[ Sum_i basis_i(x) cof_i(s) g(s) / (a_n(s) W(s)), {s,0,x} ]
 * over a fresh dummy s, where cof_i = Det of the Wronskian matrix with column i
 * replaced by the unit vector e_n, and W = detW.  The causal kernel vanishes on
 * the diagonal (K(x,x) = 0), so x_p and its first n-1 derivatives are 0 at the
 * base point x0 = 0 -- a zero-IC IVP fits its constants to 0.  Evaluating the
 * integral sifts a DiracDelta g (integrate_dirac) and closes a polynomial g,
 * while an arbitrary f is left as the convolution integral.  NULL on a
 * degenerate build (vanishing Wronskian). */
static Expr* vp_definite_convolution(Expr*** dv, Expr** basis, size_t n,
                                     const Expr* g, const Expr* leadcoef,
                                     const Expr* detW, const char* xvar) {
    /* Fresh dummy s, distinct from xvar and any symbol occurring in g. */
    const char* sname = intern_symbol("DSolve`vpS");
    if (strcmp(sname, xvar) == 0 || ds_contains(g, sname)) {
        char buf[32];
        for (int k = 1; ; k++) {
            snprintf(buf, sizeof buf, "DSolve`vpS%d", k);
            sname = intern_symbol(buf);
            if (strcmp(sname, xvar) != 0 && !ds_contains(g, sname)) break;
        }
    }
    Expr* W_s = ds_subst(expr_copy((Expr*)detW),
                         expr_new_symbol(xvar), expr_new_symbol(sname));
    /* Simplify the Wronskian so trig identities collapse (e.g. the
     * -3 Cos[3s]^2 - 3 Sin[3s]^2 denominator to -3); without this a resonant
     * cos/sin forcing convolution cannot be integrated in closed form. */
    W_s = ds_simplify(W_s);
    if (ds_is_zero(W_s)) { expr_free(W_s); return NULL; }
    Expr* g_s  = ds_subst(expr_copy((Expr*)g),
                          expr_new_symbol(xvar), expr_new_symbol(sname));
    Expr* lc_s = ds_subst(expr_copy((Expr*)leadcoef),
                          expr_new_symbol(xvar), expr_new_symbol(sname));
    Expr* one = expr_new_integer(1);
    Expr** summ = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        Expr* Wi   = vp_matrix(dv, n, (long)i, one);
        Expr* cof  = eval_and_free(ds_call1("Det", Wi));
        Expr* cofs = ds_subst(cof, expr_new_symbol(xvar), expr_new_symbol(sname));
        Expr* factors[5] = {
            expr_copy(basis[i]), cofs, expr_copy(g_s),
            vp_recip(expr_copy(lc_s)), vp_recip(expr_copy(W_s))
        };
        summ[i] = eval_and_free(
            expr_new_function(expr_new_symbol(SYM_Times), factors, 5));
    }
    expr_free(one);
    Expr* integrand = eval_and_free(
        expr_new_function(expr_new_symbol(SYM_Plus), summ, n));
    free(summ);
    /* TrigReduce linearises trig products/powers (Cos[a s]^2, Cos Sin, ...) so a
     * resonant cos/sin forcing closes to the clean t Sin[w t] form instead of an
     * unwieldy half-angle antiderivative; Expand then re-distributes, since
     * TrigReduce factors an exponential kernel as E^(a(t-s))(...), a form the
     * integrator churns on -- the expanded sum of E^(a(t-s)) f(s) terms
     * integrates termwise and fast. */
    integrand = eval_and_free(ds_call1("Expand",
                    eval_and_free(ds_call1("TrigReduce", integrand))));
    Expr* spec = expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){ expr_new_symbol(sname), expr_new_integer(0),
                   expr_new_symbol(xvar) }, 3);
    /* For an arbitrary forcing the integral stays an (unevaluated) convolution;
     * integrate_definite skips the improper/parametric methods on an undefined-
     * function integrand, so this closes fast for polynomial / trig / impulse
     * forcing and returns quickly unevaluated otherwise. */
    Expr* xp = eval_and_free(ds_call2(SYM_Integrate, integrand, spec));
    expr_free(W_s); expr_free(g_s); expr_free(lc_s);
    return xp;
}

Expr* dsolve_variation_of_parameters(Expr** basis, size_t n, const Expr* g,
                                     const Expr* leadcoef, const char* xvar) {
    Expr*** dv = malloc(n * sizeof(Expr**));
    for (size_t k = 0; k < n; k++) {
        dv[k] = malloc(n * sizeof(Expr*));
        for (size_t j = 0; j < n; j++)
            dv[k][j] = (k == 0) ? expr_copy(basis[j]) : ds_d(expr_copy(dv[k - 1][j]), expr_new_symbol(xvar));
    }
    Expr* W = vp_matrix(dv, n, -1, NULL);
    Expr* detW = eval_and_free(ds_call1("Det", W));
    Expr* yp = NULL;
    if (!ds_is_zero(detW)) {
        /* Collapse a trig/rational Wronskian (Cos^2+Sin^2 -> 1, ...) so the
         * per-term integrals close in elementary form and the recovered constant
         * factors are clean. */
        detW = ds_simplify(detW);
        Expr* gn = eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)g),
                        expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ expr_copy((Expr*)leadcoef), expr_new_integer(-1) }, 2)));
        /* DiracDelta (impulse) forcing genuinely needs the causal Green's-function
         * convolution: the impulse must be sifted at its base point and collapses
         * to 0 under indefinite integration.  EVERY other forcing -- elementary or
         * an arbitrary f(x) -- is handled by the per-term INDEFINITE Wronskian
         * integral, keeping any non-closing term as an inert Integrate.  That is
         * correct (D[Integrate]=integrand, and the inert terms carry the
         * coefficient L[basis_i]=0 in the residual, so verification is
         * Integrate-free), matches Mathematica's own integral-form answer for a
         * non-elementary / arbitrary forcing, and -- crucially -- never enters the
         * symbolic-limit definite integral, whose parametric DiffUnderInt
         * escalation blows up (Exp->Cosh/Sinh) on any forcing that is not an
         * undefined function (the §2.2.14 1337/1341/1350/1354 hang). */
        if (ds_contains(g, intern_symbol("DiracDelta"))) {
            expr_free(gn);
            yp = vp_definite_convolution(dv, basis, n, g, leadcoef, detW, xvar);
        } else {
            Expr** ut = malloc(n * sizeof(Expr*));
            bool any_inert = false;
            g_integrate_quiet++;   /* a non-closing Wronskian integral is kept inert */
            for (size_t i = 0; i < n; i++) {
                Expr* Wi = vp_matrix(dv, n, (long)i, gn);
                Expr* detWi = eval_and_free(ds_call1("Det", Wi));
                Expr* uip = eval_and_free(ds_call2(SYM_Times, detWi,
                                expr_new_function(expr_new_symbol(SYM_Power),
                                    (Expr*[]){ expr_copy(detW), expr_new_integer(-1) }, 2)));
                Expr* ui = ds_integrate(uip, expr_new_symbol(xvar));
                if (ds_has_head(ui, SYM_Integrate)) any_inert = true;
                ut[i] = eval_and_free(ds_call2(SYM_Times, expr_copy(basis[i]), ui));
            }
            g_integrate_quiet--;
            expr_free(gn);
            yp = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), ut, n));
            free(ut);
            /* Simplify only a fully-closed elementary answer; an inert-Integrate
             * body is left as-is (Simplify cannot help and could churn). */
            if (!any_inert) yp = ds_simplify(yp);
        }
    }
    expr_free(detW);
    for (size_t k = 0; k < n; k++) { for (size_t j = 0; j < n; j++) expr_free(dv[k][j]); free(dv[k]); }
    free(dv);
    return yp;
}

Expr* dsolve_run(DSolveProblem* P, DSolveTryFn fn) {
    size_t nb = 0;
    Expr** bodies = fn(P, &nb);
    if (!bodies) return NULL;
    if (nb == 0) { free(bodies); return NULL; }

    /* Fit every verified branch, recording HOW each fit turned out, then decide
     * per branch WITH the sibling context (dsolve_fit_constants cannot): an
     * Undefined fit is always dropped; a scalar unsatisfiable-condition branch
     * (FIT_EMPTY) is dropped only when a sibling actually fits (the wrong +/- or
     * Root branch of an IVP), and kept when none do (a lone basis singularity). */
    Expr** finals = malloc(nb * sizeof(Expr*));
    int*   fstate = malloc(nb * sizeof(int));
    size_t nf = 0, n_verified = 0, n_nosol = 0, n_ok = 0;
    for (size_t b = 0; b < nb; b++) {
        if (!bodies[b]) continue;
        if (!dsolve_verify_body(P, bodies[b])) { expr_free(bodies[b]); continue; }
        n_verified++;
        bool nosol = false; int st = FIT_OK;
        Expr* fitted = dsolve_fit_constants(P, bodies[b], &nosol, &st);
        expr_free(bodies[b]);
        if (nosol) { n_nosol++; if (fitted) expr_free(fitted); continue; }
        if (!fitted) continue;
        if (st == FIT_UNDEF) { expr_free(fitted); continue; }   /* never a solution */
        /* Drop a first-order explicit branch the corpus harness itself would score BAD
         * (a wrong root, or a form valid only on a sub-interval the grid overshoots):
         * declining lets a later cascade method's verifiable (often implicit) form win,
         * guaranteeing no explicit branch is shipped that back-substitutes nonzero on the
         * harness grid.  Prelude-matching + first-order-gated, so it can only remove
         * would-be-FAIL branches, never one the harness would PASS. */
        if (!ds_branch_corpus_verifiable(P, fitted)) { expr_free(fitted); continue; }
        fstate[nf] = st;
        finals[nf++] = fitted;
        if (st == FIT_OK) n_ok++;
    }
    free(bodies);
    /* Conditions present but NO branch achieved a real fit, and at least one branch
     * is FIT_UNDECIDED (Solve bubbled back, the constant unfitted): this method did
     * not solve the IVP -- its "solution" is the unfitted general form.  Decline so
     * the cascade continues to a method that CAN fit the condition (e.g. the exact /
     * homogeneous overlap 2.2.13-1205/1231, where Homogeneous's transcendental
     * log-form leaves C[1] and Exact's polynomial first integral fits it).
     * Extended to a SECOND-ORDER scalar IVP too: a special-function general
     * solution whose basis is SINGULAR at the IC point cannot be fitted there
     * (C[1]√x BesselJ[1/4,x²/2]+C[2]√x BesselY[1/4,x²/2] at x=0, §2.2.14-1381 --
     * BesselY[1/4,0] is infinite), so the fit bubbles FIT_UNDECIDED; declining lets
     * the cascade reach the Frobenius ordinary-point series, which fits the ICs
     * cleanly (now that dsolve_fit_constants Normal-izes a SeriesData body).  An
     * under-determined BVP is unaffected: its Solve SUCCEEDS (FIT_OK => n_ok>0), so
     * the n_ok==0 guard never fires for it. */
    if (P->ncond > 0 && n_ok == 0 && P->nfun == 1 && P->max_order[0] <= 2) {
        bool any_undecided = false;
        for (size_t b = 0; b < nf; b++) if (fstate[b] == FIT_UNDECIDED) any_undecided = true;
        if (any_undecided) {
            for (size_t b = 0; b < nf; b++) expr_free(finals[b]);
            free(finals); free(fstate);
            return NULL;
        }
    }
    if (n_ok > 0) {                 /* a sibling fit: drop the unsatisfiable branches */
        size_t keep = 0;
        for (size_t b = 0; b < nf; b++) {
            if (fstate[b] == FIT_EMPTY) { expr_free(finals[b]); }
            else finals[keep++] = finals[b];
        }
        nf = keep;
    }
    free(fstate);
    if (nf == 0) {
        free(finals);
        /* A verified general solution whose (boundary) conditions Solve proves
         * inconsistent is a well-posed problem with NO solution: return {} — the
         * concrete empty solution list, distinct from NULL ("decline / unevaluated"). */
        if (n_verified > 0 && n_nosol > 0)
            return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        return NULL;
    }

    Expr* result = dsolve_assemble(P, finals, nf);
    for (size_t b = 0; b < nf; b++) expr_free(finals[b]);
    free(finals);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Implicit (first-integral) solutions of a first-order ODE           *
 *                                                                     *
 *  When the general integral G(x, y) == C[1] cannot be solved         *
 *  explicitly for y (e.g. the homogeneous log-spiral                  *
 *  y' == (x+y)/(x-y), whose integral is ArcTan[y/x] - Log[x^2+y^2]/2),*
 *  DSolve returns that equation as the solution branch.  A try-fn on  *
 *  this path returns the left-hand side G (the relation is G == C[1]).*
 *  Verification differentiates the relation implicitly — from         *
 *  d/dx[G == C] = 0 the ODE forces y' == -G_x/G_y — and substitutes   *
 *  that y' into the residual; an IVP y[x0]==y0 fits the constant as   *
 *  G(x0, y0).  Output is {{ G(x, y[x]) == C[1] }}, not a y[x] -> rule.*/
static bool dsolve_verify_implicit(const DSolveProblem* P, const Expr* G) {
    if (P->nfun != 1) return false;
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    const char* Yn = intern_symbol("DSolve`impY");
    /* partials via G(x, Y):  y' = -D[G,x] / D[G,Y] */
    Expr* GY = ds_subst(expr_copy((Expr*)G), ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));
    Expr* Gx = ds_d(expr_copy(GY), expr_new_symbol(xvar));
    Expr* Gy = ds_d(GY, expr_new_symbol(Yn));                  /* consumes GY */
    Expr* yp = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                   ds_call2(SYM_Times, Gx,
                       expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ Gy, expr_new_integer(-1) }, 2))));
    yp = ds_subst(yp, expr_new_symbol(Yn), ds_make_funcapp(yname, 0, xvar));
    bool ok = true;
    for (size_t e = 0; e < P->neq && ok; e++) {
        Expr* sub = ds_subst(expr_copy(P->eq_residuals[e]),
                             ds_make_funcapp(yname, 1, xvar), expr_copy(yp));
        /* An implicit residual is F_x + N*(-F_x/F_y), which telescopes to 0 over a
         * common denominator.  The numeric zero-test can false-NEGATIVE on its
         * uncancelled Csc/Cot poles (the mu = Sin y exact family, e.g.
         * E^x + (E^x Cot y + 2 y Csc y) y' == 0) and so REJECT a correct branch.
         * Combine over a common denominator first: value-preserving (cannot turn a
         * genuinely nonzero residual into 0) and it renders the telescope a
         * syntactic 0 that the decider settles cleanly. */
        sub = eval_and_free(ds_call1("Together", sub));
        if (zero_test_decide(sub) == ZERO_TEST_FALSE) ok = false;
        expr_free(sub);
    }
    expr_free(yp);
    return ok;
}

/* The constant for the relation G == C: fitted to a first-order initial
 * condition y[x0]==y0 when one is present (C = G(x0, y0)), else C[1]. */
static Expr* dsolve_implicit_rhs(const DSolveProblem* P, const Expr* G) {
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    for (size_t c = 0; c < P->ncond; c++) {
        if (P->conds[c].fi != 0 || P->conds[c].order != 0) continue;
        Expr* v = expr_copy((Expr*)G);
        v = ds_subst(v, ds_make_funcapp(yname, 0, xvar), expr_copy(P->conds[c].value));
        v = ds_subst(v, expr_new_symbol(xvar), expr_copy(P->conds[c].point));
        return eval_and_free(v);
    }
    return ds_const(1);
}

Expr* dsolve_run_implicit(DSolveProblem* P, DSolveTryFn fn) {
    size_t nb = 0;
    Expr** gs = fn(P, &nb);
    if (!gs) return NULL;
    if (nb == 0) { free(gs); return NULL; }
    Expr** branches = malloc(nb * sizeof(Expr*));
    size_t nf = 0;
    for (size_t b = 0; b < nb; b++) {
        if (!gs[b]) continue;
        if (!dsolve_verify_implicit(P, gs[b])) { expr_free(gs[b]); continue; }
        Expr* rhs = dsolve_implicit_rhs(P, gs[b]);
        Expr* rhs2 = ds_rename_param(rhs, P->param_head);      /* renames the C[1] */
        expr_free(rhs);
        Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                       (Expr*[]){ gs[b], rhs2 }, 2);           /* consumes gs[b], rhs2 */
        branches[nf++] = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ eq }, 1);
    }
    free(gs);
    if (nf == 0) { free(branches); return NULL; }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), branches, nf);
    free(branches);
    return out;
}

Expr* dsolve_method_builtin_implicit(Expr* res, DSolveTryFn fn) {
    DSolveProblem P;
    if (!dsolve_parse(res, &P)) return NULL;
    if (P.is_pde || P.nfun != 1) { dsolve_problem_free(&P); return NULL; }
    Expr* r = dsolve_run_implicit(&P, fn);
    dsolve_problem_free(&P);
    return r;
}

/* ------------------------------------------------------------------ *
 *  Parametric solutions { x == X(t), y == Y(t) } (Lagrange/d'Alembert)*
 * ------------------------------------------------------------------ */

/* Verify a parametric branch: substitute x -> X, y[x] -> Y, and
 * y'[x] -> D[Y,t]/D[X,t] into each residual and require it not to be a decidably
 * non-zero value.  Elementary X,Y reduce via zero_test; for a transcendental pair
 * (where zero_test cannot prove D[Y,t]/D[X,t] == t) a PossibleZeroQ numeric
 * sampling of the parameter is the fallback before rejecting — the same permissive
 * policy as the implicit path. */
static bool dsolve_verify_parametric(const DSolveProblem* P, const Expr* X,
                                     const Expr* Y, const char* tname) {
    if (P->nfun != 1) return false;
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    Expr* Yp = ds_d(expr_copy((Expr*)Y), expr_new_symbol(tname));   /* dY/dt */
    Expr* Xp = ds_d(expr_copy((Expr*)X), expr_new_symbol(tname));   /* dX/dt */
    Expr* yprime = eval_and_free(ds_call2(SYM_Times, Yp,
                       expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ Xp, expr_new_integer(-1) }, 2)));  /* dY/dX */
    bool ok = true;
    for (size_t e = 0; e < P->neq && ok; e++) {
        Expr* sub = expr_copy(P->eq_residuals[e]);
        sub = ds_subst(sub, ds_make_funcapp(yname, 1, xvar), expr_copy(yprime));
        sub = ds_subst(sub, ds_make_funcapp(yname, 0, xvar), expr_copy((Expr*)Y));
        sub = ds_subst(sub, expr_new_symbol(xvar), expr_copy((Expr*)X));
        if (zero_test_decide(sub) == ZERO_TEST_FALSE) {
            /* transcendental fallback: numeric sampling of the parameter */
            Expr* pz = eval_and_free(ds_call1("PossibleZeroQ", expr_copy(sub)));
            if (!(pz->type == EXPR_SYMBOL && pz->data.symbol.name == SYM_True)) ok = false;
            expr_free(pz);
        }
        expr_free(sub);
    }
    expr_free(yprime);
    return ok;
}

/* Assemble one parametric branch { x -> Function[{t}, X], y -> Function[{t}, Y] },
 * renaming C[k] to the GeneratedParameters head.  X, Y are consumed. */
static Expr* dsolve_assemble_parametric(const DSolveProblem* P, Expr* X, Expr* Y,
                                        const char* tname) {
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    Expr* Xr = ds_rename_param(X, P->param_head);
    Expr* Yr = ds_rename_param(Y, P->param_head);
    expr_free(X); expr_free(Y);
    Expr* fx = expr_new_function(expr_new_symbol(SYM_Function), (Expr*[]){
        expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ expr_new_symbol(tname) }, 1), Xr }, 2);
    Expr* fy = expr_new_function(expr_new_symbol(SYM_Function), (Expr*[]){
        expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ expr_new_symbol(tname) }, 1), Yr }, 2);
    Expr* xrule = expr_new_function(expr_new_symbol(SYM_Rule),
                      (Expr*[]){ expr_new_symbol(xvar), fx }, 2);
    Expr* yrule = expr_new_function(expr_new_symbol(SYM_Rule),
                      (Expr*[]){ expr_new_symbol(yname), fy }, 2);
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ xrule, yrule }, 2);
}

/* Assemble one explicit scalar branch {y -> Function[{x}, body]} (or the applied
 * y[x] -> body form), renaming C[k].  body is borrowed. */
static Expr* dsolve_assemble_scalar_branch(const DSolveProblem* P, const Expr* body) {
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    Expr* b2 = ds_rename_param(body, P->param_head);
    Expr* lhs; Expr* rhs;
    if (P->applied) {
        lhs = expr_new_function(expr_new_symbol(yname), (Expr*[]){ expr_new_symbol(xvar) }, 1);
        rhs = b2;
    } else {
        lhs = expr_new_symbol(yname);
        rhs = expr_new_function(expr_new_symbol(SYM_Function), (Expr*[]){
            expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ expr_new_symbol(xvar) }, 1), b2 }, 2);
    }
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule), (Expr*[]){ lhs, rhs }, 2);
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ rule }, 1);
}

Expr* dsolve_run_parametric(DSolveProblem* P, DSolveTryFn fn) {
    if (P->ncond > 0) return NULL;             /* parametric IVP-fitting is future */
    size_t nb = 0;
    Expr** pairs = fn(P, &nb);
    if (!pairs) return NULL;
    if (nb == 0) { free(pairs); return NULL; }
    const char* wrapP = intern_symbol("DSolve`Param");
    const char* wrapE = intern_symbol("DSolve`Explicit");
    Expr** branches = malloc(nb * sizeof(Expr*));
    size_t nf = 0;
    for (size_t b = 0; b < nb; b++) {
        Expr* pr = pairs[b];
        if (!pr) continue;
        if (head_is(pr, wrapP) && pr->data.function.arg_count == 3
            && pr->data.function.args[2]->type == EXPR_SYMBOL) {
            Expr* X = expr_copy(pr->data.function.args[0]);
            Expr* Y = expr_copy(pr->data.function.args[1]);
            const char* tname = pr->data.function.args[2]->data.symbol.name;
            if (dsolve_verify_parametric(P, X, Y, tname))
                branches[nf++] = dsolve_assemble_parametric(P, X, Y, tname);  /* consumes X, Y */
            else { expr_free(X); expr_free(Y); }
        } else if (head_is(pr, wrapE) && pr->data.function.arg_count == 1) {
            /* explicit scalar branch (singular line): verify as y(x), assemble */
            const Expr* body = pr->data.function.args[0];
            if (dsolve_verify_body(P, body))
                branches[nf++] = dsolve_assemble_scalar_branch(P, body);
        }
        expr_free(pr);
    }
    free(pairs);
    if (nf == 0) { free(branches); return NULL; }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), branches, nf);
    free(branches);
    return out;
}

Expr* dsolve_method_builtin_parametric(Expr* res, DSolveTryFn fn) {
    DSolveProblem P;
    if (!dsolve_parse(res, &P)) return NULL;
    if (P.is_pde || P.nfun != 1) { dsolve_problem_free(&P); return NULL; }
    Expr* r = dsolve_run_parametric(&P, fn);
    dsolve_problem_free(&P);
    return r;
}

/* ------------------------------------------------------------------ *
 *  Systems (nfun > 1)                                                 *
 * ------------------------------------------------------------------ */

/* Renumber C[1..m] in `body` to C[*offset+1..*offset+m]; body consumed. */
Expr* dsolve_renumber_constants(Expr* body, int m, int* offset) {
    if (m <= 0) return body;
    Expr** rules = malloc((size_t)m * sizeof(Expr*));
    for (int j = 1; j <= m; j++)
        rules[j - 1] = expr_new_function(expr_new_symbol(SYM_Rule),
                           (Expr*[]){ ds_const(j), ds_const(*offset + j) }, 2);
    Expr* rl = expr_new_function(expr_new_symbol(SYM_List), rules, (size_t)m);
    free(rules);
    *offset += m;
    return eval_and_free(internal_replace_all((Expr*[]){ body, rl }, 2));
}

/* Extract the body from a scalar DSolve result {{fname -> Function[{x}, body]}}. */
Expr* dsolve_extract_system_body(Expr* r, const char* fname) {
    if (!head_is(r, SYM_List) || r->data.function.arg_count < 1) return NULL;
    Expr* inner = r->data.function.args[0];
    if (!head_is(inner, SYM_List)) return NULL;
    for (size_t k = 0; k < inner->data.function.arg_count; k++) {
        Expr* rule = inner->data.function.args[k];
        if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
            Expr* lhs = rule->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == fname) {
                Expr* rhs = rule->data.function.args[1];
                if (head_is(rhs, SYM_Function) && rhs->data.function.arg_count == 2)
                    return expr_copy(rhs->data.function.args[1]);
                return expr_copy(rhs);
            }
        }
    }
    return NULL;
}

Expr* dsolve_assemble_system(const DSolveProblem* P, Expr** bodies) {
    const char* xvar = P->ind_names[0];
    Expr** rules = malloc(P->nfun * sizeof(Expr*));
    for (size_t i = 0; i < P->nfun; i++) {
        Expr* body = ds_rename_param(bodies[i], P->param_head);
        Expr* lhs; Expr* rhs;
        if (P->applied) {
            lhs = expr_new_function(expr_new_symbol(P->fun_names[i]),
                                    (Expr*[]){ expr_new_symbol(xvar) }, 1);
            rhs = body;
        } else {
            lhs = expr_new_symbol(P->fun_names[i]);
            Expr* plist = expr_new_function(expr_new_symbol(SYM_List),
                                            (Expr*[]){ expr_new_symbol(xvar) }, 1);
            rhs = expr_new_function(expr_new_symbol(SYM_Function), (Expr*[]){ plist, body }, 2);
        }
        rules[i] = expr_new_function(expr_new_symbol(SYM_Rule), (Expr*[]){ lhs, rhs }, 2);
    }
    Expr* inner = expr_new_function(expr_new_symbol(SYM_List), rules, P->nfun);
    free(rules);
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ inner }, 1);
}

void dsolve_fit_system(const DSolveProblem* P, Expr** bodies, bool* no_solution) {
    if (no_solution) *no_solution = false;
    if (P->ncond == 0) return;
    const char* xvar = P->ind_names[0];
    Expr** params = NULL; size_t npar = 0;
    for (size_t i = 0; i < P->nfun; i++) ds_collect_consts(bodies[i], &params, &npar);
    if (npar == 0) return;

    Expr** eqs = malloc(P->ncond * sizeof(Expr*));
    size_t neq = 0;
    for (size_t c = 0; c < P->ncond; c++) {
        size_t fi = P->conds[c].fi;
        if (fi >= P->nfun) continue;
        Expr* bexpr = expr_copy(bodies[fi]);
        for (int d = 0; d < P->conds[c].order; d++) bexpr = ds_d(bexpr, expr_new_symbol(xvar));
        bexpr = ds_subst(bexpr, expr_new_symbol(xvar), expr_copy(P->conds[c].point));
        eqs[neq++] = expr_new_function(expr_new_symbol(SYM_Equal),
                        (Expr*[]){ bexpr, eval_and_free(expr_copy(P->conds[c].value)) }, 2);
    }
    Expr* eqlist = expr_new_function(expr_new_symbol(SYM_List), eqs, neq); free(eqs);
    Expr* varlist = expr_new_function(expr_new_symbol(SYM_List), params, npar); free(params);
    Expr* solres = ds_solve(eqlist, varlist);
    if (solres && head_is(solres, SYM_List)) {
        if (solres->data.function.arg_count == 0) {
            if (no_solution) *no_solution = true;   /* over-determined: no solution */
        } else {
            Expr* branch = solres->data.function.args[0];
            if (head_is(branch, SYM_List)) {
                for (size_t i = 0; i < P->nfun; i++) {
                    Expr* fitted = eval_and_free(internal_replace_all(
                        (Expr*[]){ expr_copy(bodies[i]), expr_copy(branch) }, 2));
                    expr_free(bodies[i]); bodies[i] = fitted;
                }
            }
        }
    }
    if (solres) expr_free(solres);
}

bool dsolve_verify_system(const DSolveProblem* P, Expr** bodies) {
    const char* xvar = P->ind_names[0];
    Expr** rules = malloc(P->nfun * sizeof(Expr*));
    for (size_t j = 0; j < P->nfun; j++) {
        Expr* plist = expr_new_function(expr_new_symbol(SYM_List),
                                        (Expr*[]){ expr_new_symbol(xvar) }, 1);
        Expr* fn = expr_new_function(expr_new_symbol(SYM_Function),
                                     (Expr*[]){ plist, expr_copy(bodies[j]) }, 2);
        rules[j] = expr_new_function(expr_new_symbol(SYM_Rule),
                                     (Expr*[]){ expr_new_symbol(P->fun_names[j]), fn }, 2);
    }
    Expr* rulelist = expr_new_function(expr_new_symbol(SYM_List), rules, P->nfun);
    free(rules);
    bool ok = true;
    for (size_t e = 0; e < P->neq && ok; e++) {
        Expr* sub = eval_and_free(internal_replace_all(
            (Expr*[]){ expr_copy(P->eq_residuals[e]), expr_copy(rulelist) }, 2));
        if (zero_test_decide(sub) == ZERO_TEST_FALSE) ok = false;
        expr_free(sub);
    }
    expr_free(rulelist);
    return ok;
}

/* ---- PDE (single function of nind variables) ---- */
static Expr* pde_varlist(const DSolveProblem* P) {
    Expr** vs = malloc(P->nind * sizeof(Expr*));
    for (size_t i = 0; i < P->nind; i++) vs[i] = expr_new_symbol(P->ind_names[i]);
    Expr* l = expr_new_function(expr_new_symbol(SYM_List), vs, P->nind);
    free(vs);
    return l;
}

/* Derivative[o1,o2][u][v1,v2] (a PDE derivative term of any order). */
static Expr* pde_deriv_lit(const char* u, int o1, int o2, const char* v1, const char* v2) {
    Expr* d = expr_new_function(expr_new_symbol(SYM_Derivative),
                  (Expr*[]){ expr_new_integer(o1), expr_new_integer(o2) }, 2);
    Expr* du = expr_new_function(d, (Expr*[]){ expr_new_symbol(u) }, 1);
    return expr_new_function(du, (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2);
}

/* If `e` is the term Derivative[i,j][u][v1,v2] for function `u`, return i+j;
 * else -1.  Used to discover the PDE's order (max_order is not populated for
 * PDEs — ds_scan only recognises single-index Derivative[m][u][x]). */
static int pde_term_order(const Expr* e, const char* u) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return -1;
    const Expr* h = e->data.function.head;               /* Derivative[i,j][u] */
    if (h->type != EXPR_FUNCTION || h->data.function.arg_count != 1
        || h->data.function.args[0]->type != EXPR_SYMBOL
        || h->data.function.args[0]->data.symbol.name != u) return -1;
    const Expr* d = h->data.function.head;               /* Derivative[i,j]    */
    if (d->type != EXPR_FUNCTION || d->data.function.arg_count != 2
        || d->data.function.head->type != EXPR_SYMBOL
        || d->data.function.head->data.symbol.name != SYM_Derivative
        || d->data.function.args[0]->type != EXPR_INTEGER
        || d->data.function.args[1]->type != EXPR_INTEGER) return -1;
    return (int)(d->data.function.args[0]->data.integer
               + d->data.function.args[1]->data.integer);
}

/* Max i+j over all Derivative[i,j][u] terms in `e` (0 if none). */
static int pde_scan_order(const Expr* e, const char* u) {
    int best = pde_term_order(e, u);
    if (best < 0) best = 0;
    if (e && e->type == EXPR_FUNCTION) {
        int ho = pde_scan_order(e->data.function.head, u);
        if (ho > best) best = ho;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            int a = pde_scan_order(e->data.function.args[i], u);
            if (a > best) best = a;
        }
    }
    return best;
}

/* C[k][z_] :> rhs — a test-function rule replacing the arbitrary function C[k]. */
static Expr* pde_arb_rule(int k, const char* z, Expr* rhs) {
    Expr* blank = expr_new_function(expr_new_symbol("Blank"), NULL, 0);
    Expr* patt = expr_new_function(expr_new_symbol("Pattern"),
                     (Expr*[]){ expr_new_symbol(z), blank }, 2);
    Expr* lhs = expr_new_function(ds_const(k), (Expr*[]){ patt }, 1);
    return expr_new_function(expr_new_symbol(SYM_RuleDelayed), (Expr*[]){ lhs, rhs }, 2);
}

/* Verify a 2-variable PDE of any order.  Two obstacles: the evaluator does not
 * reduce Derivative[i,j][Function[...]][...] (so we substitute the derivative
 * TERMS with D[body, {v1,i}, {v2,j}], which does reduce), and zero_test cannot
 * sample an arbitrary function C[k][...] (so we first replace each with a
 * distinct concrete test function — Sin, Cos, Exp, #^2 — a correct general
 * solution stays a solution for any choice, and distinct functions keep an
 * error in one branch from cancelling against another).  The order is scanned
 * from the residual (max_order is 0 for PDEs), so this serves first- and
 * second-order methods alike. */
static bool dsolve_verify_pde(const DSolveProblem* P, const Expr* body) {
    const char* u = P->fun_names[0];
    const char* v1 = P->ind_names[0];
    const char* v2 = P->ind_names[1];
    const char* z = intern_symbol("DSolve`pdez");

    Expr* r1 = pde_arb_rule(1, z, ds_call1("Sin", expr_new_symbol(z)));
    Expr* r2 = pde_arb_rule(2, z, ds_call1("Cos", expr_new_symbol(z)));
    Expr* r3 = pde_arb_rule(3, z, ds_call1("Exp", expr_new_symbol(z)));
    Expr* r4 = pde_arb_rule(4, z, expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ expr_new_symbol(z), expr_new_integer(2) }, 2));
    Expr* rl = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ r1, r2, r3, r4 }, 4);
    Expr* bodyC = eval_and_free(internal_replace_all((Expr*[]){ expr_copy((Expr*)body), rl }, 2));

    int maxord = 0;
    for (size_t e = 0; e < P->neq; e++) {
        int o = pde_scan_order(P->eq_residuals[e], u);
        if (o > maxord) maxord = o;
    }
    if (maxord < 1) maxord = 1;

    bool ok = true;
    for (size_t e = 0; e < P->neq && ok; e++) {
        Expr* r = expr_copy(P->eq_residuals[e]);
        for (int s = maxord; s >= 1; s--)
            for (int i = s; i >= 0; i--) {
                int j = s - i;
                Expr* dk = expr_copy(bodyC);
                for (int t = 0; t < i; t++) dk = ds_d(dk, expr_new_symbol(v1));
                for (int t = 0; t < j; t++) dk = ds_d(dk, expr_new_symbol(v2));
                r = ds_subst(r, pde_deriv_lit(u, i, j, v1, v2), dk);
            }
        r = ds_subst(r, expr_new_function(expr_new_symbol(u),
                        (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2), expr_copy(bodyC));
        if (zero_test_decide(r) == ZERO_TEST_FALSE) ok = false;
        expr_free(r);
    }
    expr_free(bodyC);
    return ok;
}

static Expr* dsolve_assemble_pde(const DSolveProblem* P, Expr* body) {
    const char* uname = P->fun_names[0];
    Expr* b = ds_rename_param(body, P->param_head);
    Expr* lhs; Expr* rhs;
    if (P->applied) {
        Expr** vs = malloc(P->nind * sizeof(Expr*));
        for (size_t i = 0; i < P->nind; i++) vs[i] = expr_new_symbol(P->ind_names[i]);
        lhs = expr_new_function(expr_new_symbol(uname), vs, P->nind);
        free(vs);
        rhs = b;
    } else {
        lhs = expr_new_symbol(uname);
        rhs = expr_new_function(expr_new_symbol(SYM_Function),
                  (Expr*[]){ pde_varlist(P), b }, 2);
    }
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule), (Expr*[]){ lhs, rhs }, 2);
    Expr* inner = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ rule }, 1);
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ inner }, 1);
}

Expr* dsolve_run_pde(DSolveProblem* P, DSolveSysFn fn) {
    Expr** bodies = fn(P);
    if (!bodies) return NULL;
    Expr* result = NULL;
    if (bodies[0] && dsolve_verify_pde(P, bodies[0]))
        result = dsolve_assemble_pde(P, bodies[0]);
    if (bodies[0]) expr_free(bodies[0]);
    free(bodies);
    return result;
}

/* ---- implicit (first-integral) PDE solutions --------------------------------
 * A first-order quasilinear PDE  P u_v1 + Q u_v2 == R  has, in general, no
 * explicit u -> body form: its general solution is the implicit relation
 * phi1(v1,v2,u) == C[1][phi2(v1,v2,u)] between two independent first integrals
 * of the characteristic system dv1/P = dv2/Q = du/R (Lagrange).  A PDE method on
 * this path returns, in bodies[0], one of two wrapper heads carrying its result
 * in terms of the BARE u-symbol DSolve`pdeU (so the verifier can differentiate
 * treating u as an independent coordinate):
 *   DSolve`PDEImplicit[phi1, phi2]  — the implicit relation above
 *   DSolve`PDEExplicit[body]        — an explicit u == body(v1,v2) solution
 *                                     (semilinear case; routed to the ordinary
 *                                     explicit PDE verify/assemble).
 * The implicit branch is verified by the implicit-function rule: with C[1]
 * pinned to a concrete test function F, Psi(v1,v2,u) = phi1 - F(phi2) defines
 * u(v1,v2) with u_vi = -Psi_vi / Psi_u; substituting into the residual must not
 * be a decidable non-zero (a PossibleZeroQ sampling is the transcendental
 * fallback, as in dsolve_verify_parametric). */
/* Shared implicit-diff core: with Psi(v1,v2,u) == const defining u(v1,v2)
 * implicitly, the implicit-function rule gives u_vi = -Psi_vi/Psi_u; substitute
 * into the single PDE residual and return false only if it is a DECIDABLE
 * non-zero (a PossibleZeroQ sampling is the transcendental fallback, as in
 * dsolve_verify_parametric).  Psi is in the bare u-symbol DSolve`pdeU and is
 * CONSUMED. */
static bool pde_implicit_residual_ok(const DSolveProblem* P, Expr* Psi) {
    const char* u  = P->fun_names[0];
    const char* v1 = P->ind_names[0];
    const char* v2 = P->ind_names[1];
    const char* sU = intern_symbol("DSolve`pdeU");
    Expr* Px = ds_d(expr_copy(Psi), expr_new_symbol(v1));
    Expr* Py = ds_d(expr_copy(Psi), expr_new_symbol(v2));
    Expr* Pu = ds_d(Psi, expr_new_symbol(sU));               /* consumes Psi */
    Expr* Puinv = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ Pu, expr_new_integer(-1) }, 2);
    Expr* ux = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                   ds_call2(SYM_Times, Px, expr_copy(Puinv))));
    Expr* uy = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                   ds_call2(SYM_Times, Py, expr_copy(Puinv))));
    expr_free(Puinv);

    Expr* r = expr_copy(P->eq_residuals[0]);
    r = ds_subst(r, pde_deriv_lit(u, 1, 0, v1, v2), ux);
    r = ds_subst(r, pde_deriv_lit(u, 0, 1, v1, v2), uy);
    r = ds_subst(r, expr_new_function(expr_new_symbol(u),
                    (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2),
                 expr_new_symbol(sU));
    bool ok = true;
    if (zero_test_decide(r) == ZERO_TEST_FALSE) {
        Expr* pz = eval_and_free(ds_call1("PossibleZeroQ", expr_copy(r)));
        if (!(pz->type == EXPR_SYMBOL && pz->data.symbol.name == SYM_True)) ok = false;
        expr_free(pz);
    }
    expr_free(r);
    return ok;
}

static bool dsolve_verify_pde_implicit(const DSolveProblem* P,
                                       const Expr* phi1, const Expr* phi2) {
    if (P->neq < 1) return false;
    const char* z  = intern_symbol("DSolve`pdez");

    /* distinct concrete test functions for the arbitrary C[1]: Sin, Cos, #^2 */
    Expr* tests[3];
    tests[0] = ds_call1("Sin", expr_new_symbol(z));
    tests[1] = ds_call1("Cos", expr_new_symbol(z));
    tests[2] = expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ expr_new_symbol(z), expr_new_integer(2) }, 2);

    bool ok = true;
    for (int t = 0; t < 3 && ok; t++) {
        Expr* Fphi2 = ds_subst(expr_copy(tests[t]), expr_new_symbol(z),
                               expr_copy((Expr*)phi2));
        Expr* Psi = eval_and_free(ds_call2(SYM_Subtract,
                        expr_copy((Expr*)phi1), Fphi2));     /* Psi = phi1 - F(phi2) */
        ok = pde_implicit_residual_ok(P, Psi);               /* consumes Psi */
        expr_free(tests[t]);
    }
    return ok;
}

/* Verify a Charpit complete integral given as the implicit relation Psi == C[k]
 * (Psi in the bare u-symbol, with arbitrary CONSTANTS — no arbitrary function to
 * pin).  Just the implicit-diff core once. */
static bool dsolve_verify_pde_relation(const DSolveProblem* P, const Expr* Psi) {
    if (P->neq < 1) return false;
    return pde_implicit_residual_ok(P, expr_copy((Expr*)Psi));
}

/* Assemble {{ Psi(v1,v2,u[v1,v2]) == C[2] }} (Psi carries the first constant C[1]),
 * renaming C[k] to the GeneratedParameters head.  Psi borrowed. */
static Expr* dsolve_assemble_pde_relation(const DSolveProblem* P, const Expr* Psi) {
    const char* uname = P->fun_names[0];
    const char* v1 = P->ind_names[0];
    const char* v2 = P->ind_names[1];
    const char* sU = intern_symbol("DSolve`pdeU");
    Expr* uapp = expr_new_function(expr_new_symbol(uname),
                     (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2);
    Expr* psi_app = ds_subst(expr_copy((Expr*)Psi), expr_new_symbol(sU), uapp);
    Expr* rel = expr_new_function(expr_new_symbol(SYM_Equal),
                    (Expr*[]){ psi_app, ds_const(2) }, 2);
    Expr* rel2 = ds_rename_param(rel, P->param_head);
    expr_free(rel);
    Expr* inner = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ rel2 }, 1);
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ inner }, 1);
}

/* Assemble {{ phi1(v1,v2,u[v1,v2]) == C[1][phi2(...)] }} from the bare-u first
 * integrals, renaming C[1] to the GeneratedParameters head.  phi1/phi2 borrowed. */
static Expr* dsolve_assemble_pde_implicit(const DSolveProblem* P,
                                          const Expr* phi1, const Expr* phi2) {
    const char* uname = P->fun_names[0];
    const char* v1 = P->ind_names[0];
    const char* v2 = P->ind_names[1];
    const char* sU = intern_symbol("DSolve`pdeU");
    Expr* uapp = expr_new_function(expr_new_symbol(uname),
                     (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2);
    Expr* p1 = ds_subst(expr_copy((Expr*)phi1), expr_new_symbol(sU), expr_copy(uapp));
    Expr* p2 = ds_subst(expr_copy((Expr*)phi2), expr_new_symbol(sU), uapp); /* consumes uapp */
    Expr* arb = expr_new_function(ds_const(1), (Expr*[]){ p2 }, 1);
    Expr* rel = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ p1, arb }, 2);
    Expr* rel2 = ds_rename_param(rel, P->param_head);
    expr_free(rel);
    Expr* inner = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ rel2 }, 1);
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ inner }, 1);
}

/* One verified explicit branch as the inner {u -> Function[...]}: assemble via
 * dsolve_assemble_pde (which returns {{rule}}) and lift its single inner list. */
static Expr* pde_branch_inner(const DSolveProblem* P, Expr* body) {
    Expr* wrapped = dsolve_assemble_pde(P, body);        /* {{ rule }} */
    Expr* inner = expr_copy(wrapped->data.function.args[0]);  /* { rule } */
    expr_free(wrapped);
    return inner;
}

Expr* dsolve_run_pde_implicit(DSolveProblem* P, DSolveSysFn fn) {
    Expr** bodies = fn(P);
    if (!bodies) return NULL;
    Expr* body = bodies[0];
    Expr* result = NULL;
    const char* wImpl = intern_symbol("DSolve`PDEImplicit");
    const char* wExpl = intern_symbol("DSolve`PDEExplicit");
    const char* wBran = intern_symbol("DSolve`PDEBranches");
    const char* wRel  = intern_symbol("DSolve`PDERelation");
    if (body && head_is(body, wImpl) && body->data.function.arg_count == 2) {
        Expr* phi1 = body->data.function.args[0];
        Expr* phi2 = body->data.function.args[1];
        if (dsolve_verify_pde_implicit(P, phi1, phi2))
            result = dsolve_assemble_pde_implicit(P, phi1, phi2);
    } else if (body && head_is(body, wRel) && body->data.function.arg_count == 1) {
        /* Charpit implicit complete integral: relation Psi == C[2] (arbitrary
         * constants, no arbitrary function). */
        Expr* Psi = body->data.function.args[0];
        if (dsolve_verify_pde_relation(P, Psi))
            result = dsolve_assemble_pde_relation(P, Psi);
    } else if (body && head_is(body, wExpl) && body->data.function.arg_count == 1) {
        Expr* b = body->data.function.args[0];
        if (dsolve_verify_pde(P, b))
            result = dsolve_assemble_pde(P, b);
    } else if (body && head_is(body, wBran)) {
        /* multiple explicit branches (complete integral + singular envelopes);
         * verify each and collect the survivors into {{...},...}. */
        size_t n = body->data.function.arg_count;
        Expr** inners = malloc((n ? n : 1) * sizeof(Expr*));
        size_t nf = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* b = body->data.function.args[i];
            if (dsolve_verify_pde(P, b))
                inners[nf++] = pde_branch_inner(P, b);
        }
        if (nf > 0) result = expr_new_function(expr_new_symbol(SYM_List), inners, nf);
        free(inners);
    }
    if (body) expr_free(body);
    free(bodies);
    return result;
}

Expr* dsolve_method_builtin_pde_implicit(Expr* res, DSolveSysFn fn) {
    DSolveProblem P;
    if (!dsolve_parse(res, &P)) return NULL;
    if (!P.is_pde) { dsolve_problem_free(&P); return NULL; }
    Expr* r = dsolve_run_pde_implicit(&P, fn);
    dsolve_problem_free(&P);
    return r;
}

Expr* dsolve_run_system(DSolveProblem* P, DSolveSysFn fn) {
    Expr** bodies = fn(P);
    if (!bodies) return NULL;
    Expr* result = NULL;
    if (dsolve_verify_system(P, bodies)) {
        bool nosol = false;
        dsolve_fit_system(P, bodies, &nosol);
        result = nosol ? expr_new_function(expr_new_symbol(SYM_List), NULL, 0)
                       : dsolve_assemble_system(P, bodies);
    }
    for (size_t i = 0; i < P->nfun; i++) if (bodies[i]) expr_free(bodies[i]);
    free(bodies);
    return result;
}

Expr* dsolve_method_builtin(Expr* res, DSolveTryFn fn) {
    DSolveProblem P;
    if (!dsolve_parse(res, &P)) return NULL;
    if (P.is_pde) { dsolve_problem_free(&P); return NULL; }
    /* Count a pinned method as one recursion level so its own sub-DSolve calls
     * (e.g. OperatorFactor peeling a quotient) run at depth >= 2 — matching how
     * they nest under the Automatic cascade — and so the depth-gated coefficient
     * normalization does not fire inside them. */
    g_dsolve_depth++;
    Expr* r = dsolve_run(&P, fn);
    g_dsolve_depth--;
    dsolve_problem_free(&P);
    return r;
}

Expr* dsolve_method_builtin_system(Expr* res, DSolveSysFn fn) {
    DSolveProblem P;
    if (!dsolve_parse(res, &P)) return NULL;
    if (P.is_pde) { dsolve_problem_free(&P); return NULL; }  /* PDE is a different method */
    Expr* r = dsolve_run_system(&P, fn);
    dsolve_problem_free(&P);
    return r;
}

Expr* dsolve_method_builtin_pde(Expr* res, DSolveSysFn fn) {
    DSolveProblem P;
    if (!dsolve_parse(res, &P)) return NULL;
    if (!P.is_pde) { dsolve_problem_free(&P); return NULL; }  /* ODE is a different method */
    Expr* r = dsolve_run_pde(&P, fn);
    dsolve_problem_free(&P);
    return r;
}

Expr* dsolve_algebraic_residual(DSolveProblem* P, const char* Yname, const char* Pname) {
    if (P->neq < 1 || P->nfun < 1 || P->nind < 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];
    Expr* R = expr_copy(P->eq_residuals[0]);
    R = ds_subst(R, ds_make_funcapp(yname, 1, xvar), expr_new_symbol(Pname)); /* y'[x] -> p */
    R = ds_subst(R, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yname)); /* y[x]  -> Y */
    return R;
}

Expr* dsolve_linear_factor_solve(Expr* Pcoef, Expr* Qcoef, const char* xvar) {
    g_integrate_quiet++;   /* a non-elementary integrating-factor integral is kept
                            * unevaluated below, not surfaced as a speculative message */
    Expr* Pint = ds_integrate(Pcoef, expr_new_symbol(xvar));       /* consumes Pcoef */
    if (ds_has_head(Pint, SYM_Integrate)) { g_integrate_quiet--; expr_free(Pint); expr_free(Qcoef); return NULL; }
    Expr* mu = eval_and_free(ds_call1("Exp", Pint));
    /* Simplify the integrating factor before integrating mu*q.  ∫p can come back
     * in an un-collapsed form (e.g. ∫Tan[x] as a sum/difference of three Logs
     * instead of -Log[Cos[x]]); Exp of that is an ugly rational-trig product that
     * makes Integrate[mu q] spin (the y' + Tan[x] y == q hang).  mu is defined
     * only up to a nonzero constant, so any simplified equivalent is valid, and
     * the constant cancels in the final mu^-1 (∫mu q + C[1]).  Elementary clean
     * cases (mu = x, Exp[x^2/2], …) are unchanged.  PowerExpand after Simplify
     * collapses the Sqrt[Sec[x]^2] / (1+Tan^2)^(1/2) artifacts that Simplify
     * alone leaves (branch caution) down to Cos[x]; mu is defined up to a
     * nonzero constant/branch and the final solution is verified by
     * back-substitution, so choosing a branch here is sound. */
    mu = eval_and_free(ds_call1("PowerExpand", ds_simplify(mu)));
    Expr* integrand = eval_and_free(ds_call2(SYM_Times, expr_copy(mu), Qcoef)); /* consumes Qcoef */
    Expr* Qint = ds_integrate(integrand, expr_new_symbol(xvar));
    g_integrate_quiet--;
    /* When ∫mu q dx is NON-elementary, Qint is a still-unevaluated Integrate[...].
     * Do NOT decline: the integrating-factor solution y = mu^-1 (∫mu q dx + C[1])
     * is exact regardless, and the unevaluated-integral form is exactly what
     * Mathematica/Maple return (and upgrades to a closed form automatically once
     * Integrate learns the case, e.g. Erfi).  Declining here previously dropped
     * y' + x y == Exp[3 x] through the whole cascade and spun the evaluator to
     * $IterationLimit, and let y' + y == Q[x] be fabricated wrong downstream by
     * UndeterminedCoefficients.  We keep Qint as-is. */
    Expr* num = eval_and_free(ds_call2(SYM_Plus, Qint, ds_const(1)));
    Expr* body = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
        num,
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ mu, expr_new_integer(-1) }, 2)
    }, 2));
    return body;
}

/* Collapse the integer-family parameters of a stripped ConditionalExpression.
 * A multivalued inverse (ArcSin / ArcCos / Log ...) returns v + 2 Pi C[k] gated
 * by Element[C[k], Integers]; for an ODE GENERAL solution that discrete family
 * is already subsumed by the continuous integration constant, so the principal
 * branch (each constrained C[k] -> 0) is the intended solution.  `out` is
 * consumed; `cond` is scanned for every C[integer] it mentions (the integration
 * constant never appears in the Element[...] condition, so it is left intact).
 * Mirrors the IC-path collapse in dsolve_fit_constants (M21). */
static Expr* collapse_cond_families(Expr* out, const Expr* cond) {
    if (!out || !cond || cond->type != EXPR_FUNCTION) return out;
    /* Collapse ONLY a C[k] that the condition constrains as an integer index
     * (Element[C[k], Integers]).  A range condition such as -Pi/2 < x^3 + C[1]
     * <= Pi/2 (from a Tan/ArcTan inversion) legitimately mentions the continuous
     * integration constant C[1]; collapsing that would zero the constant and
     * break the IVP fit (regression on y'=3x^2(1+y^2), y(0)=1). */
    if (head_is(cond, intern_symbol("Element")) && cond->data.function.arg_count == 2) {
        const Expr* p   = cond->data.function.args[0];
        const Expr* dom = cond->data.function.args[1];
        if (dom->type == EXPR_SYMBOL && dom->data.symbol.name == intern_symbol("Integers")
            && head_is(p, intern_symbol("C")) && p->data.function.arg_count == 1
            && p->data.function.args[0]->type == EXPR_INTEGER)
            return ds_subst(out, expr_copy((Expr*)p), expr_new_integer(0));
        return out;
    }
    /* Recurse through boolean combinators (And / Or / ...) to reach every
     * Element[...] atom; never collapse a bare C[k] found outside one. */
    for (size_t i = 0; i < cond->data.function.arg_count; i++)
        out = collapse_cond_families(out, cond->data.function.args[i]);
    return out;
}

Expr** dsolve_extract_solutions(Expr* solres, const char* varname, size_t* n) {
    *n = 0;
    if (!solres || !head_is(solres, SYM_List)) return NULL;
    size_t m = solres->data.function.arg_count;
    Expr** out = malloc((m ? m : 1) * sizeof(Expr*));
    size_t c = 0;
    for (size_t i = 0; i < m; i++) {
        Expr* br = solres->data.function.args[i];
        if (!head_is(br, SYM_List)) continue;
        for (size_t j = 0; j < br->data.function.arg_count; j++) {
            Expr* rule = br->data.function.args[j];
            if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
                Expr* lhs = rule->data.function.args[0];
                if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == varname) {
                    Expr* val = rule->data.function.args[1];
                    /* strip ConditionalExpression[v, cond] -> v so the body
                     * verifies, collapsing any integer periodicity family the
                     * condition constrains (2 Pi C[k]) to its principal branch */
                    if (val->type == EXPR_FUNCTION && val->data.function.arg_count == 2
                        && head_is(val, intern_symbol("ConditionalExpression")))
                        out[c++] = collapse_cond_families(
                                       expr_copy(val->data.function.args[0]),
                                       val->data.function.args[1]);
                    else
                        out[c++] = expr_copy(val);
                }
            }
        }
    }
    *n = c;
    if (c == 0) { free(out); return NULL; }
    return out;
}

Expr** dsolve_extract_applied_bodies(Expr* r, const char* fname, size_t* n) {
    *n = 0;
    if (!r || !head_is(r, SYM_List)) return NULL;
    size_t m = r->data.function.arg_count;
    Expr** out = malloc((m ? m : 1) * sizeof(Expr*));
    size_t c = 0;
    for (size_t i = 0; i < m; i++) {
        Expr* inner = r->data.function.args[i];
        if (!head_is(inner, SYM_List)) continue;
        for (size_t k = 0; k < inner->data.function.arg_count; k++) {
            Expr* rule = inner->data.function.args[k];
            if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
                Expr* lhs = rule->data.function.args[0];
                if (lhs->type == EXPR_FUNCTION
                    && lhs->data.function.head->type == EXPR_SYMBOL
                    && lhs->data.function.head->data.symbol.name == fname)
                    out[c++] = expr_copy(rule->data.function.args[1]);
            }
        }
    }
    *n = c;
    if (c == 0) { free(out); return NULL; }
    return out;
}

/* ---- shared fundamental-set builder (LinearConstantCoefficients + UC) ---- */
static Expr* hb_xpow(const char* xvar, int j) {
    return expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ expr_new_symbol(xvar), expr_new_integer(j) }, 2);
}
static Expr* hb_exp_lin(Expr* coef, const char* xvar) {
    return eval_and_free(ds_call1("Exp", ds_call2(SYM_Times, coef, expr_new_symbol(xvar))));
}
static Expr* hb_basis_exp(const char* xvar, int j, const Expr* r) {
    Expr* e = hb_exp_lin(expr_copy((Expr*)r), xvar);
    if (j == 0) return e;
    return eval_and_free(ds_call2(SYM_Times, hb_xpow(xvar, j), e));
}
static Expr* hb_basis_trig(const char* xvar, int j, const Expr* a, const Expr* b, const char* trig) {
    Expr* e = hb_exp_lin(expr_copy((Expr*)a), xvar);
    Expr* t = ds_call1(trig, ds_call2(SYM_Times, expr_copy((Expr*)b), expr_new_symbol(xvar)));
    Expr* prod = eval_and_free(ds_call2(SYM_Times, e, t));
    if (j == 0) return prod;
    return eval_and_free(ds_call2(SYM_Times, hb_xpow(xvar, j), prod));
}

Expr** dsolve_homog_basis(const Expr* charpoly, const char* lam, const char* xvar,
                          int n, size_t* count) {
    *count = 0;
    DSolveRoots R;
    if (!dsolve_analyze_roots(charpoly, lam, n, &R)) return NULL;
    Expr** basis = NULL;
    if (R.total == n) {
        bool* used = calloc(R.ndist, sizeof(bool));
        basis = malloc((size_t)n * sizeof(Expr*));
        size_t bc = 0;
        for (size_t i = 0; i < R.ndist; i++) {
            if (used[i]) continue;
            if (R.isreal[i]) {
                for (int j = 0; j < R.mult[i]; j++) basis[bc++] = hb_basis_exp(xvar, j, R.roots[i]);
                used[i] = true;
            } else {
                Expr* conj = eval_and_free(ds_call1("Conjugate", expr_copy(R.roots[i])));
                long c = -1;
                for (size_t k = 0; k < R.ndist && c < 0; k++) {
                    if (k == i || used[k]) continue;
                    Expr* diff = eval_and_free(ds_call2(SYM_Subtract, expr_copy(R.roots[k]), expr_copy(conj)));
                    if (ds_is_zero(diff)) c = (long)k;
                    expr_free(diff);
                }
                expr_free(conj);
                if (c >= 0) {
                    used[i] = used[(size_t)c] = true;
                    Expr* a = eval_and_free(ds_call1("Re", expr_copy(R.roots[i])));
                    Expr* imv = expr_copy(R.im[i]);
                    /* Concretize Re/Im of a NUMERIC complex root before it enters the
                     * basis exponent: Re/Im do not auto-evaluate on a radical power
                     * (e.g. the cube root -(-1)^(1/3) of y'''==y leaves
                     * Re[-(-1)^(1/3)]/Im[-(-1)^(1/3)] in the exponent), which a later
                     * IVP constant-fit cannot solve to a number -> the IVP is scored
                     * unfitted (§2.2.4-312).  ComplexExpand computes them for a pure
                     * number (safe: no free symbol, so no Abs/Sign is introduced); a
                     * root carrying a symbolic parameter is left untouched (its general
                     * solution still back-substitutes). */
                    Expr* isnum = eval_and_free(ds_call1("NumericQ", expr_copy(R.roots[i])));
                    if (isnum && isnum->type == EXPR_SYMBOL && isnum->data.symbol.name == SYM_True) {
                        a   = eval_and_free(ds_call1("ComplexExpand", a));
                        imv = eval_and_free(ds_call1("ComplexExpand", imv));
                    }
                    expr_free(isnum);
                    for (int j = 0; j < R.mult[i]; j++) {
                        basis[bc++] = hb_basis_trig(xvar, j, a, imv, "Cos");
                        basis[bc++] = hb_basis_trig(xvar, j, a, imv, "Sin");
                    }
                    expr_free(a); expr_free(imv);
                } else {
                    for (int j = 0; j < R.mult[i]; j++) basis[bc++] = hb_basis_exp(xvar, j, R.roots[i]);
                    used[i] = true;
                }
            }
        }
        *count = bc;
        free(used);
    }
    dsolve_roots_free(&R);
    return basis;
}

Expr* dsolve_solve_top_derivative(DSolveProblem* P, int n) {
    if (P->neq < 1 || P->nfun < 1 || P->nind < 1 || n < 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];
    Expr* R = P->eq_residuals[0];
    const char* pn = intern_symbol("DSolve`dp");

    Expr* topLit = ds_make_funcapp(yname, n, xvar);
    Expr* Rp = ds_subst(expr_copy(R), topLit, expr_new_symbol(pn));   /* consumes topLit + new pn */
    Expr* a = ds_d(expr_copy(Rp), expr_new_symbol(pn));
    if (ds_contains(a, pn) || ds_is_zero(a)) { expr_free(a); expr_free(Rp); return NULL; }
    Expr* b = ds_subst(Rp, expr_new_symbol(pn), expr_new_integer(0));  /* consumes Rp */
    /* F = -b / a */
    Expr* F = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
        expr_new_integer(-1), b,
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ a, expr_new_integer(-1) }, 2)
    }, 3));
    return F;
}
