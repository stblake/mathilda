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

#include "../sym_names.h"
#include "../sym_intern.h"
#include "../eval.h"
#include "../common.h"
#include "../internal.h"
#include "../zero_test.h"
#include "../ndarray.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

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
static Expr* ds_rename_param(const Expr* e, const char* head) {
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
        ZeroTestResult zt = zero_test_decide(sub);
        expr_free(sub);
        if (zt == ZERO_TEST_FALSE) return false;
    }
    return true;
}

/* Fit generated constants to the initial/boundary conditions; returns a fresh
 * body (the general body copied when there is nothing to fit). */
static Expr* dsolve_fit_constants(const DSolveProblem* P, const Expr* body) {
    if (P->ncond == 0) return expr_copy((Expr*)body);
    Expr** params = NULL; size_t npar = 0;
    ds_collect_consts(body, &params, &npar);
    if (npar == 0) return expr_copy((Expr*)body);
    const char* xvar = P->ind_names[0];

    Expr** eqs = malloc(P->ncond * sizeof(Expr*));
    size_t neq = 0;
    for (size_t c = 0; c < P->ncond; c++) {
        if (P->conds[c].fi != 0) continue;               /* single-function M0 */
        Expr* bexpr = expr_copy((Expr*)body);
        for (int d = 0; d < P->conds[c].order; d++)
            bexpr = ds_d(bexpr, expr_new_symbol(xvar));
        bexpr = ds_subst(bexpr, expr_new_symbol(xvar), expr_copy(P->conds[c].point));
        eqs[neq++] = expr_new_function(expr_new_symbol(SYM_Equal),
                        (Expr*[]){ bexpr, eval_and_free(expr_copy(P->conds[c].value)) }, 2);
    }
    Expr* eqlist = expr_new_function(expr_new_symbol(SYM_List), eqs, neq);
    free(eqs);
    Expr* varlist = expr_new_function(expr_new_symbol(SYM_List), params, npar);
    free(params);

    Expr* solres = ds_solve(eqlist, varlist);
    Expr* fitted = NULL;
    if (solres && head_is(solres, SYM_List) && solres->data.function.arg_count >= 1) {
        Expr* branch = solres->data.function.args[0];   /* List[Rule[C[k],val],...] */
        if (head_is(branch, SYM_List))
            fitted = eval_and_free(internal_replace_all(
                (Expr*[]){ expr_copy((Expr*)body), expr_copy(branch) }, 2));
    }
    if (solres) expr_free(solres);
    if (!fitted) fitted = expr_copy((Expr*)body);   /* could not fit: keep general */
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

bool dsolve_second_order_PQ(DSolveProblem* P, Expr** Pc, Expr** Qc) {
    if (P->nfun != 1 || P->neq != 1) return false;
    if (P->max_order[0] != 2) return false;
    Expr** c; Expr* g; int n;
    if (!dsolve_linear_coeffs(P, &c, &g, &n)) return false;
    bool homog = ds_is_zero(g);
    expr_free(g);
    if (n != 2 || !homog) { for (int k = 0; k <= n; k++) expr_free(c[k]); free(c); return false; }
    /* normalized P = c1/c2, Q = c0/c2 */
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
        Expr* integ = ds_integrate(expr_copy(half), expr_new_symbol(xvar));
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
        Expr* gn = eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)g),
                        expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ expr_copy((Expr*)leadcoef), expr_new_integer(-1) }, 2)));
        Expr** ut = malloc(n * sizeof(Expr*));
        size_t uc = 0;
        bool fail = false;
        for (size_t i = 0; i < n && !fail; i++) {
            Expr* Wi = vp_matrix(dv, n, (long)i, gn);
            Expr* detWi = eval_and_free(ds_call1("Det", Wi));
            Expr* uip = eval_and_free(ds_call2(SYM_Times, detWi,
                            expr_new_function(expr_new_symbol(SYM_Power),
                                (Expr*[]){ expr_copy(detW), expr_new_integer(-1) }, 2)));
            Expr* ui = ds_integrate(uip, expr_new_symbol(xvar));
            if (ds_has_head(ui, SYM_Integrate)) { expr_free(ui); fail = true; break; }
            ut[uc++] = eval_and_free(ds_call2(SYM_Times, expr_copy(basis[i]), ui));
        }
        expr_free(gn);
        if (!fail) {
            yp = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), ut, uc));
            yp = ds_simplify(yp);
        } else { for (size_t i = 0; i < uc; i++) expr_free(ut[i]); }
        free(ut);
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

    Expr** finals = malloc(nb * sizeof(Expr*));
    size_t nf = 0;
    for (size_t b = 0; b < nb; b++) {
        if (!bodies[b]) continue;
        if (!dsolve_verify_body(P, bodies[b])) { expr_free(bodies[b]); continue; }
        finals[nf++] = dsolve_fit_constants(P, bodies[b]);
        expr_free(bodies[b]);
    }
    free(bodies);
    if (nf == 0) { free(finals); return NULL; }

    Expr* result = dsolve_assemble(P, finals, nf);
    for (size_t b = 0; b < nf; b++) expr_free(finals[b]);
    free(finals);
    return result;
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

void dsolve_fit_system(const DSolveProblem* P, Expr** bodies) {
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
    if (solres && head_is(solres, SYM_List) && solres->data.function.arg_count >= 1) {
        Expr* branch = solres->data.function.args[0];
        if (head_is(branch, SYM_List)) {
            for (size_t i = 0; i < P->nfun; i++) {
                Expr* fitted = eval_and_free(internal_replace_all(
                    (Expr*[]){ expr_copy(bodies[i]), expr_copy(branch) }, 2));
                expr_free(bodies[i]); bodies[i] = fitted;
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

/* Derivative[o1,o2][u][v1,v2] (first-order PDE terms). */
static Expr* pde_deriv_lit(const char* u, int o1, int o2, const char* v1, const char* v2) {
    Expr* d = expr_new_function(expr_new_symbol(SYM_Derivative),
                  (Expr*[]){ expr_new_integer(o1), expr_new_integer(o2) }, 2);
    Expr* du = expr_new_function(d, (Expr*[]){ expr_new_symbol(u) }, 1);
    return expr_new_function(du, (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2);
}

/* Verify a 2-variable first-order PDE.  Two obstacles: the evaluator does not
 * reduce Derivative[i,j][Function[...]][...] (so we substitute the derivative
 * TERMS with D[body, v], which does reduce), and zero_test cannot sample an
 * arbitrary function C[1][...] (so we first replace it with a concrete test
 * function C[1][z_] :> Sin[z] — a correct general solution stays a solution for
 * any choice, and Sin makes the residual concrete). */
static bool dsolve_verify_pde(const DSolveProblem* P, const Expr* body) {
    const char* u = P->fun_names[0];
    const char* v1 = P->ind_names[0];
    const char* v2 = P->ind_names[1];
    const char* z = intern_symbol("DSolve`pdez");
    Expr* blank = expr_new_function(expr_new_symbol("Blank"), NULL, 0);
    Expr* patt = expr_new_function(expr_new_symbol("Pattern"),
                     (Expr*[]){ expr_new_symbol(z), blank }, 2);
    Expr* lhs = expr_new_function(ds_const(1), (Expr*[]){ patt }, 1);          /* C[1][z_] */
    Expr* rhs = ds_call1("Sin", expr_new_symbol(z));                           /* Sin[z]   */
    Expr* rule = expr_new_function(expr_new_symbol(SYM_RuleDelayed), (Expr*[]){ lhs, rhs }, 2);
    Expr* bodyC = eval_and_free(internal_replace_all((Expr*[]){ expr_copy((Expr*)body), rule }, 2));

    bool ok = true;
    for (size_t e = 0; e < P->neq && ok; e++) {
        Expr* r = expr_copy(P->eq_residuals[e]);
        r = ds_subst(r, pde_deriv_lit(u, 1, 0, v1, v2), ds_d(expr_copy(bodyC), expr_new_symbol(v1)));
        r = ds_subst(r, pde_deriv_lit(u, 0, 1, v1, v2), ds_d(expr_copy(bodyC), expr_new_symbol(v2)));
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

Expr* dsolve_run_system(DSolveProblem* P, DSolveSysFn fn) {
    Expr** bodies = fn(P);
    if (!bodies) return NULL;
    Expr* result = NULL;
    if (dsolve_verify_system(P, bodies)) {
        dsolve_fit_system(P, bodies);
        result = dsolve_assemble_system(P, bodies);
    }
    for (size_t i = 0; i < P->nfun; i++) if (bodies[i]) expr_free(bodies[i]);
    free(bodies);
    return result;
}

Expr* dsolve_method_builtin(Expr* res, DSolveTryFn fn) {
    DSolveProblem P;
    if (!dsolve_parse(res, &P)) return NULL;
    if (P.is_pde) { dsolve_problem_free(&P); return NULL; }
    Expr* r = dsolve_run(&P, fn);
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
    Expr* Pint = ds_integrate(Pcoef, expr_new_symbol(xvar));       /* consumes Pcoef */
    if (ds_has_head(Pint, SYM_Integrate)) { expr_free(Pint); expr_free(Qcoef); return NULL; }
    Expr* mu = eval_and_free(ds_call1("Exp", Pint));
    Expr* integrand = eval_and_free(ds_call2(SYM_Times, expr_copy(mu), Qcoef)); /* consumes Qcoef */
    Expr* Qint = ds_integrate(integrand, expr_new_symbol(xvar));
    if (ds_has_head(Qint, SYM_Integrate)) { expr_free(Qint); expr_free(mu); return NULL; }
    Expr* num = eval_and_free(ds_call2(SYM_Plus, Qint, ds_const(1)));
    Expr* body = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
        num,
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ mu, expr_new_integer(-1) }, 2)
    }, 2));
    return body;
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
                    /* strip ConditionalExpression[v, cond] -> v so the body verifies */
                    if (val->type == EXPR_FUNCTION && val->data.function.arg_count == 2
                        && head_is(val, intern_symbol("ConditionalExpression")))
                        val = val->data.function.args[0];
                    out[c++] = expr_copy(val);
                }
            }
        }
    }
    *n = c;
    if (c == 0) { free(out); return NULL; }
    return out;
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
