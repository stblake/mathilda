/* dsolve_pde_ivp.c — DSolve`PDEInitialValue: initial-value problems for an
 * evolution PDE by INVARIANT FAMILIES of the initial data, plus the two
 * classical half-line boundary problems and Liouville's equation.
 *
 * An evolution PDE R[u] == 0 in time t with initial data u(x, t0) = u0(x)
 * (and u_t(x, t0) = u1(x) when R is second order in t) is solved by a family
 * Phi(x; theta(t)) that contains the data at theta(t0).  Substituting the
 * family into R, the spatial dependence must cancel; what is left is an ODE
 * (system) in t for theta, solved by the ODE cascade with theta(t0) read off
 * the data.  Three families:
 *
 *   A1  amplitude       theta(t) u0(x)                 Riemann simple wave
 *                                                      x/(1 + t); Taylor-Green
 *   A2  shift           u0(x - c (t - t0)), c const    KdV, sine-Gordon,
 *                                                      Fisher-KPP, Burgers
 *   A3  exp-polynomial  c0 Exp[sum theta_m(t) x^m]     heat and Schrodinger
 *       (the monomials of the data's exponent)         Gaussians, Lamb-Oseen
 *
 * Every answer is verified by back-substitution into the PDE and EVERY
 * condition (ds_is_zero); a family that does not close declines.  Native C;
 * Solve, DSolve, D, Simplify and ReplaceAll are called as builtins.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../numeric.h"   /* expr_is_numeric_quantity */
#include "../arithmetic.h"   /* make_rational */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define PV_MAX_IND 4
#define PV_MAX_COND 4

static Expr* pv_ev(Expr* e) { return eval_and_free(e); }
static Expr* pv_f2(const char* h, Expr* a, Expr* b) { return ds_call2(h, a, b); }
static Expr* pv_sym(const char* s) { return expr_new_symbol(s); }
static Expr* pv_list(Expr** a, size_t n) { return expr_new_function(pv_sym(SYM_List), a, n); }

static Expr* pv_assume(const DSolveProblem* P);

/* ------------------------------------------------------------- conditions */

typedef struct {
    int   fixed;     /* index of the fixed variable                       */
    Expr* at;        /* its value                                          */
    int   order;     /* derivative order in the fixed variable             */
    Expr* rhs;       /* the condition's value, in the free variables       */
} PVCond;

static bool pv_head_is(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == h;
}

/* u or Derivative[o...][u] applied to args; orders to o[] (0 for plain u). */
static bool pv_ufun(const Expr* e, const char* u, size_t nind, int* o) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != nind) return false;
    const Expr* h = e->data.function.head;
    for (size_t i = 0; i < nind; i++) o[i] = 0;
    if (h->type == EXPR_SYMBOL) return h->data.symbol.name == u;
    if (h->type != EXPR_FUNCTION || h->data.function.arg_count != 1
        || h->data.function.args[0]->type != EXPR_SYMBOL
        || h->data.function.args[0]->data.symbol.name != u) return false;
    const Expr* d = h->data.function.head;
    if (!pv_head_is(d, SYM_Derivative) || d->data.function.arg_count != nind) return false;
    for (size_t i = 0; i < nind; i++) {
        if (d->data.function.args[i]->type != EXPR_INTEGER) return false;
        o[i] = (int)d->data.function.args[i]->data.integer;
    }
    return true;
}

/* A u-funcapp in e with exactly one argument not its own variable (a point). */
static const Expr* pv_find_point(const Expr* e, const DSolveProblem* P, int* fixed, int* o) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (pv_ufun(e, P->fun_names[0], P->nind, o)) {
        int nf = 0, fi = -1;
        for (size_t i = 0; i < P->nind; i++) {
            const Expr* a = e->data.function.args[i];
            if (a->type == EXPR_SYMBOL && a->data.symbol.name == P->ind_names[i]) continue;
            bool free = true;
            for (size_t j = 0; j < P->nind; j++) free = free && ds_free_of(a, P->ind_names[j]);
            if (!free) return NULL;
            nf++; fi = (int)i;
        }
        if (nf == 1) {
            for (size_t i = 0; i < P->nind; i++) if ((int)i != fi && o[i]) return NULL;
            *fixed = fi;
            return e;
        }
        return NULL;
    }
    const Expr* r = pv_find_point(e->data.function.head, P, fixed, o);
    for (size_t i = 0; !r && i < e->data.function.arg_count; i++)
        r = pv_find_point(e->data.function.args[i], P, fixed, o);
    return r;
}

/* Split the equations into the PDE and its point conditions. */
static bool pv_classify(const DSolveProblem* P, int* pde, PVCond* c, int* nc) {
    *pde = -1; *nc = 0;
    for (size_t e = 0; e < P->neq; e++) {
        const Expr* R = P->eq_residuals[e];
        int fixed, o[PV_MAX_IND];
        const Expr* T = pv_find_point(R, P, &fixed, o);
        if (!T) {
            if (*pde >= 0) return false;              /* one PDE only */
            *pde = (int)e;
            continue;
        }
        if (*nc == PV_MAX_COND) return false;
        /* R = a T + b  ->  rhs = -b/a, a free of the variables */
        const char* s = intern_symbol("DSolve`pvT");
        Expr* Rs = ds_subst(expr_copy((Expr*)R), expr_copy((Expr*)T), pv_sym(s));
        Expr* a = ds_d(expr_copy(Rs), pv_sym(s));
        Expr* b = ds_subst(expr_copy(Rs), pv_sym(s), expr_new_integer(0));
        bool lin = ds_free_of(a, s) && !ds_is_zero(a);
        expr_free(Rs);
        if (!lin) { expr_free(a); expr_free(b); return false; }
        PVCond* k = &c[(*nc)++];
        k->fixed = fixed;
        k->at = expr_copy(T->data.function.args[fixed]);
        k->order = o[fixed];
        k->rhs = ds_simplify(pv_f2(SYM_Times, expr_new_integer(-1), pv_f2(SYM_Times, b,
                             pv_f2(SYM_Power, a, expr_new_integer(-1)))));
    }
    return *pde >= 0;
}

static void pv_free_conds(PVCond* c, int n) {
    for (int i = 0; i < n; i++) { expr_free(c[i].at); expr_free(c[i].rhs); }
}

/* ---------------------------------------------------------- substitution -- */

static Expr* pv_vars(const DSolveProblem* P) {
    Expr* v[PV_MAX_IND];
    for (size_t i = 0; i < P->nind; i++) v[i] = pv_sym(P->ind_names[i]);
    return pv_list(v, P->nind);
}

/* R /. u -> Function[vars, body], evaluated (derivatives taken). */
static Expr* pv_apply(const DSolveProblem* P, const Expr* R, const Expr* body) {
    Expr* fn = expr_new_function(pv_sym(SYM_Function),
                                 (Expr*[]){ pv_vars(P), expr_copy((Expr*)body) }, 2);
    Expr* rule = pv_f2(SYM_Rule, pv_sym(P->fun_names[0]), fn);
    return pv_ev(pv_f2(SYM_ReplaceAll, expr_copy((Expr*)R), rule));
}

/* Does body satisfy every equation (the PDE and each condition)? */
static bool pv_nonfinite(const Expr* e) {
    return ds_contains(e, SYM_ComplexInfinity) || ds_contains(e, SYM_Indeterminate)
        || ds_contains(e, SYM_DirectedInfinity) || ds_contains(e, SYM_Infinity);
}

/* A condition at a point where body is not finite (Stokes I at t = 0, where
 * the similarity variable y/Sqrt[t] is singular) holds as the one-sided limit
 * of the condition residual: body taken at fixed = at + eps, eps -> 0+. */
static bool pv_cond_limit(const DSolveProblem* P, size_t e, const Expr* body) {
    const Expr* R = P->eq_residuals[e];
    int fixed, o[PV_MAX_IND];
    const Expr* T = pv_find_point(R, P, &fixed, o);
    if (!T) return false;
    const char* eps = intern_symbol("DSolve`pvEps");
    const char* fv = P->ind_names[fixed];
    /* R with the point moved to at + eps (the residual is linear in T) */
    Expr* moved = expr_copy((Expr*)T);
    Expr* at = expr_copy(T->data.function.args[fixed]);
    expr_free(moved->data.function.args[fixed]);
    moved->data.function.args[fixed] = pv_f2(SYM_Plus, at, pv_sym(eps));
    Expr* Rm = ds_subst(expr_copy((Expr*)R), expr_copy((Expr*)T), moved);
    Expr* r = pv_apply(P, Rm, body);
    expr_free(Rm);
    Expr* opts = pv_f2(SYM_Rule, pv_sym(intern_symbol("Direction")), expr_new_string("FromAbove"));
    Expr* lim = expr_new_function(pv_sym(SYM_Limit), (Expr*[]){ r,
        pv_f2(SYM_Rule, pv_sym(eps), expr_new_integer(0)), opts,
        pv_f2(SYM_Rule, pv_sym(SYM_Assumptions), pv_assume(P)) }, 4);
    Expr* lv = pv_ev(lim);
    bool z = lv && !ds_contains(lv, SYM_Limit) && ds_is_zero(lv);
    if (lv) expr_free(lv);
    (void)fv;
    return z;
}

static bool pv_verify(const DSolveProblem* P, const Expr* body) {
    for (size_t e = 0; e < P->neq; e++) {
        Expr* r = pv_apply(P, P->eq_residuals[e], body);
        Expr* t = pv_ev(ds_call1(SYM_Together, r));
        bool z = ds_is_zero(t);
        bool nf = !z && pv_nonfinite(t);
        expr_free(t);
        if (!z && !(nf && pv_cond_limit(P, e, body))) return false;
    }
    return true;
}

/* {{u[vars] -> body}} or {{u -> Function[vars, body]}}; body consumed. */
static Expr* pv_answer(const DSolveProblem* P, Expr* body) {
    Expr* lhs; Expr* rhs;
    if (P->applied) {
        Expr* v[PV_MAX_IND];
        for (size_t i = 0; i < P->nind; i++) v[i] = pv_sym(P->ind_names[i]);
        lhs = expr_new_function(pv_sym(P->fun_names[0]), v, P->nind);
        rhs = body;
    } else {
        lhs = pv_sym(P->fun_names[0]);
        rhs = expr_new_function(pv_sym(SYM_Function), (Expr*[]){ pv_vars(P), body }, 2);
    }
    Expr* rule = pv_f2(SYM_Rule, lhs, rhs);
    return pv_list((Expr*[]){ pv_list((Expr*[]){ rule }, 1) }, 1);
}

/* Simplify under the problem's assumptions (or plainly). */
static Expr* pv_simplify(const DSolveProblem* P, Expr* e) {
    if (P->assumptions)
        return pv_ev(pv_f2(SYM_Simplify, e, expr_copy(P->assumptions)));
    return ds_simplify(e);
}

/* The ODE IVP  Derivative[n][th][t] == G,  th(t0) == a0 (, th'(t0) == a1),
 * solved by the scalar cascade; the body of th(t), or NULL. */
static Expr* pv_ode_ivp(const char* th, const char* t, int n, Expr* G, const Expr* t0,
                        const Expr* a0, const Expr* a1) {
    Expr* lhs = ds_make_funcapp(th, n, t);
    Expr* ode = pv_f2(SYM_Equal, lhs, G);
    Expr* c0 = pv_f2(SYM_Equal, expr_new_function(pv_sym(th), (Expr*[]){ expr_copy((Expr*)t0) }, 1),
                     expr_copy((Expr*)a0));
    Expr* eqs;
    if (a1) {
        Expr* d1 = expr_new_function(expr_new_function(expr_new_function(pv_sym(SYM_Derivative),
                       (Expr*[]){ expr_new_integer(1) }, 1), (Expr*[]){ pv_sym(th) }, 1),
                       (Expr*[]){ expr_copy((Expr*)t0) }, 1);
        Expr* c1 = pv_f2(SYM_Equal, d1, expr_copy((Expr*)a1));
        eqs = pv_list((Expr*[]){ ode, c0, c1 }, 3);
    } else eqs = pv_list((Expr*[]){ ode, c0 }, 2);
    Expr* sol = pv_ev(expr_new_function(pv_sym(SYM_DSolve), (Expr*[]){
        eqs, expr_new_function(pv_sym(th), (Expr*[]){ pv_sym(t) }, 1), pv_sym(t) }, 3));
    Expr* out = NULL;
    /* {{th[t] -> body}} */
    if (pv_head_is(sol, SYM_List) && sol->data.function.arg_count >= 1
        && pv_head_is(sol->data.function.args[0], SYM_List)
        && sol->data.function.args[0]->data.function.arg_count == 1) {
        const Expr* r = sol->data.function.args[0]->data.function.args[0];
        if (pv_head_is(r, SYM_Rule) && !ds_contains(r->data.function.args[1], intern_symbol("C")))
            out = expr_copy(r->data.function.args[1]);
    }
    expr_free(sol);
    return out;
}

/* The time variable and the initial data: every condition fixes the same
 * variable at the same point, with orders 0 (and 1 for second order). */
typedef struct { int tv; Expr* t0; Expr* u0; Expr* u1; int tord; } PVIvp;

static int pv_time_order(const DSolveProblem* P, const Expr* R, int tv) {
    /* highest order in tv among the u-funcapps of R */
    int best = 0, o[PV_MAX_IND];
    if (!R || R->type != EXPR_FUNCTION) return 0;
    if (pv_ufun(R, P->fun_names[0], P->nind, o)) return o[tv];
    int b = pv_time_order(P, R->data.function.head, tv);
    if (b > best) best = b;
    for (size_t i = 0; i < R->data.function.arg_count; i++) {
        b = pv_time_order(P, R->data.function.args[i], tv);
        if (b > best) best = b;
    }
    return best;
}

static bool pv_ivp(const DSolveProblem* P, int pde, const PVCond* c, int nc, PVIvp* I) {
    memset(I, 0, sizeof *I);
    if (nc < 1) return false;
    I->tv = c[0].fixed;
    for (int i = 0; i < nc; i++)
        if (c[i].fixed != I->tv || !expr_eq(c[i].at, c[0].at)) return false;
    I->tord = pv_time_order(P, P->eq_residuals[pde], I->tv);
    if (I->tord < 1 || I->tord > 2 || nc != I->tord) return false;
    for (int i = 0; i < nc; i++) {
        if (c[i].order == 0) I->u0 = c[i].rhs;
        else if (c[i].order == 1) I->u1 = c[i].rhs;
    }
    if (!I->u0 || (I->tord == 2 && !I->u1)) return false;
    I->t0 = c[0].at;
    return true;
}

/* ------------------------------------------------------------- families -- */

/* A1: theta(t) u0(x). */
static Expr* pv_amplitude(const DSolveProblem* P, const Expr* R, const PVIvp* I) {
    const char* t = P->ind_names[I->tv];
    const char* th = intern_symbol("DSolve`pvA");
    Expr* thT = expr_new_function(pv_sym(th), (Expr*[]){ pv_sym(t) }, 1);
    Expr* body = pv_ev(pv_f2(SYM_Times, thT, expr_copy(I->u0)));
    Expr* Rs = pv_apply(P, R, body);
    /* solve for the top derivative of theta */
    const char* Dn = intern_symbol("DSolve`pvD");
    Expr* top = ds_make_funcapp(th, I->tord, t);
    Expr* Rd = ds_subst(Rs, top, pv_sym(Dn));
    Expr* sol = ds_solve(pv_f2(SYM_Equal, Rd, expr_new_integer(0)), pv_sym(Dn));
    size_t ng = 0;
    Expr** g = dsolve_extract_solutions(sol, Dn, &ng);
    expr_free(sol);
    Expr* out = NULL;
    for (size_t i = 0; i < ng && !out; i++) {
        Expr* G = ds_simplify(expr_copy(g[i]));
        bool xfree = true;
        for (size_t v = 0; v < P->nind && xfree; v++)
            if ((int)v != I->tv) xfree = ds_free_of(G, P->ind_names[v]);
        if (!xfree) { expr_free(G); continue; }
        Expr* a1 = NULL;
        if (I->tord == 2) {                       /* u1 = theta'(t0) u0 */
            a1 = ds_simplify(pv_f2(SYM_Times, expr_copy(I->u1), pv_f2(SYM_Power, expr_copy(I->u0), expr_new_integer(-1))));
            bool cst = true;
            for (size_t v = 0; v < P->nind && cst; v++)
                if ((int)v != I->tv) cst = ds_free_of(a1, P->ind_names[v]);
            if (!cst) { expr_free(a1); expr_free(G); continue; }
        }
        Expr* one = expr_new_integer(1);
        Expr* thb = pv_ode_ivp(th, t, I->tord, G, I->t0, one, a1);
        expr_free(one);
        if (a1) expr_free(a1);
        if (!thb) continue;
        out = pv_ev(pv_f2(SYM_Times, thb, expr_copy(I->u0)));
    }
    for (size_t i = 0; i < ng; i++) expr_free(g[i]);
    free(g);
    expr_free(body);
    return out;
}

/* A2: u0(x - c (t - t0)) with one spatial variable. */
static Expr* pv_shift(const DSolveProblem* P, const Expr* R, const PVIvp* I) {
    if (P->nind != 2) return NULL;
    int xv = 1 - I->tv;
    const char* x = P->ind_names[xv];
    const char* t = P->ind_names[I->tv];
    if (ds_free_of(I->u0, x)) return NULL;
    const char* cs = intern_symbol("DSolve`pvC");
    Expr* arg = pv_f2(SYM_Plus, pv_sym(x), pv_f2(SYM_Times, expr_new_integer(-1),
                      pv_f2(SYM_Times, pv_sym(cs), pv_f2(SYM_Plus, pv_sym(t),
                            pv_f2(SYM_Times, expr_new_integer(-1), expr_copy(I->t0))))));
    Expr* body = ds_subst(expr_copy(I->u0), pv_sym(x), arg);
    Expr* c = NULL;
    if (I->tord == 2) {
        /* u_t(x, t0) = -c u0'(x) */
        c = pv_simplify(P, pv_f2(SYM_Times, expr_new_integer(-1),
                pv_f2(SYM_Times, expr_copy(I->u1), pv_f2(SYM_Power,
                      ds_d(expr_copy(I->u0), pv_sym(x)), expr_new_integer(-1)))));
        if (!ds_free_of(c, x)) { expr_free(c); c = NULL; }
    } else {
        /* in the moving coordinate xi = x - c (t - t0) the profile's argument
         * is c-free, so c is left only as the coefficient u_t produces */
        const char* xi = intern_symbol("DSolve`pvXi");
        Expr* Rs = pv_apply(P, R, body);
        Rs = ds_subst(Rs, pv_sym(x), pv_f2(SYM_Plus, pv_sym(xi), pv_f2(SYM_Times, pv_sym(cs),
                      pv_f2(SYM_Plus, pv_sym(t), pv_f2(SYM_Times, expr_new_integer(-1), expr_copy(I->t0))))));
        Expr* sol = ds_solve(pv_f2(SYM_Equal, Rs, expr_new_integer(0)), pv_sym(cs));
        size_t ng = 0;
        Expr** g = dsolve_extract_solutions(sol, cs, &ng);
        expr_free(sol);
        for (size_t i = 0; i < ng; i++) {
            Expr* G = g[i];
            if (!c && ds_free_of(G, xi) && ds_free_of(G, t)) c = pv_simplify(P, expr_copy(G));
            expr_free(G);
        }
        free(g);
        if (c && (!ds_free_of(c, xi) || !ds_free_of(c, t))) { expr_free(c); c = NULL; }
    }
    Expr* out = NULL;
    if (c) out = ds_subst(body, pv_sym(cs), c);
    else expr_free(body);
    return out;
}

/* A3: c0 Exp[sum theta_m(t) mono_m] over the monomials of u0's exponent. */
static bool pv_exp_split(const Expr* u0, Expr** c0, Expr** expo) {
    /* u0 = c0 E^P (c0 free of the spatial variables is checked by the caller) */
    if (pv_head_is(u0, SYM_Power) && u0->data.function.args[0]->type == EXPR_SYMBOL
        && u0->data.function.args[0]->data.symbol.name == SYM_E) {
        *c0 = expr_new_integer(1);
        *expo = expr_copy(u0->data.function.args[1]);
        return true;
    }
    if (!pv_head_is(u0, SYM_Times)) return false;
    size_t n = u0->data.function.arg_count;
    Expr** rest = malloc(sizeof(Expr*) * n);
    Expr* P = NULL; size_t w = 0;
    for (size_t i = 0; i < n; i++) {
        const Expr* f = u0->data.function.args[i];
        if (!P && pv_head_is(f, SYM_Power) && f->data.function.args[0]->type == EXPR_SYMBOL
            && f->data.function.args[0]->data.symbol.name == SYM_E)
            P = expr_copy(f->data.function.args[1]);
        else rest[w++] = expr_copy((Expr*)f);
    }
    if (!P) { for (size_t i = 0; i < w; i++) expr_free(rest[i]); free(rest); return false; }
    *c0 = pv_ev(expr_new_function(pv_sym(SYM_Times), rest, w));
    free(rest);
    *expo = P;
    return true;
}

#define PV_MAX_MONO 8
static Expr* pv_exppoly(const DSolveProblem* P, const Expr* R, const PVIvp* I) {
    if (I->tord != 1) return NULL;
    const char* t = P->ind_names[I->tv];
    Expr *c0 = NULL, *ex = NULL;
    if (!pv_exp_split(I->u0, &c0, &ex)) return NULL;
    Expr* out = NULL;
    /* monomials of Expand[exponent]: each term = coefficient * spatial monomial */
    Expr* E = pv_ev(ds_call1(SYM_Expand, ex));
    size_t nt = pv_head_is(E, SYM_Plus) ? E->data.function.arg_count : 1;
    Expr* const* tm = pv_head_is(E, SYM_Plus) ? E->data.function.args : &E;
    Expr* mono[PV_MAX_MONO]; Expr* init[PV_MAX_MONO]; int nm = 0;
    bool ok = nt + 1 <= PV_MAX_MONO;
    /* the constant monomial first (it carries no data: theta_0(t0) = 0) */
    mono[nm] = expr_new_integer(1); init[nm] = expr_new_integer(0); nm++;
    for (size_t i = 0; i < nt && ok; i++) {
        /* coefficient = term with every spatial variable set to 1 */
        Expr* coef = expr_copy(tm[i]);
        Expr* m = expr_copy(tm[i]);
        bool spatial = false;
        for (size_t v = 0; v < P->nind; v++) {
            if ((int)v == I->tv) continue;
            if (!ds_free_of(tm[i], P->ind_names[v])) spatial = true;
            coef = ds_subst(coef, pv_sym(P->ind_names[v]), expr_new_integer(1));
        }
        if (!spatial) {                       /* a constant term: into theta_0 */
            Expr* nv = pv_ev(pv_f2(SYM_Plus, init[0], coef));
            init[0] = nv;
            expr_free(m);
            continue;
        }
        Expr* mm = pv_ev(pv_f2(SYM_Times, m, pv_f2(SYM_Power, expr_copy(coef), expr_new_integer(-1))));
        /* the monomial must be free of t and of parameters: a product of powers */
        if (!ds_free_of(mm, t)) { expr_free(mm); expr_free(coef); ok = false; break; }
        mono[nm] = mm; init[nm] = coef; nm++;
    }
    expr_free(E);
    if (ok) {
        /* Phi = c0 Exp[sum th_m[t] mono_m] */
        const char* th[PV_MAX_MONO];
        Expr** sum = malloc(sizeof(Expr*) * (size_t)nm);
        for (int m = 0; m < nm; m++) {
            char nm_[32]; snprintf(nm_, sizeof nm_, "DSolve`pvE%d", m);
            th[m] = intern_symbol(nm_);
            sum[m] = pv_f2(SYM_Times, expr_new_function(pv_sym(th[m]), (Expr*[]){ pv_sym(t) }, 1),
                           expr_copy(mono[m]));
        }
        Expr* Phi = pv_ev(pv_f2(SYM_Times, expr_copy(c0),
                         pv_f2(SYM_Power, pv_sym(SYM_E), expr_new_function(pv_sym(SYM_Plus), sum, (size_t)nm))));
        free(sum);
        Expr* Rs = pv_apply(P, R, Phi);
        Expr* Q = pv_ev(ds_call1(SYM_Expand, pv_f2(SYM_Times, expr_copy(Rs),
                        pv_f2(SYM_Power, expr_copy(Phi), expr_new_integer(-1)))));
        /* coefficient of mono_m: Q is a sum of (coefficient) * (spatial monomial);
         * the coefficient of mono_m is read by substituting a symbol for it via
         * the theta_m' term it must carry: solve Q's monomial components by
         * matching each term's spatial part. */
        size_t qn = pv_head_is(Q, SYM_Plus) ? Q->data.function.arg_count : 1;
        Expr* const* qt = pv_head_is(Q, SYM_Plus) ? Q->data.function.args : &Q;
        Expr** grp = calloc((size_t)nm, sizeof(Expr*));
        bool closes = true;
        for (size_t i = 0; i < qn && closes; i++) {
            Expr* coef = expr_copy(qt[i]);
            for (size_t v = 0; v < P->nind; v++)
                if ((int)v != I->tv) coef = ds_subst(coef, pv_sym(P->ind_names[v]), expr_new_integer(1));
            if (ds_is_zero(coef)) { expr_free(coef); continue; }
            Expr* sp = pv_ev(pv_f2(SYM_Times, expr_copy(qt[i]), pv_f2(SYM_Power, expr_copy(coef), expr_new_integer(-1))));
            int hit = -1;
            for (int m = 0; m < nm && hit < 0; m++) if (expr_eq(sp, mono[m])) hit = m;
            /* a monomial outside the family (a Laurent r^-1 term, an x^3): the
             * family closes only if such groups cancel -- left to the final
             * verification; a non-monomial spatial part declines. */
            bool ismono = true;
            for (size_t v = 0; v < P->nind && ismono; v++)
                if ((int)v != I->tv) ismono = ds_free_of(coef, P->ind_names[v]);
            if (!ismono || ds_contains(sp, SYM_E)) closes = false;
            else if (hit >= 0) {
                grp[hit] = grp[hit] ? pv_ev(pv_f2(SYM_Plus, grp[hit], coef)) : coef;
                coef = NULL;
            }
            expr_free(sp);
            if (coef) expr_free(coef);
        }
        /* solve the groups triangularly: each grp[m] == 0 for theta_m' */
        Expr** solved = calloc((size_t)nm, sizeof(Expr*));
        int nsolved = 0;
        for (int pass = 0; pass < nm && closes && nsolved < nm; pass++) {
            for (int m = 0; m < nm && closes; m++) {
                if (solved[m]) continue;
                if (!grp[m]) { closes = false; break; }
                /* substitute the solved thetas */
                Expr* g = expr_copy(grp[m]);
                for (int k = 0; k < nm; k++) if (solved[k]) {
                    Expr* fn = expr_new_function(pv_sym(SYM_Function), (Expr*[]){
                        pv_list((Expr*[]){ pv_sym(t) }, 1), expr_copy(solved[k]) }, 2);
                    g = pv_ev(pv_f2(SYM_ReplaceAll, g, pv_f2(SYM_Rule, pv_sym(th[k]), fn)));
                }
                /* depends only on theta_m (and t)? */
                bool own = true;
                for (int k = 0; k < nm && own; k++) if (k != m && !solved[k] && ds_contains(g, th[k])) own = false;
                if (!own) { expr_free(g); continue; }
                const char* Dn = intern_symbol("DSolve`pvD");
                Expr* top = ds_make_funcapp(th[m], 1, t);
                Expr* gd = ds_subst(g, top, pv_sym(Dn));
                Expr* sol = ds_solve(pv_f2(SYM_Equal, gd, expr_new_integer(0)), pv_sym(Dn));
                size_t ng = 0;
                Expr** gs = dsolve_extract_solutions(sol, Dn, &ng);
                expr_free(sol);
                if (ng == 1) solved[m] = pv_ode_ivp(th[m], t, 1, expr_copy(gs[0]), I->t0, init[m], NULL);
                for (size_t q = 0; q < ng; q++) expr_free(gs[q]);
                free(gs);
                if (!solved[m]) { closes = false; break; }
                nsolved++;
            }
        }
        if (closes && nsolved == nm) {
            Expr** s2 = malloc(sizeof(Expr*) * (size_t)nm);
            for (int m = 0; m < nm; m++) s2[m] = pv_f2(SYM_Times, expr_copy(solved[m]), expr_copy(mono[m]));
            out = pv_ev(pv_f2(SYM_Times, expr_copy(c0), pv_f2(SYM_Power, pv_sym(SYM_E),
                        expr_new_function(pv_sym(SYM_Plus), s2, (size_t)nm))));
            free(s2);
        }
        for (int m = 0; m < nm; m++) { if (grp[m]) expr_free(grp[m]); if (solved[m]) expr_free(solved[m]); }
        free(grp); free(solved);
        expr_free(Q); expr_free(Rs); expr_free(Phi);
    }
    for (int m = 0; m < nm; m++) { expr_free(mono[m]); expr_free(init[m]); }
    expr_free(c0);
    return out;
}

/* --------------------------------------------- boundary-driven problems -- */

/* The assumptions for the boundary-driven methods: the user's, or -- with none
 * given -- every parameter of the problem positive: the generic, well-posed
 * case (a forward diffusion needs nu > 0; a frequency is taken positive), as
 * a general solution is generic in its parameters -- and the variables on the
 * half-line domain y > 0, t > 0.  Owned. */
/* The parameter symbols of e: leaves that are not numeric constants (Pi, E),
 * not the unknown and not a variable (Variables would treat Cos[omega t] as
 * an atom and miss omega). */
static void pv_params(const DSolveProblem* P, const Expr* e, const char** out, size_t* n, size_t cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        if (nm == P->fun_names[0] || expr_is_numeric_quantity(e)) return;
        for (size_t k = 0; k < P->nind; k++) if (nm == P->ind_names[k]) return;
        for (size_t i = 0; i < *n; i++) if (out[i] == nm) return;
        if (*n < cap) out[(*n)++] = nm;
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    for (size_t i = 0; i < e->data.function.arg_count; i++) pv_params(P, e->data.function.args[i], out, n, cap);
}

static Expr* pv_assume(const DSolveProblem* P) {
    if (P->assumptions) return expr_copy(P->assumptions);
    const char* ps[32]; size_t np = 0;
    for (size_t e = 0; e < P->neq; e++) pv_params(P, P->eq_residuals[e], ps, &np, 32);
    Expr** c = malloc(sizeof(Expr*) * (np + P->nind + 1));
    size_t w = 0;
    /* the variables on the half-line domain, then the parameters */
    for (size_t k = 0; k < P->nind; k++)
        c[w++] = pv_f2(SYM_Greater, pv_sym(P->ind_names[k]), expr_new_integer(0));
    for (size_t i = 0; i < np; i++) c[w++] = pv_f2(SYM_Greater, pv_sym(ps[i]), expr_new_integer(0));
    Expr* out = expr_new_function(pv_sym(SYM_And), c, w);
    free(c);
    return out;
}

/* A ConditionalExpression-free value of the general solution's rule list. */
static Expr* pv_dsolve_general(Expr* ode, const char* F, const char* s) {
    Expr* sol = pv_ev(expr_new_function(pv_sym(SYM_DSolve), (Expr*[]){
        ode, expr_new_function(pv_sym(F), (Expr*[]){ pv_sym(s) }, 1), pv_sym(s) }, 3));
    Expr* out = NULL;
    if (pv_head_is(sol, SYM_List) && sol->data.function.arg_count == 1
        && pv_head_is(sol->data.function.args[0], SYM_List)
        && sol->data.function.args[0]->data.function.arg_count == 1) {
        const Expr* r = sol->data.function.args[0]->data.function.args[0];
        if (pv_head_is(r, SYM_Rule)) out = expr_copy(r->data.function.args[1]);
    }
    expr_free(sol);
    return out;
}

static Expr* pv_Ck(int k) { return expr_new_function(pv_sym(intern_symbol("C")), (Expr*[]){ expr_new_integer(k) }, 1); }

/* B: similarity.  u_t = L[u] invariant under (y, t) -> (lam y, lam^2 t), with
 * constant data: u(y, t0) = a, u(y0, t) = b.  u = F((y - y0)/Sqrt[t - t0]);
 * the PDE becomes an ODE in eta; F(0) = b, F(Infinity) = a. */
static Expr* pv_similarity(const DSolveProblem* P, const Expr* R, const PVCond* c, int nc) {
    if (P->nind != 2 || nc != 2) return NULL;
    const PVCond *ic = NULL, *bc = NULL;
    for (int i = 0; i < 2; i++) {
        if (c[i].order != 0) return NULL;
        if (!ic && pv_time_order(P, R, c[i].fixed) == 1) ic = &c[i];
        else bc = &c[i];
    }
    if (!ic || !bc || ic->fixed == bc->fixed) return NULL;
    const char* t = P->ind_names[ic->fixed];
    const char* y = P->ind_names[bc->fixed];
    if (!ds_free_of(ic->rhs, y) || !ds_free_of(ic->rhs, t) || !ds_free_of(bc->rhs, t)) return NULL;
    const char* F = intern_symbol("DSolve`pvF");
    const char* eta = intern_symbol("DSolve`pvEta");
    Expr* tau = pv_f2(SYM_Plus, pv_sym(t), pv_f2(SYM_Times, expr_new_integer(-1), expr_copy(ic->at)));
    Expr* et = pv_f2(SYM_Times, pv_f2(SYM_Plus, pv_sym(y), pv_f2(SYM_Times, expr_new_integer(-1), expr_copy(bc->at))),
                     pv_f2(SYM_Power, expr_copy(tau), make_rational(-1, 2)));
    Expr* body = expr_new_function(pv_sym(F), (Expr*[]){ et }, 1);
    Expr* Rs = pv_apply(P, R, body);
    /* y = y0 + eta Sqrt[tau]; then tau^k Rs must be tau-free */
    Rs = ds_subst(Rs, pv_sym(y), pv_f2(SYM_Plus, expr_copy(bc->at),
                  pv_f2(SYM_Times, pv_sym(eta), pv_f2(SYM_Power, expr_copy(tau), make_rational(1, 2)))));
    Expr* G = NULL;
    static const int64_t kn[4] = { 1, 1, 3, 2 }, kd[4] = { 1, 2, 2, 1 };
    for (int k = 0; k < 4 && !G; k++) {
        Expr* g = pv_ev(ds_call1(SYM_Together, pv_f2(SYM_Times, expr_copy(Rs),
                        pv_f2(SYM_Power, expr_copy(tau), make_rational(kn[k], kd[k])))));
        if (ds_free_of(g, t) && ds_contains(g, F)) G = g; else expr_free(g);
    }
    expr_free(Rs);
    Expr* out = NULL;
    if (G) {
        Expr* gen = pv_dsolve_general(pv_f2(SYM_Equal, expr_copy(G), expr_new_integer(0)), F, eta);
        if (gen && ds_contains(gen, intern_symbol("C"))) {
            /* F(0) = b,  F(Infinity) = a  (linear in C[1], C[2]) */
            Expr* f0 = ds_subst(expr_copy(gen), pv_sym(eta), expr_new_integer(0));
            Expr* finf = pv_ev(expr_new_function(pv_sym(SYM_Limit), (Expr*[]){ expr_copy(gen),
                pv_f2(SYM_Rule, pv_sym(eta), pv_sym(SYM_Infinity)),
                pv_f2(SYM_Rule, pv_sym(SYM_Assumptions), pv_assume(P)) }, 3));
            if (!ds_contains(finf, SYM_Limit)) {
                Expr* eqs = pv_list((Expr*[]){ pv_f2(SYM_Equal, f0, expr_copy(bc->rhs)),
                                               pv_f2(SYM_Equal, finf, expr_copy(ic->rhs)) }, 2);
                finf = NULL;
                Expr* cs = pv_list((Expr*[]){ pv_Ck(1), pv_Ck(2) }, 2);
                Expr* sol = ds_solve(eqs, cs);
                if (pv_head_is(sol, SYM_List) && sol->data.function.arg_count == 1) {
                    Expr* fit = pv_ev(pv_f2(SYM_ReplaceAll, expr_copy(gen), expr_copy(sol->data.function.args[0])));
                    if (!ds_contains(fit, intern_symbol("C"))) {
                        Expr* e2 = pv_f2(SYM_Times, pv_f2(SYM_Plus, pv_sym(y), pv_f2(SYM_Times, expr_new_integer(-1), expr_copy(bc->at))),
                                         pv_f2(SYM_Power, expr_copy(tau), make_rational(-1, 2)));
                        out = ds_subst(fit, pv_sym(eta), e2);
                    } else expr_free(fit);
                }
                expr_free(sol);
            } else { expr_free(f0); }
            if (finf) expr_free(finf);
        }
        if (gen) expr_free(gen);
        expr_free(G);
    }
    expr_free(body); expr_free(tau);
    return out;
}

/* Every free symbol of e other than the variables is proved positive under
 * the assumptions (then PowerExpand is exact: Sqrt[4 I nu w] = 2 Sqrt[I]
 * Sqrt[nu] Sqrt[w]). */
static bool pv_params_positive(const DSolveProblem* P, const Expr* e) {
    Expr* asm_ = pv_assume(P);
    Expr* vs = pv_ev(ds_call1(intern_symbol("Variables"), expr_copy((Expr*)e)));
    bool ok = pv_head_is(vs, SYM_List);
    for (size_t i = 0; ok && i < vs->data.function.arg_count; i++) {
        const Expr* v = vs->data.function.args[i];
        bool isvar = false;
        for (size_t k = 0; k < P->nind; k++)
            isvar = isvar || (v->type == EXPR_SYMBOL && v->data.symbol.name == P->ind_names[k]);
        if (isvar) continue;
        if (v->type != EXPR_SYMBOL) { ok = false; break; }
        Expr* q = pv_ev(pv_f2(SYM_Simplify, pv_f2(SYM_Greater, expr_copy((Expr*)v), expr_new_integer(0)),
                              expr_copy(asm_)));
        ok = q->type == EXPR_SYMBOL && q->data.symbol.name == SYM_True;
        expr_free(q);
    }
    expr_free(vs); expr_free(asm_);
    return ok;
}

/* Re of a constant exponent rate, its parameters set to 13/10, 17/10, ...:
 * the sign decides decay.  2 for undecided. */
static int pv_re_sign(const Expr* lam) {
    Expr* vs = pv_ev(ds_call1(intern_symbol("Variables"), expr_copy((Expr*)lam)));
    Expr* v = pv_ev(ds_call1(SYM_PowerExpand, expr_copy((Expr*)lam)));
    if (pv_head_is(vs, SYM_List))
        for (size_t i = 0; i < vs->data.function.arg_count; i++)
            v = ds_subst(v, expr_copy(vs->data.function.args[i]), make_rational(13 + 4 * (int64_t)i, 10));
    expr_free(vs);
    Expr* n = pv_ev(ds_call1(SYM_N, ds_call1(SYM_Re, v)));
    int sg = 2;
    if (n && n->type == EXPR_REAL) sg = n->data.real < -1e-12 ? -1 : (n->data.real > 1e-12 ? 1 : 0);
    if (n) expr_free(n);
    return sg;
}

/* C: time-harmonic boundary data, u(y0, t) = Re[A e^(i w t)], no initial
 * condition: the periodic steady state u = Re[e^(i w t) F(y)], F the bounded
 * solution of the reduced ODE with F(y0) = A. */
static Expr* pv_harmonic(const DSolveProblem* P, const Expr* R, const PVCond* c, int nc) {
    if (P->nind != 2 || nc != 1 || c[0].order != 0) return NULL;
    const char* y = P->ind_names[c[0].fixed];
    const char* t = P->ind_names[1 - c[0].fixed];
    const Expr* g = c[0].rhs;
    if (!ds_free_of(g, y) || ds_free_of(g, t)) return NULL;
    /* w^2 = -g''/g constant */
    Expr* asmH = pv_assume(P);
    Expr* w2 = pv_ev(pv_f2(SYM_Simplify, pv_f2(SYM_Times, expr_new_integer(-1),
        pv_f2(SYM_Times, ds_d(ds_d(expr_copy((Expr*)g), pv_sym(t)), pv_sym(t)),
              pv_f2(SYM_Power, expr_copy((Expr*)g), expr_new_integer(-1)))), expr_copy(asmH)));
    if (!ds_free_of(w2, t) || ds_is_zero(w2)) { expr_free(w2); expr_free(asmH); return NULL; }
    Expr* w = pv_ev(pv_f2(SYM_Simplify, pv_f2(SYM_Power, w2, make_rational(1, 2)), expr_copy(asmH)));
    /* A = g(0) - i g'(0)/w */
    Expr* g0 = ds_subst(expr_copy((Expr*)g), pv_sym(t), expr_new_integer(0));
    Expr* g1 = ds_subst(ds_d(expr_copy((Expr*)g), pv_sym(t)), pv_sym(t), expr_new_integer(0));
    Expr* A = pv_ev(pv_f2(SYM_Plus, g0, pv_f2(SYM_Times, pv_f2(SYM_Times, expr_new_integer(-1), pv_sym(SYM_I)),
                     pv_f2(SYM_Times, g1, pv_f2(SYM_Power, expr_copy(w), expr_new_integer(-1))))));
    const char* F = intern_symbol("DSolve`pvH");
    Expr* ex = pv_f2(SYM_Power, pv_sym(SYM_E), pv_f2(SYM_Times, pv_sym(SYM_I), pv_f2(SYM_Times, expr_copy(w), pv_sym(t))));
    Expr* body = pv_f2(SYM_Times, expr_copy(ex), expr_new_function(pv_sym(F), (Expr*[]){ pv_sym(y) }, 1));
    Expr* Rs = pv_apply(P, R, body);
    Expr* G = pv_ev(ds_call1(SYM_Together, pv_f2(SYM_Times, Rs, pv_f2(SYM_Power, expr_copy(ex), expr_new_integer(-1)))));
    Expr* out = NULL;
    if (ds_free_of(G, t)) {
        Expr* gen = pv_dsolve_general(pv_f2(SYM_Equal, expr_copy(G), expr_new_integer(0)), F, y);
        if (gen) {
            /* keep the basis members that decay as y -> Infinity */
            Expr* keep = expr_copy(gen);
            for (int k = 1; k <= 2; k++) {
                Expr* only = pv_ev(pv_f2(SYM_ReplaceAll, expr_copy(gen), pv_list((Expr*[]){
                    pv_f2(SYM_Rule, pv_Ck(k), expr_new_integer(1)),
                    pv_f2(SYM_Rule, pv_Ck(3 - k), expr_new_integer(0)) }, 2)));
                /* a basis member E^(lam y): it decays iff Re[lam] < 0 */
                Expr* lam = ds_simplify(pv_f2(SYM_Times, ds_d(expr_copy(only), pv_sym(y)),
                                              pv_f2(SYM_Power, expr_copy(only), expr_new_integer(-1))));
                bool decays = ds_free_of(lam, y) && pv_params_positive(P, lam) && pv_re_sign(lam) == -1;
                expr_free(lam); expr_free(only);
                if (!decays) keep = pv_ev(pv_f2(SYM_ReplaceAll, keep, pv_f2(SYM_Rule, pv_Ck(k), expr_new_integer(0))));
            }
            /* F(y0) = A for the remaining constant */
            Expr* f0 = ds_subst(expr_copy(keep), pv_sym(y), expr_copy(c[0].at));
            for (int k = 1; k <= 2 && !out; k++) {
                if (!ds_contains(keep, intern_symbol("C"))) break;
                Expr* sol = ds_solve(pv_f2(SYM_Equal, expr_copy(f0), expr_copy(A)), pv_Ck(k));
                if (pv_head_is(sol, SYM_List) && sol->data.function.arg_count == 1) {
                    Expr* fit = pv_ev(pv_f2(SYM_ReplaceAll, expr_copy(keep), expr_copy(sol->data.function.args[0])));
                    if (!ds_contains(fit, intern_symbol("C"))) {
                        Expr* uc = pv_ev(ds_call1(SYM_PowerExpand, pv_f2(SYM_Times, expr_copy(ex), fit)));
                        out = pv_ev(ds_call1(SYM_ComplexExpand, ds_call1(SYM_Re, uc)));
                    } else expr_free(fit);
                }
                expr_free(sol);
            }
            expr_free(f0); expr_free(keep); expr_free(gen);
        }
    }
    expr_free(G); expr_free(body); expr_free(ex); expr_free(A); expr_free(w);
    if (out) out = pv_ev(pv_f2(SYM_Simplify, out, expr_copy(asmH)));
    expr_free(asmH);
    return out;
}

/* D: Liouville's equation a u_xy + b E^(m u) == 0 (any constants):
 * u = (1/m) Log[2 C[1]'(x) C[2]'(y) / (lam m (C[1](x) + C[2](y))^2)], lam = -b/a. */
static Expr* pv_liouville(const DSolveProblem* P) {
    if (P->nind != 2 || P->neq != 1) return NULL;
    const char* u = P->fun_names[0];
    const char* x = P->ind_names[0];
    const char* y = P->ind_names[1];
    const Expr* R = P->eq_residuals[0];
    Expr* uxy = expr_new_function(expr_new_function(expr_new_function(pv_sym(SYM_Derivative),
                    (Expr*[]){ expr_new_integer(1), expr_new_integer(1) }, 2), (Expr*[]){ pv_sym(u) }, 1),
                    (Expr*[]){ pv_sym(x), pv_sym(y) }, 2);
    const char* s = intern_symbol("DSolve`pvT");
    Expr* Rs = ds_subst(expr_copy((Expr*)R), uxy, pv_sym(s));
    Expr* a = ds_d(expr_copy(Rs), pv_sym(s));
    Expr* rest = ds_subst(Rs, pv_sym(s), expr_new_integer(0));
    Expr* out = NULL;
    /* rest = b E^(m u[x, y]) */
    Expr* uf = expr_new_function(pv_sym(u), (Expr*[]){ pv_sym(x), pv_sym(y) }, 2);
    const char* v = intern_symbol("DSolve`pvU");
    Expr* rv = ds_subst(rest, uf, pv_sym(v));
    if (!ds_contains(rv, u) && !ds_is_zero(a) && ds_free_of(a, x) && ds_free_of(a, y)
        && ds_free_of(rv, x) && ds_free_of(rv, y)) {
        /* m = d/dv log(rv) must be a nonzero constant */
        Expr* m = ds_simplify(pv_f2(SYM_Times, ds_d(expr_copy(rv), pv_sym(v)),
                              pv_f2(SYM_Power, expr_copy(rv), expr_new_integer(-1))));
        if (ds_free_of(m, v) && !ds_is_zero(m)) {
            Expr* b = ds_simplify(pv_f2(SYM_Times, expr_copy(rv), pv_f2(SYM_Power, pv_sym(SYM_E),
                                  pv_f2(SYM_Times, expr_new_integer(-1), pv_f2(SYM_Times, expr_copy(m), pv_sym(v))))));
            if (ds_free_of(b, v)) {
                Expr* lam = pv_ev(pv_f2(SYM_Times, expr_new_integer(-1), pv_f2(SYM_Times, b,
                                  pv_f2(SYM_Power, expr_copy(a), expr_new_integer(-1)))));
                Expr* F1 = expr_new_function(pv_Ck(1), (Expr*[]){ pv_sym(x) }, 1);
                Expr* F2 = expr_new_function(pv_Ck(2), (Expr*[]){ pv_sym(y) }, 1);
                Expr* d1 = ds_d(expr_copy(F1), pv_sym(x));
                Expr* d2 = ds_d(expr_copy(F2), pv_sym(y));
                Expr* arg = pv_f2(SYM_Times, pv_f2(SYM_Times, expr_new_integer(2), pv_f2(SYM_Times, d1, d2)),
                    pv_f2(SYM_Power, pv_f2(SYM_Times, pv_f2(SYM_Times, lam, expr_copy(m)),
                        pv_f2(SYM_Power, pv_f2(SYM_Plus, F1, F2), expr_new_integer(2))), expr_new_integer(-1)));
                out = pv_ev(pv_f2(SYM_Times, pv_f2(SYM_Power, expr_copy(m), expr_new_integer(-1)), ds_call1(SYM_Log, arg)));
            }
        }
        expr_free(m);
    }
    expr_free(rv); expr_free(a);
    return out;
}

/* ------------------------------------------------------------- the runner -- */

Expr* dsolve_pde_ivp_run(DSolveProblem* P);
Expr* dsolve_pde_ivp_run(DSolveProblem* P) {
    if (!P->is_pde || P->nfun != 1 || P->nind < 2 || P->nind > PV_MAX_IND) return NULL;
    if (P->neq == 1) {
        Expr* body = pv_liouville(P);
        if (!body) return NULL;
        if (pv_verify(P, body)) return pv_answer(P, body);
        expr_free(body);
        return NULL;
    }
    PVCond c[PV_MAX_COND]; int nc = 0, pde = -1;
    if (!pv_classify(P, &pde, c, &nc)) { pv_free_conds(c, nc); return NULL; }
    Expr* result = NULL;
    PVIvp I;
    if (pv_ivp(P, pde, c, nc, &I)) {
        const Expr* R = P->eq_residuals[pde];
        Expr* (*fam[3])(const DSolveProblem*, const Expr*, const PVIvp*) = { pv_exppoly, pv_amplitude, pv_shift };
        for (int f = 0; f < 3 && !result; f++) {
            Expr* body = fam[f](P, R, &I);
            if (!body) continue;
            Expr* sb = pv_simplify(P, body);
            if (pv_verify(P, sb)) result = pv_answer(P, sb);
            else expr_free(sb);
        }
    } else {
        const Expr* R = P->eq_residuals[pde];
        Expr* body = pv_similarity(P, R, c, nc);
        if (!body) body = pv_harmonic(P, R, c, nc);
        if (body) {
            Expr* sb = pv_ev(pv_f2(SYM_Simplify, body, pv_assume(P)));
            if (pv_verify(P, sb)) result = pv_answer(P, sb);
            else expr_free(sb);
        }
    }
    pv_free_conds(c, nc);
    return result;
}

static Expr* builtin_dsolve_pde_ivp(Expr* res) {
    DSolveProblem P;
    if (!dsolve_parse(res, &P)) { dsolve_problem_free(&P); return NULL; }
    Expr* r = dsolve_pde_ivp_run(&P);
    dsolve_problem_free(&P);
    return r;
}

void dsolve_pde_ivp_init(void);
void dsolve_pde_ivp_init(void) {
    symtab_add_builtin("DSolve`PDEInitialValue", builtin_dsolve_pde_ivp);
    symtab_get_def("DSolve`PDEInitialValue")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`PDEInitialValue",
        "DSolve`PDEInitialValue[{pde, ic..}, u, {x.., t}] solves the initial-value problem of an "
        "evolution PDE by an invariant family of the data -- amplitude theta(t) u0(x), shift "
        "u0(x - c t), or exponential polynomial Exp[sum theta_m(t) x^m] -- reducing it to an ODE "
        "in t; every answer is verified by back-substitution.");
}
