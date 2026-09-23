/*
 * dsolve_solvefor.c — DSolve`SolvableForY / DSolve`SolvableForX.
 *
 * The differentiation ("dp") method (Maple's dp/dp2): an ODE F(x,y,y')=0 that is
 * algebraic in y (SolvableForY) or in x (SolvableForX) but not directly
 * separable/linear/exact/... is solved by isolating that variable and
 * differentiating.  With p = y':
 *
 *   SolvableForY.  Solve F==0 for y  ->  y = G(x,p).  Differentiate w.r.t. x
 *   (p = dy/dx):
 *        p = G_x + G_p p'      =>      dx/dp = G_p / (p - G_x).
 *   Recurse DSolve on this first-order ODE for x(p) -> x = X(p,C); then
 *   Y = G(X(p,C), p), giving the PARAMETRIC general solution {x=X(p), y=Y(p)}
 *   with the slope p = y' as the parameter.
 *
 *   SolvableForX.  Solve F==0 for x  ->  x = H(y,p) (dx/dy = 1/p).  Differentiate
 *   w.r.t. y:
 *        1/p = H_y + H_p p'(y)  =>      dy/dp = H_p / (1/p - H_y).
 *   Recurse DSolve for y(p) -> y = Y(p,C); then X = H(Y(p,C), p).
 *
 * Both emit the branch through DSolve`Param[X, Y, p]; the parametric substrate
 * (dsolve_run_parametric) verifies it by substituting y' = D[Y,p]/D[X,p] into the
 * ORIGINAL residual and assembles {{x -> Function[{p}, X], y -> Function[{p}, Y]}}
 * -- so a wrong solve/recursion branch is dropped (0 FAIL by construction).
 *
 * This GENERALIZES DSolve`Lagrange, which is the special case where y = x phi(p)
 * + psi(p) makes the induced x(p) ODE LINEAR (solved directly by the
 * integrating-factor helper rather than by recursing the cascade).  Degenerate
 * branches decline rather than mis-answer: p - G_x == 0 is the Clairaut/singular
 * family (owned earlier); G free of p is not this method; a Root[...] / implicit
 * solve branch has no closed form.
 *
 * Runs LATE in the first-order cascade -- after every deterministic specialist and
 * before the heuristic Lie backstop -- because it recurses DSolve and is therefore
 * expensive (matching Maple's late dp ordering).  All recursion is bounded exactly
 * as M14 (TimeConstrained sub-solves + a wall-clock deadline + a re-entry guard +
 * a per-top-level decline memo), since the induced-ODE solve re-enters DSolve.
 * Parametric IVP-fitting is future (an IVP is declined by the substrate), as with
 * DSolve`Lagrange.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

/* ---- recursion bounding kit (mirrors dsolve_changevar.c) ---- */
static time_t g_sfy_deadline;
static bool sfy_expired(void) { return time(NULL) >= g_sfy_deadline; }
static int  g_sfy_active = 0;   /* re-entry guard: no solve-for inside solve-for */

/* decline memo (the evaluator re-invokes a declining builtin ~3x/call) */
#define SFY_MEMO_SLOTS 64
static uint64_t sfy_epoch = 0;
static int sfy_memo_n = 0;
static uint64_t sfy_memo[SFY_MEMO_SLOTS];
static void sfy_memo_sync(uint64_t tid){ if(tid!=sfy_epoch){sfy_epoch=tid;sfy_memo_n=0;} }
static bool sfy_memo_seen(uint64_t h){ for(int i=0;i<sfy_memo_n;i++) if(sfy_memo[i]==h) return true; return false; }
static void sfy_memo_add(uint64_t h){ if(sfy_memo_n<SFY_MEMO_SLOTS && !sfy_memo_seen(h)) sfy_memo[sfy_memo_n++]=h; }

static void sfy_count(const Expr* e, long* a, long b){
    if(*a>b||!e) return;
    (*a)++;
    if(e->type==EXPR_FUNCTION){ sfy_count(e->data.function.head,a,b);
        for(size_t i=0;i<e->data.function.arg_count;i++) sfy_count(e->data.function.args[i],a,b); }
}
static bool sfy_too_big(const Expr* e){ long a=0; sfy_count(e,&a,40000); return a>40000; }

static Expr* sfy_beval(Expr* e, int secs){
    Expr* g = expr_new_function(expr_new_symbol("TimeConstrained"),
                  (Expr*[]){ e, expr_new_integer(secs), expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    return eval_and_free(g);
}
static bool sfy_aborted(const Expr* e){
    return e && e->type==EXPR_SYMBOL && e->data.symbol.name==intern_symbol("$Aborted");
}

/* Simplify `e` under a time bound; on abort / blow-up keep the raw form (still a
 * correct expression, just uglier -- the parametric substrate verifies either
 * way).  Unbounded Simplify on a messy recovered coordinate could spin for tens
 * of seconds (e.g. a 1 - p^3 induced-ODE integrating factor). `e` is consumed. */
static Expr* sfy_bsimplify(Expr* e, int secs){
    Expr* s = sfy_beval(expr_new_function(expr_new_symbol(SYM_Simplify),
                  (Expr*[]){ expr_copy(e) }, 1), secs);
    if (!s || sfy_aborted(s) || sfy_too_big(s)) { if (s) expr_free(s); return e; }
    expr_free(e);
    return s;
}

/* Solve DSolve[eqn, F[t], t] (TimeConstrained) -> applied-form RHS body or NULL. eqn consumed. */
static Expr* sfy_run_applied(Expr* eqn, const char* fname, const char* tvar, int secs){
    Expr* lhs = ds_call1(fname, expr_new_symbol(tvar));
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eqn, lhs, expr_new_symbol(tvar) }, 3);
    Expr* r = sfy_beval(call, secs);
    Expr* body = NULL;
    if (head_is(r, SYM_List) && r->data.function.arg_count >= 1){
        Expr* inner = r->data.function.args[0];
        if (head_is(inner, SYM_List))
            for (size_t k=0;k<inner->data.function.arg_count && !body;k++){
                Expr* rule = inner->data.function.args[k];
                if (head_is(rule, SYM_Rule) && rule->data.function.arg_count==2){
                    Expr* rl = rule->data.function.args[0];
                    if (rl->type==EXPR_FUNCTION && rl->data.function.head->type==EXPR_SYMBOL
                        && rl->data.function.head->data.symbol.name==fname)
                        body = expr_copy(rule->data.function.args[1]);
                }
            }
    }
    expr_free(r);
    return body;
}

/* base^-1; base consumed. */
static Expr* sfy_powneg1(Expr* base){
    return expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ base, expr_new_integer(-1) }, 2);
}

/* A clean bare parameter symbol (Function args must be bare symbols): the first of
 * a small pool not equal to the reserved vars and not occurring in X or Y. */
static const char* sfy_pick_param(const Expr* X, const Expr* Y, const char* xvar, const char* yname){
    static const char* pool[] = { "t","s","u","w","r","q","v" };
    for (size_t i=0;i<sizeof(pool)/sizeof(pool[0]);i++){
        const char* c = intern_symbol(pool[i]);
        if (c==xvar || c==yname) continue;
        if (ds_contains(X,c) || ds_contains(Y,c)) continue;
        return c;
    }
    return intern_symbol("t");
}

/* Shared core.  for_y: solve for y = G(x,p) and solve the induced x(p) ODE.
 * !for_y (SolvableForX): solve for x = H(y,p) and solve the induced y(p) ODE.
 * Returns a malloc'd array of DSolve`Param[X,Y,t] branch wrappers, or NULL. */
static Expr** sfy_solve_for(DSolveProblem* P, bool for_y, size_t* nbranch){
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    if (P->ncond > 0) return NULL;          /* parametric IVP-fitting is future */
    if (g_sfy_active) return NULL;          /* no solve-for inside solve-for */

    const char* xvar  = P->ind_names[0];
    const char* yname = P->fun_names[0];
    const char* Yn = intern_symbol("DSolve`Y");
    const char* Pn = intern_symbol("DSolve`p");

    /* decline memo, keyed per residual AND direction (a for-y decline must not
     * suppress the for-x attempt). */
    uint64_t h = expr_hash(P->eq_residuals[0]) ^ (for_y ? 0x9E3779B97F4A7C15ull : 0ull);
    sfy_memo_sync(eval_toplevel_id());
    if (sfy_memo_seen(h)) return NULL;

    /* R(x, Y, p): residual with y'[x] -> p, y[x] -> Y */
    Expr* R = dsolve_algebraic_residual(P, Yn, Pn);
    if (!R) return NULL;

    const char* solve_var = for_y ? Yn   : xvar;   /* isolate this variable      */
    const char* keep_var  = for_y ? xvar : Yn;     /* becomes the unknown F(p)   */

    /* Cheap pre-gate: only attempt when R is POLYNOMIAL in the isolated variable.
     * A transcendental dependence (Tan[x y], Log[Log[y]], ...) makes Solve emit
     * inverse-function branches and the induced ODE rarely closes, while burning
     * ~10 s per direction (a P->U timing regression).  A polynomial dependence is
     * exactly the tractable "solvable for y/x" class (the 347-family, the
     * f(y')-coefficient forms), and Solve on it is fast and clean.  Mirrors
     * NthAlgebraic's PolynomialQ gate. */
    Expr* pqE = eval_and_free(ds_call2("PolynomialQ", expr_copy(R), expr_new_symbol(solve_var)));
    bool ispoly = (pqE && pqE->type==EXPR_SYMBOL && pqE->data.symbol.name==SYM_True);
    if (pqE) expr_free(pqE);
    if (!ispoly) { expr_free(R); sfy_memo_add(h); return NULL; }

    g_sfy_deadline = time(NULL) + 4;

    /* Solve R==0 for the isolated variable (TimeConstrained backstop). */
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                   (Expr*[]){ R, expr_new_integer(0) }, 2);   /* consumes R */
    Expr* sol = sfy_beval(expr_new_function(expr_new_symbol(SYM_Solve),
                    (Expr*[]){ eq, expr_new_symbol(solve_var) }, 2), 3);   /* consumes eq */
    size_t nG = 0;
    Expr** Gs = dsolve_extract_solutions(sol, solve_var, &nG);
    if (sol) expr_free(sol);
    if (!Gs) { sfy_memo_add(h); return NULL; }

    const char* Ff = intern_symbol("DSolve`sfF");   /* the induced-ODE unknown F(p) */
    size_t nb = 0, cap = 4;
    Expr** out = malloc(cap * sizeof(Expr*));

    for (size_t i = 0; i < nG && !sfy_expired(); i++){
        Expr* G = Gs[i]; Gs[i] = NULL;   /* take ownership */
        /* need a genuine closed form that still depends on p and is not implicit */
        if (!G || sfy_too_big(G) || ds_has_head(G, intern_symbol("Root"))
            || ds_contains(G, solve_var)) { if (G) expr_free(G); continue; }

        Expr* Gp = ds_d(expr_copy(G), expr_new_symbol(Pn));       /* dG/dp        */
        if (ds_is_zero(Gp)) { expr_free(Gp); expr_free(G); continue; }   /* free of p */
        Expr* Gk = ds_d(expr_copy(G), expr_new_symbol(keep_var)); /* dG/d(keep)   */

        /* denominator:  (p - G_k)  [for_y]   or   (1/p - G_k)  [for_x] */
        Expr* denom = eval_and_free(for_y
            ? ds_call2(SYM_Subtract, expr_new_symbol(Pn), Gk)
            : ds_call2(SYM_Subtract, sfy_powneg1(expr_new_symbol(Pn)), Gk));
        if (ds_is_zero(denom)) { expr_free(denom); expr_free(Gp); expr_free(G); continue; }

        /* induced ODE:  F'[p] == G_p / denom, with keep_var -> F[p]. */
        Expr* rhs = eval_and_free(ds_call2(SYM_Times, Gp, sfy_powneg1(denom)));
        rhs = ds_subst(rhs, expr_new_symbol(keep_var), ds_make_funcapp(Ff, 0, Pn));
        Expr* ode = expr_new_function(expr_new_symbol(SYM_Equal),
                        (Expr*[]){ ds_make_funcapp(Ff, 1, Pn), rhs }, 2);   /* consumes rhs */

        /* recurse under the remaining wall-clock budget (hard-caps each direction) */
        long rem = (long)(g_sfy_deadline - time(NULL));
        if (rem < 1) { expr_free(ode); expr_free(G); break; }
        g_sfy_active++;
        Expr* Fsol = sfy_run_applied(ode, Ff, Pn, (int)(rem < 3 ? rem : 3));
        g_sfy_active--;
        if (!Fsol) { expr_free(G); continue; }

        /* other coordinate: G with keep_var -> F(p) (Simplify time-bounded) */
        Expr* other = sfy_bsimplify(ds_subst(expr_copy(G), expr_new_symbol(keep_var), expr_copy(Fsol)), 2);
        expr_free(G);

        Expr* X = for_y ? Fsol  : other;   /* for_y: x=F(p), y=G(F,p)             */
        Expr* Y = for_y ? other : Fsol;    /* for_x: y=F(p), x=H(F,p)             */
        if (sfy_too_big(X) || sfy_too_big(Y)) { expr_free(X); expr_free(Y); continue; }

        /* rename the parameter p -> a clean symbol t for the Function form */
        const char* tname = sfy_pick_param(X, Y, xvar, yname);
        X = ds_subst(X, expr_new_symbol(Pn), expr_new_symbol(tname));
        Y = ds_subst(Y, expr_new_symbol(Pn), expr_new_symbol(tname));

        if (nb >= cap) { cap *= 2; out = realloc(out, cap * sizeof(Expr*)); }
        out[nb++] = expr_new_function(expr_new_symbol("DSolve`Param"),
                        (Expr*[]){ X, Y, expr_new_symbol(tname) }, 3);   /* consumes X, Y */
    }
    for (size_t i = 0; i < nG; i++) if (Gs[i]) expr_free(Gs[i]);
    free(Gs);

    if (nb == 0) { free(out); sfy_memo_add(h); return NULL; }
    *nbranch = nb;
    return out;
}

Expr** dsolve_solvablefory_try(DSolveProblem* P, size_t* nbranch){ return sfy_solve_for(P, true,  nbranch); }
Expr** dsolve_solvableforx_try(DSolveProblem* P, size_t* nbranch){ return sfy_solve_for(P, false, nbranch); }

static Expr* builtin_dsolve_solvablefory(Expr* res){ return dsolve_method_builtin_parametric(res, dsolve_solvablefory_try); }
static Expr* builtin_dsolve_solvableforx(Expr* res){ return dsolve_method_builtin_parametric(res, dsolve_solvableforx_try); }

void dsolve_solvefor_init(void){
    symtab_add_builtin("DSolve`SolvableForY", builtin_dsolve_solvablefory);
    symtab_get_def("DSolve`SolvableForY")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`SolvableForY",
        "DSolve`SolvableForY[eqn, y, x] solves an ODE F(x,y,y')==0 that is algebraic "
        "in y by isolating y = G(x,y') and differentiating: the slope p = y' becomes "
        "the parameter of the first-order ODE dx/dp = G_p/(p - G_x), solved by "
        "recursion; the general solution is parametric {{x -> Function[{p}, X], "
        "y -> Function[{p}, Y]}}.  Generalizes DSolve`Lagrange (its linear case).");
    symtab_add_builtin("DSolve`SolvableForX", builtin_dsolve_solvableforx);
    symtab_get_def("DSolve`SolvableForX")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`SolvableForX",
        "DSolve`SolvableForX[eqn, y, x] solves an ODE F(x,y,y')==0 that is algebraic "
        "in x by isolating x = H(y,y') and differentiating w.r.t. y (dx/dy = 1/y'): "
        "the slope p = y' parameterizes the first-order ODE dy/dp = H_p/(1/p - H_y), "
        "solved by recursion; the general solution is parametric.");
}
