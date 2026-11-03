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
#include "dsolve.h"
#include "../core.h"
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

/* DSolve[eqn, F[t], t] under TimeConstrained -> the raw result.  eqn consumed. */
static Expr* sfy_run_raw(Expr* eqn, const char* fname, const char* tvar, int secs){
    Expr* lhs = ds_call1(fname, expr_new_symbol(tvar));
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eqn, lhs, expr_new_symbol(tvar) }, 3);
    return sfy_beval(call, secs);
}

/* The applied-form RHS body of the first explicit branch of a DSolve result for
 * F, or NULL.  r borrowed. */
static Expr* sfy_applied_body(const Expr* r, const char* fname){
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
    return body;
}

/* The left side Phi of a single-branch IMPLICIT result {{Phi == C[1]}} (Phi free
 * of the generated constant), or NULL.  r borrowed. */
static Expr* sfy_implicit_lhs(const Expr* r){
    if (!head_is(r, SYM_List) || r->data.function.arg_count != 1) return NULL;
    const Expr* b = r->data.function.args[0];
    if (!head_is(b, SYM_List) || b->data.function.arg_count != 1) return NULL;
    const Expr* eq = b->data.function.args[0];
    if (!head_is(eq, SYM_Equal) || eq->data.function.arg_count != 2) return NULL;
    Expr* Ck = ds_const(1);
    bool ok = expr_eq(eq->data.function.args[1], Ck)
              && !ds_contains(eq->data.function.args[0], intern_symbol("C"));
    expr_free(Ck);
    return ok ? expr_copy(eq->data.function.args[0]) : NULL;
}

/* ---- implicit elimination (DSolve`SolvableForYImplicit) ----
 * When the induced x(p) ODE solves only IMPLICITLY, Phi(p, x) == C[1] (e.g.
 * y'^2 == x^2 + y: y = p^2 - x^2 gives the homogeneous dx/dp = 2p/(p + 2x)), the
 * parametric form is unavailable, but the slope can be ELIMINATED instead: when
 * F(x, y, p) == 0 solves explicitly for p, each root p_k(x, y) turns the relation
 * into an implicit solution Phi(p_k(x, y[x]), x) == C[1] of the original ODE --
 * Maple's dp-method answer.  Taken only with a single y-root (so every p-root lies
 * on that branch) and only when every p-root converts, so the family is complete.
 * The relations are stashed for the implicit pass placed right after the parametric
 * slot in the cascade (no second solve), which re-verifies each against the
 * ORIGINAL ODE by the implicit-function rule. */
#define SFY_STASH_TAG 2

/* Eliminate p from Phi(p, F[p]) == C[1]: F[p] -> x, p -> each explicit root of
 * R(x, Y, p) == 0, Y -> y[x].  Returns the relations (all or nothing).  Phi and
 * R borrowed. */
static Expr** sfy_eliminate_slope(const Expr* Phi, const Expr* R, const char* Ff,
                                  const char* Pn, const char* Yn, const char* xvar,
                                  const char* yname, size_t* nout){
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                   (Expr*[]){ expr_copy((Expr*)R), expr_new_integer(0) }, 2);
    Expr* sol = sfy_beval(expr_new_function(expr_new_symbol(SYM_Solve),
                    (Expr*[]){ eq, expr_new_symbol(Pn) }, 2), 3);
    size_t np = 0;
    Expr** ps = dsolve_extract_solutions(sol, Pn, &np);
    if (sol) expr_free(sol);
    if (!ps) return NULL;
    Expr* base = ds_subst(expr_copy((Expr*)Phi), ds_make_funcapp(Ff, 0, Pn), expr_new_symbol(xvar));
    bool ok = !ds_contains(base, Ff);
    Expr** out = malloc((np ? np : 1) * sizeof(Expr*));
    size_t no = 0;
    for (size_t k = 0; k < np; k++){
        Expr* pk = ps[k]; ps[k] = NULL;
        if (!ok || !pk || ds_has_head(pk, intern_symbol("Root")) || ds_contains(pk, Pn)){
            if (pk) expr_free(pk);
            ok = false; continue;
        }
        Expr* G = ds_subst(expr_copy(base), expr_new_symbol(Pn), pk);
        G = ds_subst(G, expr_new_symbol(Yn), ds_make_funcapp(yname, 0, xvar));
        if (sfy_too_big(G) || !ds_contains(G, yname)){ expr_free(G); ok = false; continue; }
        out[no++] = G;
    }
    free(ps);
    expr_free(base);
    if (!ok || no == 0){ for (size_t k=0;k<no;k++) expr_free(out[k]); free(out); return NULL; }
    *nout = no;
    return out;
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
/* The automatic SolvableForX gate (lin_gate): the induced ODE dF/dp = rhs(p, F)
 * must be LINEAR or BERNOULLI in F -- rhs = (sum of c_k(p) F^k)/d(p) with the
 * occurring degrees {0, 1}, {1}, or {1, m} -- which the cascade closes in
 * milliseconds and whose parametric answer is rational/elementary in p.  rhs
 * borrowed, in terms of keep (the bare unknown symbol). */
static bool sfy_linear_induced(const Expr* rhs, const char* keep){
    Expr* tg  = eval_and_free(ds_call1(SYM_Together, expr_copy((Expr*)rhs)));
    Expr* num = eval_and_free(ds_call1("Numerator", expr_copy(tg)));
    Expr* den = eval_and_free(ds_call1("Denominator", tg));
    bool ok = ds_free_of(den, keep);
    expr_free(den);
    if (ok) {
        Expr* pq = eval_and_free(ds_call2("PolynomialQ", expr_copy(num), expr_new_symbol(keep)));
        ok = (pq && pq->type == EXPR_SYMBOL && pq->data.symbol.name == SYM_True);
        if (pq) expr_free(pq);
    }
    if (ok) {
        Expr* cl = eval_and_free(ds_call2("CoefficientList", expr_copy(num), expr_new_symbol(keep)));
        ok = head_is(cl, SYM_List) && cl->data.function.arg_count >= 2
             && cl->data.function.arg_count <= 16;
        size_t n_hi = 0;
        bool has1 = false;
        for (size_t k = 0; ok && k < cl->data.function.arg_count; k++) {
            if (ds_is_zero(cl->data.function.args[k])) continue;
            if (k == 1) has1 = true;
            else if (k >= 2) n_hi++;
        }
        /* linear: {0, 1} or {1};  Bernoulli: {1, m} (no constant term) */
        if (ok) ok = has1 && n_hi <= 1 && !(n_hi == 1 && !ds_is_zero(cl->data.function.args[0]));
        expr_free(cl);
    }
    expr_free(num);
    return ok;
}

static Expr** sfy_solve_for(DSolveProblem* P, bool for_y, bool lin_gate, size_t* nbranch){
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    if (P->ncond > 0) return NULL;          /* parametric IVP-fitting is future */
    if (g_sfy_active) return NULL;          /* no solve-for inside solve-for */
    /* Not below the first nesting level.  Each attempt recurses the whole engine
     * on the induced ODE under its own 4 s window, so a solve-for two levels deep
     * spends the budget of every method above it: in ExactODE's reduction of the
     * m62 stress equation 2 y y(3) + 2 (y + 3 y') y(2) + 2 (y')^2 == 0 it burned
     * 2.95 s at depth 3 before Lie (0.15 s) was reached, and the outer 3 s window
     * then expired -- HEAD had passed only through an accidental full re-solve on
     * a warm cache.  Declined there, the equation solves cold in 0.7 s. */
    if (g_dsolve_depth > 2) return NULL;

    const char* xvar  = P->ind_names[0];
    const char* yname = P->fun_names[0];
    const char* Yn = intern_symbol("DSolve`Y");
    const char* Pn = intern_symbol("DSolve`p");

    /* decline memo, keyed per residual AND direction (a for-y decline must not
     * suppress the for-x attempt). */
    uint64_t h = expr_hash(P->eq_residuals[0]) ^ (for_y ? 0x9E3779B97F4A7C15ull : 0ull)
                 ^ (lin_gate ? 0xD1B54A32D192ED03ull : 0ull);
    sfy_memo_sync(eval_toplevel_id());
    if (sfy_memo_seen(h)) return NULL;

    /* R(x, Y, p): residual with y'[x] -> p, y[x] -> Y */
    Expr* R = dsolve_algebraic_residual(P, Yn, Pn);
    if (!R) return NULL;
    Expr* R0 = expr_copy(R);          /* kept for the implicit slope elimination */

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
    if (!ispoly) { expr_free(R); expr_free(R0); sfy_memo_add(h); return NULL; }

    g_sfy_deadline = time(NULL) + 4;

    /* Solve R==0 for the isolated variable (TimeConstrained backstop). */
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                   (Expr*[]){ R, expr_new_integer(0) }, 2);   /* consumes R */
    Expr* sol = sfy_beval(expr_new_function(expr_new_symbol(SYM_Solve),
                    (Expr*[]){ eq, expr_new_symbol(solve_var) }, 2), 3);   /* consumes eq */
    size_t nG = 0;
    Expr** Gs = dsolve_extract_solutions(sol, solve_var, &nG);
    if (sol) expr_free(sol);
    if (!Gs) { expr_free(R0); sfy_memo_add(h); return NULL; }
    if (lin_gate && nG != 1) {          /* the automatic gate: one isolated root only */
        for (size_t i = 0; i < nG; i++) if (Gs[i]) expr_free(Gs[i]);
        free(Gs); expr_free(R0); sfy_memo_add(h); return NULL;
    }

    Expr** imp = NULL; size_t nimp = 0;   /* implicit relations (single y-root only) */
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
        if (lin_gate && !sfy_linear_induced(rhs, keep_var)) { expr_free(rhs); expr_free(G); continue; }
        rhs = ds_subst(rhs, expr_new_symbol(keep_var), ds_make_funcapp(Ff, 0, Pn));
        Expr* ode = expr_new_function(expr_new_symbol(SYM_Equal),
                        (Expr*[]){ ds_make_funcapp(Ff, 1, Pn), rhs }, 2);   /* consumes rhs */

        /* recurse under the remaining wall-clock budget (hard-caps each direction) */
        long rem = (long)(g_sfy_deadline - time(NULL));
        if (rem < 1) { expr_free(ode); expr_free(G); break; }
        g_sfy_active++;
        Expr* Fr = sfy_run_raw(ode, Ff, Pn, (int)(rem < 3 ? rem : 3));
        g_sfy_active--;
        Expr* Fsol = sfy_applied_body(Fr, Ff);
        if (!Fsol && for_y && nG == 1 && !getenv("MATHILDA_NO_SFY_IMPLICIT")) {
            Expr* Phi = sfy_implicit_lhs(Fr);
            if (Phi) {
                imp = sfy_eliminate_slope(Phi, R0, Ff, Pn, Yn, xvar, yname, &nimp);
                expr_free(Phi);
            }
        }
        expr_free(Fr);
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
    expr_free(R0);

    /* stash the eliminated relations only when the parametric pass found nothing */
    if (imp) {
        if (nb == 0) dsolve_stash_put(SFY_STASH_TAG, h, imp, nimp);
        else { for (size_t k=0;k<nimp;k++) expr_free(imp[k]); free(imp); }
    }

    if (nb == 0) { free(out); sfy_memo_add(h); return NULL; }
    *nbranch = nb;
    return out;
}

Expr** dsolve_solvablefory_try(DSolveProblem* P, size_t* nbranch){ return sfy_solve_for(P, true,  false, nbranch); }
Expr** dsolve_solvableforx_try(DSolveProblem* P, size_t* nbranch){ return sfy_solve_for(P, false, false, nbranch); }
/* The automatic SolvableForX: one root x = H(y, p) whose induced y(p) ODE is linear
 * or Bernoulli (sfy_linear_induced) -- the x-analogue of Lagrange, e.g.
 * y^2 y'^3 - x y' + y == 0: x = y^2 p^2 + y/p, dy/dp = -y/p + 1/(2 p^4).  The
 * ungated method stays pinned-only (its general induced ODE can be costly to
 * verify); MATHILDA_NO_SFX_LINEAR disables this slot for A/B. */
Expr** dsolve_solvableforx_linear_try(DSolveProblem* P, size_t* nbranch){
    if (getenv("MATHILDA_NO_SFX_LINEAR")) return NULL;
    if (g_dsolve_depth > 1) return NULL;   /* top level: a nested caller needs y(x) */
    return sfy_solve_for(P, false, true, nbranch);
}

/* The implicit pass: consume what the parametric for-y pass stashed for this
 * residual (pinned on its own, run the parametric pass first to populate it). */
static Expr** sfy_implicit_impl(DSolveProblem* P, size_t* nbranch, bool may_solve){
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 1) return NULL;
    uint64_t h = expr_hash(P->eq_residuals[0]) ^ 0x9E3779B97F4A7C15ull;
    size_t n = 0;
    Expr** G = dsolve_stash_take(SFY_STASH_TAG, h, &n);
    if (!G && may_solve){
        size_t ne = 0;
        Expr** ex = sfy_solve_for(P, true, false, &ne);
        if (ex){ for (size_t i=0;i<ne;i++) expr_free(ex[i]); free(ex); return NULL; }
        G = dsolve_stash_take(SFY_STASH_TAG, h, &n);
    }
    if (!G) return NULL;
    *nbranch = n;
    return G;
}
Expr** dsolve_solvablefory_implicit_try(DSolveProblem* P, size_t* nbranch){ return sfy_implicit_impl(P, nbranch, false); }
static Expr** sfy_implicit_pinned_try(DSolveProblem* P, size_t* nbranch){ return sfy_implicit_impl(P, nbranch, true); }

static Expr* builtin_dsolve_solvablefory(Expr* res){ return dsolve_method_builtin_parametric(res, dsolve_solvablefory_try); }
static Expr* builtin_dsolve_solvableforx(Expr* res){ return dsolve_method_builtin_parametric(res, dsolve_solvableforx_try); }
static Expr* builtin_dsolve_solvablefory_implicit(Expr* res){ return dsolve_method_builtin_implicit(res, sfy_implicit_pinned_try); }

void dsolve_solvefor_init(void){
    tc_register_unwind_counter(&g_sfy_active);   /* restored on a timeout unwind */
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
    symtab_add_builtin("DSolve`SolvableForYImplicit", builtin_dsolve_solvablefory_implicit);
    symtab_get_def("DSolve`SolvableForYImplicit")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`SolvableForYImplicit",
        "DSolve`SolvableForYImplicit[eqn, y, x] is DSolve`SolvableForY when the induced "
        "x(p) ODE solves only implicitly, Phi(p, x) == C: the slope p = y' is eliminated "
        "through each explicit root p_k(x, y) of eqn, giving one verified implicit "
        "solution Phi(p_k(x, y[x]), x) == C[1] per root.");
}
