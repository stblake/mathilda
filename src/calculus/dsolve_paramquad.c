/*
 * dsolve_paramquad.c — DSolve`ParametricQuadrature.
 *
 * A first-order ODE F(x, y, p) == 0 (p = y') from which one of x, y is MISSING is
 * a "quadrature" in the parametric sense, even when it cannot be solved for p:
 *
 *   y-free  F(x, p) == 0.  If the curve F(x, p) = 0 has a parametrisation
 *           x = U(t), p = P(t), then dy = p dx gives
 *                 x = U(t),   y = Integrate[P(t) U'(t), t] + C[1].
 *   x-free  F(y, p) == 0.  With y = U(t), p = P(t), dx = dy/p gives
 *                 y = U(t),   x = Integrate[U'(t)/P(t), t] + C[1].
 *
 * The parametrisation is searched for by isolating the present variable u (x or y)
 * after one of four substitutions for the slope (Kamke §I.4 / Maple `_quadrature`):
 *
 *     p = t          u = G(t) when F is solvable for u      (Sin[p] + p == x)
 *     p = t u        the homogeneous case                   (p^3 - a x p + x^3 == 0)
 *     p = t + u      a cubic in (p - u)                     ((p - x)^3 + p (p - x) == 1)
 *     p = t - u                                             ((p + y)^4 == p + 3 y)
 *
 * Each candidate root u = G(t) must be a closed form (no Root / inverse function /
 * Solve-generated constant) depending on t, and is screened NUMERICALLY with every
 * parameter sampled before it is used: R(U(t), P(t)) must vanish at sampled t, which
 * rejects a spurious radical root that Solve keeps.  The answer goes out through the
 * parametric substrate (dsolve_run_parametric), which re-verifies each branch by
 * substituting y' = Y'(t)/X'(t) into the ORIGINAL residual.
 *
 * Not recursive (Solve and Integrate only, each TimeConstrained), so no re-entry
 * guard is needed; a per-top-level decline memo stops the evaluator's re-invocation
 * from repeating a decline.  Runs late in the first-order cascade (after every
 * specialist that solves for p), MATHILDA_NO_PARAMQUAD disables it for A/B.
 */
#include "dsolve_common.h"
#include "dsolve.h"
#include "integrate.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../message.h"
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

bool dsolve_nth_algebraic_has_relations(const DSolveProblem* P);

#define PQ_SOLVE_SEC 2
#define PQ_INTEGRATE_SEC 2
#define PQ_MAX_BRANCHES 8

/* decline memo (the evaluator re-invokes a declining builtin) */
#define PQ_MEMO_SLOTS 32
static uint64_t pq_epoch = 0;
static int pq_memo_n = 0;
static uint64_t pq_memo[PQ_MEMO_SLOTS];
static void pq_memo_sync(uint64_t tid){ if(tid!=pq_epoch){pq_epoch=tid;pq_memo_n=0;} }
static bool pq_memo_seen(uint64_t h){ for(int i=0;i<pq_memo_n;i++) if(pq_memo[i]==h) return true; return false; }
static void pq_memo_add(uint64_t h){ if(pq_memo_n<PQ_MEMO_SLOTS && !pq_memo_seen(h)) pq_memo[pq_memo_n++]=h; }

/* f evaluated under TimeConstrained[_, secs, $Aborted], with diagnostics muted and
 * rolled back (a speculative Solve::ifun / Integrate::nonelem is probe noise).
 * Consumes f; NULL on abort. */
static Expr* pq_beval(Expr* f, int secs){
    Expr* g = expr_new_function(expr_new_symbol(SYM_TimeConstrained),
                  (Expr*[]){ f, expr_new_integer(secs), expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    unsigned long fired0 = mth_msg_fired_count();
    mth_msg_suppress_push();
    g_integrate_quiet++;
    Expr* r = eval_and_free(g);
    g_integrate_quiet--;
    mth_msg_suppress_pop();
    mth_msg_fired_restore(fired0);
    if (r && r->type == EXPR_SYMBOL && r->data.symbol.name == intern_symbol("$Aborted")) {
        expr_free(r); return NULL;
    }
    return r;
}

/* A closed-form candidate: no Root / inverse function / generated constant /
 * conditional / unevaluated Integrate. */
static bool pq_closed(const Expr* e){
    static const char* bad[] = { "Root", "InverseFunction", "ConditionalExpression",
                                 "C", "Integrate", "Solve", "$Aborted" };
    for (size_t i = 0; i < sizeof(bad)/sizeof(bad[0]); i++)
        if (ds_contains(e, intern_symbol(bad[i]))) return false;
    return true;
}

/* Distinct free parameter symbols of e (argument positions, not t, not protected). */
static void pq_params(const Expr* e, const char* t, const char** names, size_t* n, size_t cap){
    if (!e || *n >= cap) return;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        if (nm == t) return;
        SymbolDef* d = symtab_lookup(nm);
        if (d && (d->attributes & ATTR_PROTECTED)) return;
        for (size_t i = 0; i < *n; i++) if (names[i] == nm) return;
        names[(*n)++] = nm;
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        pq_params(e->data.function.args[i], t, names, n, cap);
}

/* Numeric screen of a residual r(t) that must vanish identically: every parameter
 * sampled at a distinct generic value and t swept over several reals.  True only
 * when at least two samples numericize and a majority of them are ~0. */
static bool pq_num_zero(const Expr* r, const char* t){
    static const double ts[] = { 0.37, 0.83, 1.29, 1.71, 2.43, 3.17 };
    const char* params[12]; size_t np = 0;
    pq_params(r, t, params, &np, 12);
    if (np >= 12) return false;
    Expr* base = expr_copy((Expr*)r);
    for (size_t k = 0; k < np; k++)
        base = ds_subst(base, expr_new_symbol(params[k]),
                        eval_and_free(expr_new_function(expr_new_symbol(SYM_Rational),
                            (Expr*[]){ expr_new_integer(113 + 37 * (int64_t)k), expr_new_integer(100) }, 2)));
    int nnum = 0, nsmall = 0;
    for (size_t i = 0; i < sizeof(ts)/sizeof(ts[0]); i++) {
        Expr* ev = ds_subst(expr_copy(base), expr_new_symbol(t), expr_new_real(ts[i]));
        ev = pq_beval(ds_call1("Abs", ds_call2("N", ev, expr_new_integer(30))), 1);
        if (!ev) continue;
        Expr* m = eval_and_free(ds_call1("N", ev));   /* to machine precision */
        double v = (m && m->type == EXPR_REAL) ? m->data.real
                 : (m && m->type == EXPR_INTEGER) ? (double)m->data.integer : NAN;
        if (m) expr_free(m);
        if (isnan(v) || !isfinite(v)) continue;
        nnum++;
        if (v < 1e-8) nsmall++;
    }
    expr_free(base);
    return nnum >= 2 && 2 * nsmall > nnum;
}

/* A clean bare name for the parameter of the final Function[{s}, ..] pair: the
 * first of a small pool occurring in neither coordinate (the computation itself
 * runs on the private DSolve`pqT).  A name carrying a value is skipped: the
 * parametric verifier differentiates with respect to it. */
static const char* pq_pick_param(const Expr* X, const Expr* Y){
    static const char* pool[] = { "t","s","w","r","q","v" };
    for (size_t i = 0; i < sizeof(pool)/sizeof(pool[0]); i++){
        const char* c = intern_symbol(pool[i]);
        if (ds_contains(X, c) || ds_contains(Y, c)) continue;
        SymbolDef* d = symtab_lookup(c);   /* the verifier evaluates with it: no value */
        if (d && (d->own_values || (d->attributes & ATTR_PROTECTED))) continue;
        return c;
    }
    return intern_symbol("DSolve`pqT");
}

/* The slope as a function of (t, u) for substitution kind k:
 *   0: t   1: t u   2: t + u   3: t - u */
static Expr* pq_slope(int k, const char* t, const char* u){
    switch (k) {
        case 0:  return expr_new_symbol(t);
        case 1:  return ds_call2(SYM_Times, expr_new_symbol(t), expr_new_symbol(u));
        case 2:  return ds_call2(SYM_Plus, expr_new_symbol(t), expr_new_symbol(u));
        default: return ds_call2(SYM_Subtract, expr_new_symbol(t), expr_new_symbol(u));
    }
}

/* Try one substitution kind.  Returns the number of DSolve`Param branches appended
 * to out[] (capacity PQ_MAX_BRANCHES), and sets *all when every candidate root of
 * the isolation produced a branch.  R borrowed. */
static size_t pq_try_kind(const DSolveProblem* P, const Expr* R, int kind, bool yfree, const char* u,
                          const char* Pn, const char* t, Expr** out, bool* all){
    *all = false;
    Expr* Rs = eval_and_free(ds_subst(expr_copy((Expr*)R), expr_new_symbol(Pn), pq_slope(kind, t, u)));
    if (!ds_contains(Rs, u) || !ds_contains(Rs, t)) { expr_free(Rs); return 0; }
    /* isolate u: the substituted residual (or its numerator) must be algebraic in u */
    Expr* num = eval_and_free(ds_call1("Numerator", eval_and_free(ds_call1(SYM_Together, expr_copy(Rs)))));
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                   (Expr*[]){ num, expr_new_integer(0) }, 2);
    Expr* sol = pq_beval(ds_call2(SYM_Solve, eq, expr_new_symbol(u)), PQ_SOLVE_SEC);
    size_t ng = 0;
    Expr** gs = sol ? dsolve_extract_solutions(sol, u, &ng) : NULL;
    if (sol) expr_free(sol);
    if (!gs) { expr_free(Rs); return 0; }

    size_t nb = 0, ncand = 0;
    for (size_t i = 0; i < ng; i++) {
        Expr* G = gs[i]; gs[i] = NULL;
        if (!G) continue;
        if (!ds_contains(G, t) || ds_contains(G, u) || !pq_closed(G)) { expr_free(G); continue; }
        ncand++;
        if (nb >= PQ_MAX_BRANCHES) { expr_free(G); continue; }
        /* the slope along the root, and the residual screen R(G(t), P(t)) == 0 */
        Expr* Pt = eval_and_free(ds_subst(pq_slope(kind, t, u), expr_new_symbol(u), expr_copy(G)));
        Expr* chk = ds_subst(ds_subst(expr_copy((Expr*)R), expr_new_symbol(Pn), expr_copy(Pt)),
                             expr_new_symbol(u), expr_copy(G));
        bool ok = !ds_is_zero(Pt) && pq_num_zero(chk, t);
        expr_free(chk);
        Expr* Gp = ok ? ds_d(expr_copy(G), expr_new_symbol(t)) : NULL;
        if (ok && ds_is_zero(Gp)) ok = false;
        if (!ok) { expr_free(Pt); expr_free(G); if (Gp) expr_free(Gp); continue; }

        /* the other coordinate by quadrature */
        Expr* integrand = yfree
            ? ds_call2(SYM_Times, expr_copy(Pt), Gp)
            : ds_call2(SYM_Times, Gp, ds_call2(SYM_Power, expr_copy(Pt), expr_new_integer(-1)));
        integrand = eval_and_free(ds_call1(SYM_Together, integrand));
        Expr* I = pq_beval(ds_call2(SYM_Integrate, integrand, expr_new_symbol(t)), PQ_INTEGRATE_SEC);
        expr_free(Pt);
        if (!I || !pq_closed(I)) { if (I) expr_free(I); expr_free(G); continue; }
        Expr* other = eval_and_free(ds_call2(SYM_Plus, I, ds_const(1)));
        Expr* X = yfree ? G : other;
        Expr* Y = yfree ? other : G;
        /* rename the internal parameter to a clean bare symbol for Function[{s}, ..] */
        const char* s = pq_pick_param(X, Y);
        X = ds_subst(X, expr_new_symbol(t), expr_new_symbol(s));
        Y = ds_subst(Y, expr_new_symbol(t), expr_new_symbol(s));
        /* the substrate's verify, here: an antiderivative valid on part of the
         * parameter line only (Sqrt[t^2 - t^3] read as t Sqrt[1 - t]) is dropped,
         * so the next parametrisation gets its turn */
        if (!dsolve_verify_parametric(P, X, Y, s)) { expr_free(X); expr_free(Y); continue; }
        out[nb++] = expr_new_function(expr_new_symbol("DSolve`Param"),
                        (Expr*[]){ X, Y, expr_new_symbol(s) }, 3);
    }
    for (size_t i = 0; i < ng; i++) if (gs[i]) expr_free(gs[i]);
    free(gs);
    expr_free(Rs);
    *all = (nb > 0 && nb == ncand);
    return nb;
}

Expr** dsolve_paramquad_try(DSolveProblem* P, size_t* nbranch){
    if (getenv("MATHILDA_NO_PARAMQUAD")) return NULL;      /* A/B killswitch */
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 1) return NULL;
    if (P->ncond > 0) return NULL;          /* parametric IVP-fitting is future */
    /* Top level only: a nested caller (HomogeneousG's autonomous u(t) solve, a
     * reduction's sub-ODE) needs y as a function or a first integral, and a
     * parametric pair here would pre-empt the later slots that give it one
     * (§2.2.56-5536 lost its implicit answer that way). */
    if (g_dsolve_depth > 1) return NULL;
    /* NthAlgebraic found a first integral for every y'-root: its implicit slot
     * (later in the cascade) answers in x and y directly; leave it that */
    if (dsolve_nth_algebraic_has_relations(P)) return NULL;

    const char* xvar  = P->ind_names[0];
    const char* yname = P->fun_names[0];
    uint64_t h = expr_hash(P->eq_residuals[0]) ^ 0x5851F42D4C957F2Dull;
    pq_memo_sync(eval_toplevel_id());
    if (pq_memo_seen(h)) return NULL;

    const char* Yn = intern_symbol("DSolve`pqY");
    const char* Pn = intern_symbol("DSolve`pqP");
    Expr* R = dsolve_algebraic_residual(P, Yn, Pn);
    if (!R) return NULL;
    /* only y, y' and x may occur: no other derivative, no y at another point */
    bool hasY = ds_contains(R, Yn), hasX = ds_contains(R, xvar);
    if (ds_contains(R, yname) || !ds_contains(R, Pn) || hasY == hasX) { expr_free(R); return NULL; }
    bool yfree = !hasY;
    const char* u = yfree ? xvar : Yn;
    const char* t = intern_symbol("DSolve`pqT");

    Expr** best = NULL; size_t nbest = 0;
    for (int kind = 0; kind < 4; kind++) {
        Expr** out = malloc(PQ_MAX_BRANCHES * sizeof(Expr*));
        bool all = false;
        size_t nb = pq_try_kind(P, R, kind, yfree, u, Pn, t, out, &all);
        if (nb > 0 && (all || !best)) {
            if (best) { for (size_t i = 0; i < nbest; i++) expr_free(best[i]); free(best); }
            best = out; nbest = nb;
            if (all) break;                 /* a complete parametrisation */
        } else {
            for (size_t i = 0; i < nb; i++) expr_free(out[i]);
            free(out);
        }
    }
    expr_free(R);
    if (!best) { pq_memo_add(h); return NULL; }   /* a decline is not retried */
    *nbranch = nbest;
    return best;
}

static Expr* builtin_dsolve_paramquad(Expr* res){
    return dsolve_method_builtin_parametric(res, dsolve_paramquad_try);
}

void dsolve_paramquad_init(void){
    symtab_add_builtin("DSolve`ParametricQuadrature", builtin_dsolve_paramquad);
    symtab_get_def("DSolve`ParametricQuadrature")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`ParametricQuadrature",
        "DSolve`ParametricQuadrature[eqn, y, x] solves a first-order ODE from which x "
        "or y is missing, F(x, y') == 0 or F(y, y') == 0, through a parametrisation "
        "of the curve F == 0 by p = y' = t (or t u, t + u, t - u): x = U(t), "
        "y = Integrate[P(t) U'(t), t] + C[1] (y-free), or y = U(t), x = Integrate"
        "[U'(t)/P(t), t] + C[1] (x-free).  The answer is parametric "
        "{{x -> Function[{t}, X], y -> Function[{t}, Y]}}.");
}
