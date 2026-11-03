/*
 * dsolve_nth_algebraic.c — DSolve`NthAlgebraic.
 *
 * Handles an ODE that is *algebraic* (polynomial of degree >= 2) in its highest
 * derivative y^(n): Solve[eqn, y^(n)] gives several branches y^(n) == g_k, each
 * recursed back into the scalar cascade (a branch free of y hits Quadrature) and
 * the branch bodies unioned.  Also the degenerate no-derivative case (an equation
 * purely algebraic in y[x]): Solve[eqn, y[x]] gives the constant/algebraic
 * solutions directly.  Mirrors SymPy's `nth_algebraic`.
 *
 * Runs at the FRONT of the cascade (with Factorable), so a top-derivative-in-a-power
 * form is split before the linear specialists try to match it.  Declines the normal
 * case (degree 1 in the top derivative with derivatives present) so those methods
 * claim it.  The substrate back-substitution verifies every unioned branch.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../message.h"
#include "../common.h"
#include "integrate.h"
#include <stdlib.h>
#include <time.h>

/* This method sits SECOND in the scalar cascade and its per-root recursion is a full
 * `DSolve`, so an unbounded one lets any degree->=2-in-y' equation burn arbitrary time
 * before the specialists downstream are ever reached.  Measured on
 * (y - x y')^2 == 1 + y'^2: the two roots each cost ~5.9 s and each returns an
 * IMPLICIT relation, which `dsolve_extract_applied_bodies` (explicit bodies only) then
 * discards -- 11.5 s of work, all of it thrown away -- and, with no decline memo, the
 * evaluator's fixed-point re-invocation repeated it, for ~25 s in total on an equation
 * `Clairaut` answers in 0.01 s.  Hence the standard kit: a wall-clock deadline, a
 * per-branch TimeConstrained, and a per-top-level decline memo. */
#define NA_TOTAL_BUDGET_SEC   6
#define NA_BRANCH_BUDGET_SEC  4
static time_t g_na_deadline;
static bool na_expired(void) { return time(NULL) >= g_na_deadline; }

#define NA_MEMO_SLOTS 32
static uint64_t na_epoch = 0;
static int na_memo_n = 0;
static uint64_t na_memo[NA_MEMO_SLOTS];
static void na_memo_sync(uint64_t tid){ if(tid!=na_epoch){na_epoch=tid;na_memo_n=0;} }
static bool na_memo_seen(uint64_t h){ for(int i=0;i<na_memo_n;i++) if(na_memo[i]==h) return true; return false; }
static void na_memo_add(uint64_t h){ if(na_memo_n<NA_MEMO_SLOTS && !na_memo_seen(h)) na_memo[na_memo_n++]=h; }

/* ---- implicit second pass (DSolve`NthAlgebraicImplicit) ----
 * When no root's sub-solve is EXPLICIT but every root's sub-solve answered (an
 * implicit first integral, e.g. y'^2 == (y-1) y^2 -> x/2 + ArcTan[Sqrt[y-1]] ==
 * C[1]), the explicit pass declines and the relations would be lost.  They are
 * stashed instead, as first-integral LHSs G(x, y[x]) (G == C[1] solves the
 * branch), keyed by the residual hash within the top-level evaluation, and the
 * implicit pass -- placed LATE in the cascade so every specialist still claims
 * the original equation first -- hands them to dsolve_run_implicit, which
 * re-verifies each G against the ORIGINAL ODE by the implicit-function rule.
 * The stash costs no second solve; the explicit pass's behaviour is unchanged.
 * Shared substrate: dsolve_result_relations / dsolve_stash_* (dsolve_common). */
#define NA_STASH_TAG 1

/* Recurse DSolve[subeqn, y[x], x] (applied form) and append every explicit body
 * to *acc (grown by realloc).  When `rel` is non-NULL the branch is also
 * collected as first-integral LHSs into *rel; *covered is cleared when it cannot
 * be (an unsolved root).  `subeqn` is consumed. */
static void recurse_collect(Expr* subeqn, const char* yname, const char* xvar,
                            Expr*** acc, size_t* nacc,
                            Expr*** rel, size_t* nrel, bool* covered) {
    if (na_expired()) { expr_free(subeqn); if (covered) *covered = false; return; }
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ subeqn, ds_make_funcapp(yname, 0, xvar),
                                expr_new_symbol(xvar) }, 3);
    long left = (long)(g_na_deadline - time(NULL));
    if (left > NA_BRANCH_BUDGET_SEC) left = NA_BRANCH_BUDGET_SEC;
    if (left < 1) left = 1;
    call = expr_new_function(expr_new_symbol(intern_symbol("TimeConstrained")),
               (Expr*[]){ call, expr_new_integer(left),
                          expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    Expr* r = eval_and_free(call);
    size_t nb = 0;
    Expr** bodies = dsolve_extract_applied_bodies(r, yname, &nb);
    if (rel && covered && *covered && !dsolve_result_relations(r, yname, xvar, rel, nrel))
        *covered = false;
    expr_free(r);
    if (bodies) {
        *acc = realloc(*acc, (*nacc + nb) * sizeof(Expr*));
        for (size_t i = 0; i < nb; i++) (*acc)[(*nacc)++] = bodies[i];
        free(bodies);
    }
}

/* ---- perfect-square extraction from a quadratic's root radical ----
 * Solve returns the roots of a y'-quadratic as (-b ± Sqrt[s^2 r])/(2a) with the
 * square s^2 left inside the radical (5432: Sqrt[4 Y^2 + 4 Y^2 Cot[x]^2]).  That
 * literal form hides the structure -- the linear / separable / homogeneous
 * specialists all decline it and the general methods spin -- while the rewritten
 * root -Y Cot[x] ± Y Sqrt[1 + Cot[x]^2] is a linear ODE.  Each Sqrt[R] becomes
 * sqrt(c2) * prod f^Floor[m/2] * Sqrt[c1 * prod f^Mod[m,2]] from
 * FactorSquareFreeList[R] (c = c2 c1 the numeric content, c2 an exact square).
 * s Sqrt[r] == ±Sqrt[s^2 r] pointwise, so for a QUADRATIC, rewriting both roots
 * the same way leaves the pair {root+, root-} unchanged as a set; the caller
 * applies it only at degree 2 and re-checks that every rewritten root still
 * annihilates the polynomial, keeping the raw roots otherwise.  The rewrite itself
 * is ds_sqrt_extract_squares (dsolve_common.c), shared with LowerDerivativeReduction. */
/* Rewrite the degree-2 roots gs[0..ng) in place when every rewritten root still
 * satisfies poly(Dn) == 0 identically (else leave all of them untouched). */
static void na_normalise_quadratic_roots(Expr** gs, size_t ng, const Expr* poly,
                                         const char* Dn) {
    Expr** nw = malloc(ng * sizeof(Expr*));
    bool ok = true, changed = false;
    for (size_t i = 0; i < ng; i++) {
        nw[i] = ds_sqrt_extract_squares(gs[i]);
        if (!expr_eq(nw[i], gs[i])) changed = true;
    }
    for (size_t i = 0; i < ng && ok && changed; i++) {
        Expr* chk = ds_subst(expr_copy((Expr*)poly), expr_new_symbol(Dn), expr_copy(nw[i]));
        chk = eval_and_free(ds_call1(SYM_Together, eval_and_free(ds_call1(SYM_Expand, chk))));
        if (!ds_is_zero(chk)) ok = false;
        expr_free(chk);
    }
    for (size_t i = 0; i < ng; i++) {
        if (ok && changed) { expr_free(gs[i]); gs[i] = nw[i]; }
        else expr_free(nw[i]);
    }
    free(nw);
}

/* ---- binomial separable form:  y'^n == A(x) B(y) ----
 * A first-order ODE c y'^n + Q(x, y) == 0 whose S = -Q/c factors into an x-part
 * and a y-part (5612: y'^3 + f[x] (y - a)^2 (y - b)^2 == 0; 5667: y'^4 ==
 * (y - a)^3 (y - b)^2) has the n branches y' = w_k A^(1/n) B^(1/n), w_k the n-th
 * roots of unity -- as a SET exactly the n roots of S, since (A B)^(1/n) and
 * A^(1/n) B^(1/n) differ by an n-th root of unity at every point.  Each branch is
 * separable, with the first integral
 *
 *     Integrate[B(y)^(-1/n), y] - w_k Integrate[A(x)^(1/n), x] == C[1].
 *
 * Solve's roots keep the product under ONE radical (expanded), which hides the
 * split, and the recursion then pays a full non-elementary separable attempt per
 * root (5612: 6.6 s to decline).  So when the y-integral does not close quickly the
 * n relations are built here directly -- each integral attempted ONCE, a
 * non-elementary or timed-out one kept inert (Inactive[Integrate]) -- and stashed
 * for the implicit pass, which re-verifies each against the ORIGINAL ODE.  When the
 * y-integral is elementary nothing changes: the ordinary root recursion owns it
 * (and may invert explicitly).  A, B are built from FactorList's factors, so
 * A B == S holds by construction.  MATHILDA_NO_NA_BINOMIAL disables it for A/B. */
#define NA_BINOM_INT_SEC 2
#define NA_BINOM_STASH_TAG 41

/* Integrate[f, v] once under a short budget, muted; an unevaluated, aborted or
 * otherwise unclosed result becomes the inert Inactive[Integrate][f, v].  f consumed. */
static Expr* na_integral_or_inert(Expr* f, const char* v) {
    Expr* call = ds_call2(SYM_Integrate, expr_copy(f), expr_new_symbol(v));
    Expr* g = expr_new_function(expr_new_symbol(intern_symbol("TimeConstrained")),
                  (Expr*[]){ call, expr_new_integer(NA_BINOM_INT_SEC),
                             expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    unsigned long fired0 = mth_msg_fired_count();
    mth_msg_suppress_push();
    g_integrate_quiet++;
    Expr* r = eval_and_free(g);
    g_integrate_quiet--;
    mth_msg_suppress_pop();
    mth_msg_fired_restore(fired0);
    if (r && !ds_has_head(r, SYM_Integrate) && !ds_contains(r, intern_symbol("$Aborted"))) {
        expr_free(f);
        return r;
    }
    if (r) expr_free(r);
    Expr* ih = expr_new_function(expr_new_symbol(SYM_Inactive),
                   (Expr*[]){ expr_new_symbol(SYM_Integrate) }, 1);
    return expr_new_function(ih, (Expr*[]){ f, expr_new_symbol(v) }, 2);
}

/* Returns true when R (in Dn = y', Yn = y) is binomial-separable with a y-integral
 * that does not close; then *rel / *nrel receive the n first-integral LHSs in
 * terms of y[x].  Rsub borrowed. */
static bool na_binomial_separable(const Expr* Rsub, const char* Dn, long deg,
                                  const char* yname, const char* xvar,
                                  Expr*** rel, size_t* nrel) {
    if (getenv("MATHILDA_NO_NA_BINOMIAL") || deg < 2 || deg > 12) return false;
    const char* Yn = intern_symbol("DSolve`naY");
    Expr* R = ds_subst(expr_copy((Expr*)Rsub), ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));
    if (ds_contains(R, yname)) { expr_free(R); return false; }
    Expr* cl = eval_and_free(ds_call2("CoefficientList", R, expr_new_symbol(Dn)));
    bool ok = head_is(cl, SYM_List) && cl->data.function.arg_count == (size_t)deg + 1;
    for (long k = 1; ok && k < deg; k++)
        if (!ds_is_zero(cl->data.function.args[k])) ok = false;
    if (ok && (ds_contains(cl->data.function.args[deg], Yn)
               || ds_contains(cl->data.function.args[deg], xvar)
               || ds_is_zero(cl->data.function.args[deg]))) ok = false;
    Expr* S = NULL;
    if (ok) {
        S = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                ds_call2(SYM_Times, expr_copy(cl->data.function.args[0]),
                    ds_call2(SYM_Power, expr_copy(cl->data.function.args[deg]), expr_new_integer(-1)))));
        if (!ds_contains(S, Yn)) ok = false;           /* a pure quadrature: not this form */
    }
    expr_free(cl);
    if (!ok) { if (S) expr_free(S); return false; }

    /* split S's factors into the x-part A and the y-part B, each to the power 1/n */
    Expr* fl = eval_and_free(ds_call1("FactorList", S));
    if (!head_is(fl, SYM_List)) { expr_free(fl); return false; }
    size_t nf = fl->data.function.arg_count;
    Expr** ax = malloc((nf + 1) * sizeof(Expr*)); size_t na = 0;
    Expr** by = malloc((nf + 1) * sizeof(Expr*)); size_t nbf = 0;
    for (size_t i = 0; i < nf && ok; i++) {
        Expr* pr = fl->data.function.args[i];
        if (!head_is(pr, SYM_List) || pr->data.function.arg_count != 2) { ok = false; break; }
        Expr* f = pr->data.function.args[0];
        bool hasY = ds_contains(f, Yn), hasX = ds_contains(f, xvar);
        if (hasY && hasX) { ok = false; break; }
        /* f^(m/n): the x-side to +1/n, the y-side to -1/n (the integrand 1/B^(1/n)) */
        Expr* ex = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                       expr_copy(pr->data.function.args[1]),
                       expr_new_function(expr_new_symbol(SYM_Rational),
                           (Expr*[]){ expr_new_integer(hasY ? -1 : 1), expr_new_integer(deg) }, 2) }, 2));
        Expr* pw = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy(f), ex }, 2);
        if (hasY) by[nbf++] = pw; else ax[na++] = pw;
    }
    expr_free(fl);
    if (!ok || nbf == 0) {
        for (size_t i = 0; i < na; i++) expr_free(ax[i]);
        for (size_t i = 0; i < nbf; i++) expr_free(by[i]);
        free(ax); free(by);
        return false;
    }
    Expr* invB = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), by, nbf));
    Expr* A = na ? eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), ax, na))
                 : expr_new_integer(1);
    free(ax); free(by);

    /* the y-integral first: elementary -> leave the equation to the root recursion */
    Expr* IY = na_integral_or_inert(invB, Yn);
    if (!ds_has_head(IY, SYM_Inactive)) {
        expr_free(IY); expr_free(A);
        return false;
    }
    Expr* IX = ds_free_of(A, xvar)
             ? eval_and_free(ds_call2(SYM_Times, A, expr_new_symbol(xvar)))
             : na_integral_or_inert(A, xvar);
    Expr* IYy = ds_subst(IY, expr_new_symbol(Yn), ds_make_funcapp(yname, 0, xvar));
    Expr** out = malloc((size_t)deg * sizeof(Expr*));
    for (long k = 0; k < deg; k++) {
        Expr* w = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){
                      expr_new_integer(-1),
                      expr_new_function(expr_new_symbol(SYM_Rational),
                          (Expr*[]){ expr_new_integer(2 * k), expr_new_integer(deg) }, 2) }, 2));
        out[k] = eval_and_free(ds_call2(SYM_Subtract, expr_copy(IYy),
                     ds_call2(SYM_Times, w, expr_copy(IX))));
    }
    expr_free(IYy); expr_free(IX);
    *rel = out; *nrel = (size_t)deg;
    return true;
}

Expr** dsolve_nth_algebraic_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];
    int n = P->max_order[0];
    Expr* R = P->eq_residuals[0];

    uint64_t memo_h = expr_hash(R);
    na_memo_sync(eval_toplevel_id());
    if (na_memo_seen(memo_h)) return NULL;
    g_na_deadline = time(NULL) + NA_TOTAL_BUDGET_SEC;

    Expr** acc = NULL; size_t nacc = 0;
    Expr** rel = NULL; size_t nrel = 0;     /* first-integral LHSs, for the implicit pass */
    bool covered = true;                    /* every root's sub-solve gave relations */
    bool solved_roots = false;              /* the polynomial-in-y' root split ran */

    if (n < 1) {
        /* degenerate: purely algebraic in y[x] — Solve[R==0, y[x]] */
        const char* Yn = intern_symbol("DSolve`naY");
        Expr* Rsub = ds_subst(expr_copy(R), ds_make_funcapp(yname, 0, xvar),
                              expr_new_symbol(Yn));
        Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                       (Expr*[]){ Rsub, expr_new_integer(0) }, 2);
        Expr* sol = ds_solve(eq, expr_new_symbol(Yn));
        size_t ng = 0;
        Expr** gs = dsolve_extract_solutions(sol, Yn, &ng);   /* roots are the bodies */
        if (sol) expr_free(sol);
        if (gs) { acc = gs; nacc = ng; }
    } else {
        /* substitute the top derivative y^(n)[x] -> plain symbol, measure degree */
        const char* Dn = intern_symbol("DSolve`naD");
        Expr* topLit = ds_make_funcapp(yname, n, xvar);
        Expr* Rsub = ds_subst(expr_copy(R), expr_copy(topLit), expr_new_symbol(Dn));
        Expr* degE = eval_and_free(ds_call2("Exponent", expr_copy(Rsub),
                                            expr_new_symbol(Dn)));
        long deg = (degE->type == EXPR_INTEGER) ? degE->data.integer : -1;
        expr_free(degE);

        /* only solve when genuinely polynomial in the top derivative (guards
         * Solve against a radical/transcendental form the specialists own) */
        Expr* pq = eval_and_free(ds_call2("PolynomialQ", expr_copy(Rsub),
                                          expr_new_symbol(Dn)));
        bool ispoly = (pq->type == EXPR_SYMBOL && pq->data.symbol.name == SYM_True);
        expr_free(pq);

        Expr** brel = NULL; size_t nbrel = 0;
        if (deg >= 2 && ispoly && n == 1
            && na_binomial_separable(Rsub, Dn, deg, yname, xvar, &brel, &nbrel)) {
            /* binomial separable with a non-elementary y-integral: the n relations
             * are final -- stash them for the binomial implicit slot that follows
             * this one in the cascade, and skip the root recursion */
            dsolve_stash_put(NA_BINOM_STASH_TAG, memo_h, brel, nbrel);
            expr_free(Rsub);
        } else if (deg >= 2 && ispoly) {
            Expr* poly = (deg == 2 && !getenv("MATHILDA_NO_NA_SQRT")) ? expr_copy(Rsub) : NULL;
            Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                           (Expr*[]){ Rsub, expr_new_integer(0) }, 2);  /* consumes Rsub */
            Expr* sol = ds_solve(eq, expr_new_symbol(Dn));
            size_t ng = 0;
            Expr** gs = dsolve_extract_solutions(sol, Dn, &ng);
            if (sol) expr_free(sol);
            if (poly) {
                if (gs && ng == 2) na_normalise_quadratic_roots(gs, ng, poly, Dn);
                expr_free(poly);
            }
            /* A branch given as an implicit Root[...] object means Solve found no
             * explicit closed form for the top derivative (e.g. a general cubic in
             * y' with symbolic coefficients).  Recursing DSolve on y^(n) == Root[...]
             * is futile and can HANG (the implicit-derivative ODE has no method that
             * closes it, and some rewrite loops), so decline the whole method and let
             * a later specialist (Lagrange for a Lagrange/d'Alembert form, ...) take
             * the original equation. */
            bool has_root = false;
            const char* rootn = intern_symbol("Root");
            if (gs) for (size_t i = 0; i < ng; i++)
                if (ds_has_head(gs[i], rootn)) { has_root = true; break; }
            if (gs && !has_root) {
                solved_roots = true;
                /* first order only: the stashed relations are first integrals */
                bool want_rel = (n == 1);
                for (size_t i = 0; i < ng; i++) {
                    /* y^(n)[x] == g_i */
                    Expr* subeqn = expr_new_function(expr_new_symbol(SYM_Equal),
                                       (Expr*[]){ expr_copy(topLit), gs[i] }, 2);
                    recurse_collect(subeqn, yname, xvar, &acc, &nacc,
                                    want_rel ? &rel : NULL, &nrel, &covered);
                }
                if (!want_rel) covered = false;
            } else if (gs) {
                for (size_t i = 0; i < ng; i++) expr_free(gs[i]);   /* declined */
            }
            if (gs) free(gs);
        } else if (!ispoly) {
            /* The top derivative appears RATIONALLY (in a denominator) rather than
             * polynomially, e.g. A/y' == B.  Clear it: Numerator[Together[A/Dn - B]]
             * = A - B Dn is polynomial in Dn.  When that cleared numerator genuinely
             * carries the top derivative (degree >= 1 in Dn), recurse DSolve on the
             * cleared ODE (top derivative restored) and let the full cascade own it --
             * the linear normalizer handles a y-dependent leading coefficient (y^3 x'
             * = 1 - 4 y^2 x) and Lagrange/NthAlgebraic-poly a higher-degree cleared
             * form.  This closes the 12000.org "x = x(y)" spelling (A/x'[y] == B).
             * The cleared form has the derivative polynomially, so re-entry can never
             * clear again (no infinite recursion), and a normal polynomial-in-y' ODE
             * (ispoly true) never reaches here, so nothing is hijacked. */
            Expr* numDn = eval_and_free(ds_call1("Numerator",
                              eval_and_free(ds_call1("Together", expr_copy(Rsub)))));
            Expr* npq = eval_and_free(ds_call2("PolynomialQ", expr_copy(numDn),
                                               expr_new_symbol(Dn)));
            bool npoly = (npq->type == EXPR_SYMBOL && npq->data.symbol.name == SYM_True);
            expr_free(npq);
            /* The cleared numerator must still carry the top derivative (else clearing
             * produced no ODE).  Tested with FreeQ, NOT Exponent: Exponent[expr, Dn]
             * mis-returns 0 whenever expr contains any function application (Sin[y],
             * x[y], ...), which every "x = x(y)" residual does. */
            Expr* fq = eval_and_free(ds_call2("FreeQ", expr_copy(numDn),
                                              expr_new_symbol(Dn)));
            bool has_deriv = (fq->type == EXPR_SYMBOL && fq->data.symbol.name == SYM_False);
            expr_free(fq);
            if (npoly && has_deriv) {
                /* restore the top derivative (Dn -> y^(n)[x]) and recurse on cleared==0 */
                Expr* clearedR = ds_subst(numDn, expr_new_symbol(Dn), expr_copy(topLit));
                Expr* subeqn = expr_new_function(expr_new_symbol(SYM_Equal),
                                   (Expr*[]){ clearedR, expr_new_integer(0) }, 2);
                recurse_collect(subeqn, yname, xvar, &acc, &nacc, NULL, NULL, NULL);
            } else {
                expr_free(numDn);
            }
            expr_free(Rsub);
        } else {
            expr_free(Rsub);
        }
        expr_free(topLit);
    }

    /* Stash the relations only when the explicit pass finds nothing: an explicit
     * answer (even a partial union) is returned exactly as before. */
    if (nacc == 0 && solved_roots && covered && nrel > 0) {
        dsolve_stash_put(NA_STASH_TAG, memo_h, rel, nrel);
        rel = NULL; nrel = 0;
    }
    for (size_t i = 0; i < nrel; i++) expr_free(rel[i]);
    free(rel);

    if (nacc == 0) { if (acc) free(acc); na_memo_add(memo_h); return NULL; }
    *nbranch = nacc;
    return acc;
}

/* DSolve`NthAlgebraicImplicit: the implicit second pass (see the stash above).  In
 * the cascade it consumes what the explicit pass stashed for this residual; pinned
 * on its own it runs the explicit pass first to populate the stash. */
static Expr** na_implicit_impl(DSolveProblem* P, size_t* nbranch, bool may_solve) {
    if (getenv("MATHILDA_NO_NA_IMPLICIT")) return NULL;   /* A/B killswitch */
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 1) return NULL;
    uint64_t h = expr_hash(P->eq_residuals[0]);
    size_t n = 0;
    Expr** G = dsolve_stash_take(NA_STASH_TAG, h, &n);
    if (!G) G = dsolve_stash_take(NA_BINOM_STASH_TAG, h, &n);
    if (!G && may_solve) {
        size_t ne = 0;
        Expr** ex = dsolve_nth_algebraic_try(P, &ne);   /* populates the stash on decline */
        if (ex) { for (size_t i = 0; i < ne; i++) expr_free(ex[i]); free(ex); return NULL; }
        G = dsolve_stash_take(NA_STASH_TAG, h, &n);
        if (!G) G = dsolve_stash_take(NA_BINOM_STASH_TAG, h, &n);
    }
    if (!G) return NULL;
    *nbranch = n;
    return G;
}

/* The binomial-separable relations (na_binomial_separable), consumed by the slot
 * right after the explicit pass: they are the whole answer for that form, and every
 * specialist between here and the late implicit slot would only re-attempt the same
 * non-elementary separation (Separable's sampled split: ~3 s on 5612). */
Expr** dsolve_nth_algebraic_binomial_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 1) return NULL;
    size_t n = 0;
    Expr** G = dsolve_stash_take(NA_BINOM_STASH_TAG, expr_hash(P->eq_residuals[0]), &n);
    if (!G) return NULL;
    *nbranch = n;
    return G;
}

/* True when the explicit pass stashed a first integral for every y'-root of this
 * equation -- the late implicit slot will answer it, so a parametric method placed
 * before that slot defers. */
bool dsolve_nth_algebraic_has_relations(const DSolveProblem* P) {
    return P->neq == 1 && dsolve_stash_has(NA_STASH_TAG, expr_hash(P->eq_residuals[0]));
}

Expr** dsolve_nth_algebraic_implicit_try(DSolveProblem* P, size_t* nbranch) {
    return na_implicit_impl(P, nbranch, false);
}
static Expr** na_implicit_pinned_try(DSolveProblem* P, size_t* nbranch) {
    return na_implicit_impl(P, nbranch, true);
}

static Expr* builtin_dsolve_nth_algebraic(Expr* res) {
    return dsolve_method_builtin(res, dsolve_nth_algebraic_try);
}
static Expr* builtin_dsolve_nth_algebraic_implicit(Expr* res) {
    return dsolve_method_builtin_implicit(res, na_implicit_pinned_try);
}

void dsolve_nth_algebraic_init(void) {
    symtab_add_builtin("DSolve`NthAlgebraic", builtin_dsolve_nth_algebraic);
    symtab_get_def("DSolve`NthAlgebraic")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`NthAlgebraic",
        "DSolve`NthAlgebraic[eqn, y, x] solves an ODE algebraic (degree >= 2) in its "
        "highest derivative by solving for that derivative and recursing on each "
        "branch; also the degenerate purely-algebraic case.");
    symtab_add_builtin("DSolve`NthAlgebraicImplicit", builtin_dsolve_nth_algebraic_implicit);
    symtab_get_def("DSolve`NthAlgebraicImplicit")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`NthAlgebraicImplicit",
        "DSolve`NthAlgebraicImplicit[eqn, y, x] solves a first-order ODE algebraic "
        "(degree >= 2) in y' whose y'-branches solve only implicitly, returning one "
        "verified first integral G[x, y[x]] == C[1] per branch.");
}
