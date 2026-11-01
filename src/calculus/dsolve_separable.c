/*
 * dsolve_separable.c — DSolve`Separable (explicit + implicit twin).
 *
 * Both paths are reached through the single registered builtin `DSolve`Separable`
 * (there is no `DSolve`SeparableImplicit` symbol — an earlier version of this
 * comment claimed one, which is a trap when bisecting a case with pinned methods).
 *
 * Solves the first-order separable ODE  y'[x] == g(x) h(y).  The RHS F(x, Y)
 * (Y a plain symbol standing for y[x]) is separated by evaluating it at sample
 * points: for X0, Y0 with F(X0,Y0) != 0,
 *
 *     g(x) = F(x, Y0),   h(Y) = F(X0, Y) / F(X0, Y0),
 *
 * which reconstructs F exactly when F is genuinely separable (verified by
 * PossibleZeroQ[F - g h]).  Then Integrate[1/h, Y] == Integrate[g, x] + C[1].
 *
 * Two dispatch entries share the separation search (`sep_find_split`):
 *   - dsolve_separable_try           solves the relation EXPLICITLY for y.
 *   - dsolve_separable_implicit_try  returns the first integral
 *         G = Integrate[1/h, y] - Integrate[g, x]   (relation G == C[1])
 *     when the relation is elementary but not invertible for y (Log, Root, or a
 *     transcendental spelling), so it is emitted as an implicit solution rather
 *     than declined.  This is the first-order analogue of the Exact /
 *     ExactImplicit twin, and also carries the IVP fit (C = G(x0, y0), no
 *     inversion needed) for the Root-form cubic separables.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include "integrate.h"          /* g_integrate_quiet: silence speculative nonelem */
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Fast numeric pre-filter for the separability check  chk = F - g*h  (in x and
 * Y): returns false only when chk evaluates to a clearly NON-zero finite number
 * at a generic real sample -- a cheap reject of a non-separable F that avoids the
 * expensive symbolic zero-test on a transcendental residual (the split search's
 * whole cost on e.g. the E^(x y) / E^(2y) exact families).  Returns true (defer
 * to the symbolic ds_is_zero) when chk is ~0 or does not numericize (a residual
 * carrying an extra free parameter), so it NEVER rejects a genuine split. */
/* Distinct free PARAMETER symbols of e (argument positions only, never heads;
 * excluding x, Y and protected system symbols such as Pi, E, I) into names[0..*n),
 * capped at `cap`. */
static void sep_collect_params(const Expr* e, const char* xvar, const char* Yname,
                               const char** names, size_t* n, size_t cap) {
    if (!e || *n >= cap) return;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        if (nm == xvar || nm == Yname) return;
        SymbolDef* d = symtab_lookup(nm);
        if (d && (d->attributes & ATTR_PROTECTED)) return;
        for (size_t i = 0; i < *n; i++) if (names[i] == nm) return;
        names[(*n)++] = nm;
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        sep_collect_params(e->data.function.args[i], xvar, Yname, names, n, cap);
}

/* Cheap numeric pre-screen for the split check F - g h.  Evaluated at a few
 * (x, Y) points AND with every free parameter at a moderate distinct value: the
 * exact zero test alone is not enough, because its sampler can certify an
 * expression that merely agrees asymptotically -- with a symbolic exponent n,
 * F - g h for x (a + x y^n) y' + b y == 0 (§2.2.54-5331) has both terms tending to
 * -b/(x^2 Y^(n-1)) at large n and was judged identically zero.  That accepted a
 * wrong split whose integrals then spent the whole shared separable budget,
 * starving every nested Separable call of the same top-level solve. */
static bool sep_chk_maybe_zero(const Expr* chk, const char* xvar, const char* Yname) {
    static const double xs[] = { 1.3, 0.7 };
    static const double ys[] = { 0.9, 1.7 };
    const char* params[12]; size_t np = 0;
    if (!getenv("MATHILDA_NO_SEP_PARAMCHK"))                     /* A/B killswitch */
        sep_collect_params(chk, xvar, Yname, params, &np, 12);
    for (int i = 0; i < 2; i++) {
        Expr* e = expr_copy((Expr*)chk);
        e = ds_subst(e, expr_new_symbol(xvar),  expr_new_real(xs[i]));
        e = ds_subst(e, expr_new_symbol(Yname), expr_new_real(ys[i]));
        for (int with_params = 0; with_params < 2; with_params++) {
            Expr* ev = expr_copy(e);
            if (with_params) {
                if (np == 0) { expr_free(ev); break; }
                for (size_t k = 0; k < np; k++)
                    ev = ds_subst(ev, expr_new_symbol(params[k]),
                                  expr_new_real(1.13 + 0.37 * (double)k + 0.11 * i));
            }
            ev = eval_and_free(ds_call1("Abs", ev));
            double m = (ev && ev->type == EXPR_REAL)    ? ev->data.real
                     : (ev && ev->type == EXPR_INTEGER) ? (double)ev->data.integer : NAN;
            expr_free(ev);
            if (!isnan(m) && isfinite(m) && m > 1e-6) { expr_free(e); return false; }  /* clearly non-zero */
        }
        expr_free(e);
    }
    return true;                                                  /* small / non-numeric */
}

/* Find a separation F(x,Y) = g(x) h(Y).  On success returns true and fills *g_out
 * / *h_out (owned by the caller) and *Yname_out (the interned dummy for y[x]).
 * The RHS F is solved for from P and copied out via *FY_out (also owned) so the
 * implicit path can reuse it.  Returns false (nothing allocated) when P is not a
 * first-order single ODE or no separable split is found. */
static bool sep_find_split(DSolveProblem* P, Expr** g_out, Expr** h_out,
                           const char** Yname_out) {
    if (P->nfun != 1 || P->neq != 1) return false;
    if (P->max_order[0] != 1) return false;
    const char* xvar = P->ind_names[0];
    const char* yname = P->fun_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return false;

    const char* Yname = intern_symbol("DSolve`Y");
    Expr* yx = ds_make_funcapp(yname, 0, xvar);
    Expr* FY = ds_subst(F, yx, expr_new_symbol(Yname));   /* consumes F, yx */

    Expr* g = NULL; Expr* h = NULL;

    /* STRUCTURAL split first: when F is literally a product whose every factor is
     * free of x or free of Y, g = the Y-free factors and h = the x-free ones, with
     * no sampling.  The sampled split below reconstructs g = F(x, Y0) and
     * h = F(X0, Y)/F(X0, Y0), which leaves the sample's trig values as unreduced
     * constants -- y' + Sin[2x] Csc[2y] == 0 (§2.2.48-4725) splits into
     * -Sin[2x] Csc[4] and Sin[4] Csc[2Y], the Csc[4] Sin[4] never cancels, and the
     * explicit Solve inversion of the resulting relation fails.  The structural
     * split gives the clean -Sin[2x] and Csc[2Y].  (A factor depending on both
     * variables leaves it to the sampling search, unchanged.) */
    {
        size_t nf = head_is(FY, SYM_Times) ? FY->data.function.arg_count : 1;
        Expr** gx = malloc(nf * sizeof(Expr*)); size_t ngx = 0;
        Expr** hy = malloc(nf * sizeof(Expr*)); size_t nhy = 0;
        bool ok = true;
        for (size_t i = 0; i < nf && ok; i++) {
            Expr* fi = head_is(FY, SYM_Times) ? FY->data.function.args[i] : FY;
            if (!ds_contains(fi, Yname)) gx[ngx++] = expr_copy(fi);
            else if (!ds_contains(fi, xvar)) hy[nhy++] = expr_copy(fi);
            else ok = false;
        }
        if (ok && ngx > 0 && nhy > 0) {
            g = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), gx, ngx));
            h = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), hy, nhy));
        } else {
            for (size_t i = 0; i < ngx; i++) expr_free(gx[i]);
            for (size_t i = 0; i < nhy; i++) expr_free(hy[i]);
        }
        free(gx); free(hy);
    }

    static const int samples[] = { 2, 3, 5, 7, -2, -3 };
    size_t nsamp = sizeof(samples) / sizeof(samples[0]);
    for (size_t ai = 0; ai < nsamp && !g; ai++) {
        for (size_t bi = 0; bi < nsamp && !g; bi++) {
            int X0 = samples[ai], Y0 = samples[bi];
            Expr* denom = ds_subst(
                ds_subst(expr_copy(FY), expr_new_symbol(xvar), expr_new_integer(X0)),
                expr_new_symbol(Yname), expr_new_integer(Y0));
            /* Skip a sample only when the denominator F(X0,Y0) is PROVABLY zero.
             * The old gate demanded provably NON-zero, which rejected every sample
             * of a symbolic-parameter RHS such as (a y+b)/(c y+d) (F(X0,Y0) is not
             * a decidable nonzero) and declined a genuinely separable equation.
             * A generically-nonzero symbolic denominator reconstructs the split,
             * and the F - g h == 0 check below rejects a wrong one regardless. */
            if (ds_is_zero(denom)) { expr_free(denom); continue; }
            Expr* gcand = ds_subst(expr_copy(FY), expr_new_symbol(Yname), expr_new_integer(Y0));
            Expr* hx0 = ds_subst(expr_copy(FY), expr_new_symbol(xvar), expr_new_integer(X0));
            Expr* hcand = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                hx0,
                expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy(denom), expr_new_integer(-1) }, 2)
            }, 2));
            expr_free(denom);
            Expr* prod = ds_call2(SYM_Times, expr_copy(gcand), expr_copy(hcand));
            Expr* check = eval_and_free(ds_call2(SYM_Subtract, expr_copy(FY), prod));
            if (sep_chk_maybe_zero(check, xvar, Yname) && ds_is_zero(check)) {
                g = gcand; h = hcand;
            } else { expr_free(gcand); expr_free(hcand); }
            expr_free(check);
        }
    }
    expr_free(FY);
    if (!g) return false;
    *g_out = g; *h_out = h; *Yname_out = Yname;
    return true;
}

/* The two integrals ∫1/h dY and ∫g dx; consumes g and h.  Returns true and fills
 * *lhs / *rhs (owned).  Muted so a speculative non-elementary integrand does not
 * leak Integrate::nonelem on the DSolve`Y dummy.  When `require_elem` is true
 * (the EXPLICIT path, which must invert the relation for y) a non-elementary
 * integral fails; when false (the IMPLICIT path) a non-elementary integral is
 * KEPT unevaluated, so the first integral G == C is still emitted (and verified
 * by the implicit-function rule, since D[Integrate[f,y],y] == f) -- e.g. the
 * autonomous y' == -2 ArcTan[y]/(1+y^2), whose y-side integral is non-elementary. */
/* Per-antiderivative budget, in seconds.
 *
 * WHY THIS EXISTS.  Separable sits at cascade slot 7 of 52 and had NO bound on its
 * own antiderivatives, so one integrand that does not close spent the entire DSolve
 * budget on behalf of every method behind it.  And the integrands here are not the
 * ones the user wrote: the split search samples F at integer points, which turns a
 * trig kernel with a compound argument into a SHIFTED kernel carrying a constant trig
 * coefficient — `Csc[2] Csc[x] Cos[x-2]` for g, and the quotient
 * `Cot[2]^2/(Csc[2] Csc[Y] Cos[2-Y] - 1)` for 1/h — each of which Integrate searches
 * for ~20 s before coming back unevaluated.
 *
 * Measured on §2.2.36-3521 `y' == Cos[x-y]/(Sin[x] Sin[y]) - 1` and its two siblings
 * 3527/3599, which this repairs: before, `DSolve`Separable` ran 30 s without deciding
 * and the solve died at the 8 s corpus bound; after, it declines in ~4 s and
 * `DSolve`Linearizable` — three slots later, and able to solve the ORIGINAL equation
 * all along — answers `ArcCos[C[1] Csc[x]]` in 0.1 s.  So these records were never
 * missing a capability.  One method was spending the budget of the method that had
 * the answer, and the fix is to bound the attempt rather than to teach Separable a
 * new trick: a bounded attempt helps every integrand nobody has thought of, and a
 * rewrite (TrigExpand before sampling) was built, measured against this, and
 * DISCARDED — it was faster on these three (0.08 s) but claimed them for Separable's
 * implicit twin, whose answer carries the sampling artefact `Cot[2]` and is an
 * implicit relation where Mathematica and Linearizable give the explicit ArcCos.
 *
 * It is ONE DEADLINE PER TOP-LEVEL SOLVE, shared by the explicit entry and the
 * implicit twin, and not a budget per integral — which is the difference between a
 * fix and a flaky one.  A per-integral bound has to serve two measured requirements
 * at once and cannot: §2.1.2-1134's answer legitimately carries an unevaluated
 * `Integrate` whose sampled integrand needs ~1-2 s to come back undecided, so a 1 s
 * bound LOSES that record (measured: an answer at 5.1 s became a 16 s abort) and a
 * 2 s bound leaves it sitting on the boundary, passing or failing with the load.
 * Raising the bound instead pushes the records this repairs toward the wall, because
 * each of the two paths pays its own: at 2 s §2.2.36-3521 closes at 4.1 s, at 3 s at
 * 6.2 s, against an 8 s budget.  A shared deadline separates the two cases on the
 * quantity that actually distinguishes them — 1134's integrals cost seconds, 3521's
 * cost twenty — and charges a hopeless integrand ONCE: at six seconds 1134 is
 * unchanged at 5.3 s (measured twice, where the deadline never binds it) and 3521
 * declines at the deadline and closes at 6.15 s, after which Linearizable answers.
 * Six is therefore the smallest value that clears 1134's genuine ~5 s need, and the
 * sharing is what keeps 3521 inside the 8 s bound at that value rather than paying
 * 12 s for two paths.
 *
 * A TIMEOUT IS TREATED AS A DECLINE ON BOTH PATHS, including the implicit twin that
 * normally keeps a non-elementary integral unevaluated.  The distinction matters: a
 * non-elementary integrand has been *decided*, and an inert `Integrate` for it is an
 * answer; a timed-out one has been decided of nothing, and handing back a raw
 * unevaluated `Integrate` re-enters the integration cascade on every later
 * re-evaluation of the body — the trap that turned a 3 s budget into 375 s in M62. */
#define SEP_BUDGET_SEC 6

static uint64_t g_sep_epoch = 0;
static bool     g_sep_epoch_set = false;
static time_t   g_sep_deadline = 0;

/* Seconds left in this top-level solve's separable budget; 0 when it is spent.  The
 * deadline is (re)armed on the first call of each top-level evaluation, so the
 * explicit and implicit entries of one DSolve share it. */
static int sep_remaining(void) {
    uint64_t ep = eval_toplevel_id();
    if (!g_sep_epoch_set || ep != g_sep_epoch) {
        g_sep_epoch = ep; g_sep_epoch_set = true;
        g_sep_deadline = time(NULL) + SEP_BUDGET_SEC;
    }
    time_t now = time(NULL);
    return (g_sep_deadline > now) ? (int)(g_sep_deadline - now) : 0;
}

/* Integrate[integrand, var] inside what is left of the budget.  Consumes
 * `integrand`; returns NULL and sets *spent when the deadline is reached. */
static Expr* sep_integrate_bounded(Expr* integrand, const char* var, bool* spent) {
    int left = sep_remaining();
    if (left <= 0) { expr_free(integrand); *spent = true; return NULL; }
    Expr* call = expr_new_function(expr_new_symbol(SYM_Integrate),
                     (Expr*[]){ integrand, expr_new_symbol(var) }, 2);
    Expr* guarded = expr_new_function(expr_new_symbol(SYM_TimeConstrained),
                     (Expr*[]){ call, expr_new_integer(left),
                                expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    Expr* out = eval_and_free(guarded);
    if (out && out->type == EXPR_SYMBOL
        && out->data.symbol.name == intern_symbol("$Aborted")) {
        expr_free(out);
        *spent = true;
        return NULL;
    }
    return out;
}

static bool sep_integrals(Expr* g, Expr* h, const char* Yname, const char* xvar,
                          Expr** lhs, Expr** rhs, bool require_elem) {
    Expr* invh = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ h, expr_new_integer(-1) }, 2));
    bool spent = false;
    g_integrate_quiet++;
    Expr* lhsInt = sep_integrate_bounded(invh, Yname, &spent);
    Expr* rhsInt = NULL;
    if (spent) expr_free(g);            /* the second attempt never ran: g is still ours */
    else       rhsInt = sep_integrate_bounded(g, xvar, &spent);   /* consumes g regardless */
    g_integrate_quiet--;
    if (spent) {
        expr_free(lhsInt); expr_free(rhsInt);
        return false;
    }
    if (require_elem && (ds_has_head(lhsInt, SYM_Integrate) || ds_has_head(rhsInt, SYM_Integrate))) {
        expr_free(lhsInt); expr_free(rhsInt); return false;
    }
    *lhs = lhsInt; *rhs = rhsInt;
    return true;
}

/* Collect the DISTINCT arguments of Log[..] subexpressions of e into args[0..*n)
 * (borrowed pointers, de-duplicated by expr_eq), capped at `cap`.  Does not recurse
 * into a Log's own argument (a nested Log is not expected in a separable integral). */
static void sep_collect_log_args(const Expr* e, Expr** args, size_t* n, size_t cap) {
    if (!e || e->type != EXPR_FUNCTION || *n >= cap) return;
    Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL
        && head->data.symbol.name == SYM_Log
        && e->data.function.arg_count == 1) {
        Expr* a = e->data.function.args[0];
        for (size_t i = 0; i < *n; i++) if (expr_eq(args[i], a)) return;
        args[(*n)++] = a;
        return;
    }
    sep_collect_log_args(head, args, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        sep_collect_log_args(e->data.function.args[i], args, n, cap);
}

/* True when the LHS antiderivative Integrate[1/h, Y] carries THREE OR MORE distinct
 * Log[..] arguments -- the partial-fraction integral of a cubic-or-higher h(Y),
 * e.g. -1/2 Log[Y-1] + 1/6 Log[Y+1] + 1/3 Log[Y-2] (however the coefficients are
 * spelled: outside the Log, or absorbed as Log[(Y-1)^(-1/2)] powers).  The relation
 * Sum c_i Log[Y-r_i] == RHS has no elementary explicit inverse for Y that ds_solve
 * closes -- it churns -- so the caller declines and the implicit first-integral twin
 * emits G == C[1] (2.2.16-1590).  A ONE- or TWO-log relation is left on the explicit
 * path: a single Log inverts by Exp, and two logs merge to a Mobius Log[(Y-a)/(Y-b)]
 * that inverts (Tanh / logistic). */
static bool sep_noninvertible_logsum(const Expr* lhsInt) {
    Expr* args[8]; size_t n = 0;
    sep_collect_log_args(lhsInt, args, &n, 8);
    return n >= 3;
}

/* True iff e carries a genuine ROOT -- Sqrt[u] or u^(p/q) with q > 1 and u not a
 * number.  (ds_has_radical_power also counts E^Y and symbolic exponents, which are
 * single-valued: inverting E^Y == Log[t] + C to Y == Log[Log[t] + C] is exact.) */
static bool sep_has_root(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && e->data.function.arg_count == 2
        && h->data.symbol.name == SYM_Power) {
        const Expr* b  = e->data.function.args[0];
        const Expr* ex = e->data.function.args[1];
        bool numbase = b->type == EXPR_INTEGER || b->type == EXPR_BIGINT || b->type == EXPR_REAL;
        if (!numbase && head_is((Expr*)ex, SYM_Rational)) return true;
    }
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Sqrt) return true;
    if (sep_has_root(h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (sep_has_root(e->data.function.args[i])) return true;
    return false;
}

/* True when the antiderivative carries TWO OR MORE distinct Log[..] arguments in Y,
 * at least one of them wrapping a root of Y -- e.g. Log[u] - Log[4 + 2 u + 2 Sqrt[
 * u (4 + u)]] (the separable branch of an isobaric quadratic, §2.2.55-5429).  Solve
 * can only invert it by exponentiating and rationalising the radical away, which
 * took 5 s to reach a degree-16 Root in E^(t + C) -- a worse answer than the first
 * integral, and long enough to time out every caller.  A single Log wrapping a root
 * (Log[Y + Sqrt[1 + Y^2]], the ArcSinh shape) still inverts cleanly and is left
 * alone.  The caller declines; the implicit first-integral twin answers. */
static bool sep_log_radical_mix(const Expr* lhsInt, const char* Yname) {
    Expr* args[8]; size_t n = 0;
    sep_collect_log_args(lhsInt, args, &n, 8);
    size_t ny = 0; bool rooted = false;
    for (size_t i = 0; i < n; i++) {
        if (!ds_contains(args[i], Yname)) continue;
        ny++;
        if (sep_has_root(args[i])) rooted = true;
    }
    return ny >= 2 && rooted;
}

Expr** dsolve_separable_try(DSolveProblem* P, size_t* nbranch) {
    const char* xvar = P->ind_names[0];
    Expr* g = NULL; Expr* h = NULL; const char* Yname = NULL;
    if (!sep_find_split(P, &g, &h, &Yname)) return NULL;

    /* A RADICAL in h(Y) whose explicit inversion removes it (y' x == Sqrt[1 - y^2]
     * -> ArcSin[y] == Log[x] + C -> y == Sin[Log[x] + C]) can yield a form valid on
     * ONE branch of the root only -- there y' == Cos[..]/x while the equation demands
     * |Cos[..]|/x.  Such an explicit form is dropped (test below) so the implicit
     * first integral (ArcSin[y] - Log[x] == C), valid on every branch, answers.  (A
     * body that still carries the radical keeps the existing per-branch numeric
     * check.)  Measured: §2.2.12-1134, §2.2.35-3412, §2.2.24-2360/§2.2.26-2536,
     * where the explicit IVP branch was a wrong answer. */
    bool h_radical = sep_has_root(h);

    Expr* lhsInt = NULL; Expr* rhsInt = NULL;
    if (!sep_integrals(g, h, Yname, xvar, &lhsInt, &rhsInt, true)) return NULL;

    /* Skip a non-invertible >=3-distinct-log relation (cubic+ h(Y)): ds_solve
     * churns on it, so decline here and let the implicit first-integral twin
     * (dsolve_separable_implicit_try, next in the cascade) return G == C[1]. */
    if (sep_noninvertible_logsum(lhsInt)
        || (!getenv("MATHILDA_NO_SEP_LOGRAD") && sep_log_radical_mix(lhsInt, Yname))) {
        expr_free(lhsInt); expr_free(rhsInt);
        return NULL;
    }

    /* Integrate[1/h, Y] == Integrate[g, x] + C[1], solved explicitly for Y. */
    Expr* rhs = eval_and_free(ds_call2(SYM_Plus, rhsInt, ds_const(1)));
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ lhsInt, rhs }, 2);

    Expr* solres = ds_solve(eq, expr_new_symbol(Yname));
    size_t nb = 0;
    Expr** bodies = dsolve_extract_solutions(solres, Yname, &nb);
    if (solres) expr_free(solres);
    if (!bodies) return NULL;

    /* Drop a spurious inversion body that is free of the independent variable: a
     * genuine general solution of a non-autonomous separable ODE always depends on
     * x, so an x-free body (e.g. the bare integration constant C[1] a collapsed
     * multivalued inverse can produce) is not a solution.  Falling through to the
     * implicit fallback yields the correct first integral instead. */
    size_t keep = 0;
    bool ivp = (P->ncond > 0);
    for (size_t i = 0; i < nb; i++) {
        if (!bodies[i]) continue;
        /* Drop a spurious x-free inversion body (a collapsed constant is not a
         * general solution).  On an IVP, also drop a Root-form body: the explicit
         * fitter cannot Solve[Root[..C..]==y0, C], so leaving them all out (keep==0
         * -> decline) lets the implicit twin fit the constant on the first integral
         * G(x0,y0) with no inversion (the Root-form cubic separables 1149/1150). */
        /* Also drop a branch whose residual is CONFIDENTLY nonzero against the original
         * ODE: the inversion of a fractional-power relation (Integrate[1/h,Y]==... with
         * h carrying a cube/square root) emits both +/- roots, and the branch invalid
         * under the principal root is nonzero (often complex) even though its symbolic
         * residual is a branch-cut expression zero_test keeps -- y'==3x(y-1)^(1/3) emits
         * a spurious complex twin beside the correct 1+x^3 (§2.2.17-1622/1624).  Kept
         * conservative in ds_branch_num_ok so a partial-domain-valid branch is never lost. */
        bool drop = ds_free_of(bodies[i], xvar) ||
                    (ivp && ds_contains(bodies[i], intern_symbol("Root"))) ||
                    (ds_has_radical_power(bodies[i]) && !ds_branch_num_ok(P, bodies[i]));
        if (drop) { expr_free(bodies[i]); bodies[i] = NULL; }
        else bodies[keep++] = bodies[i];
    }
    if (keep > 0 && h_radical) {
        bool any_radical_body = false;
        for (size_t i = 0; i < keep; i++)
            if (sep_has_root(bodies[i])) any_radical_body = true;
        /* The inversion removed the radical: keep the explicit form only if it is
         * valid on EVERY branch, i.e. its residual vanishes for every value of the
         * constant -- sampled at a spread of C[1] values through the numeric reject
         * filter.  Sinh[x + C] for p' == Sqrt[1 + p^2] passes (Cosh is never
         * negative; §2.2.33-3257 needs the explicit p); Sin[Log[x] + C] fails at
         * any C putting Cos negative on the grid, and defers to the implicit form. */
        if (!any_radical_body) {
            static const double cs[] = { 0.37, 1.9, 2.8, -1.4 };
            bool global = true;
            for (size_t i = 0; i < keep && global; i++)
                for (size_t k = 0; k < sizeof(cs) / sizeof(cs[0]) && global; k++) {
                    Expr* bc = ds_subst(expr_copy(bodies[i]), ds_const(1), expr_new_real(cs[k]));
                    if (!ds_branch_num_ok(P, bc)) global = false;
                    expr_free(bc);
                }
            if (!global) {
                for (size_t i = 0; i < keep; i++) expr_free(bodies[i]);
                keep = 0;
            }
        }
    }
    if (keep == 0) { free(bodies); return NULL; }
    *nbranch = keep;
    return bodies;
}

/* Implicit first integral G = Integrate[1/h, y] - Integrate[g, x]; relation
 * G == C[1].  Reached (after the explicit try) when the relation is elementary
 * but not invertible for y (Log / Root / transcendental).  Returned in terms of
 * y[x] so dsolve_run_implicit's implicit-function-rule verify and (x0,y0) fit
 * apply directly. */
Expr** dsolve_separable_implicit_try(DSolveProblem* P, size_t* nbranch) {
    const char* xvar = P->ind_names[0];
    const char* yname = P->fun_names[0];
    Expr* g = NULL; Expr* h = NULL; const char* Yname = NULL;
    if (!sep_find_split(P, &g, &h, &Yname)) return NULL;

    Expr* lhsInt = NULL; Expr* rhsInt = NULL;
    /* IMPLICIT path: keep a non-elementary integral unevaluated (require_elem=false)
     * so the first integral is still returned as an implicit solution. */
    if (!sep_integrals(g, h, Yname, xvar, &lhsInt, &rhsInt, false)) return NULL;

    /* G(x, y) = (∫1/h dY)|_{Y -> y[x]} - ∫g dx */
    Expr* lhsY = ds_subst(lhsInt, expr_new_symbol(Yname), ds_make_funcapp(yname, 0, xvar));
    Expr* G = eval_and_free(ds_call2(SYM_Subtract, lhsY, rhsInt));

    Expr** out = malloc(sizeof(Expr*));
    out[0] = G;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_separable(Expr* res) {
    Expr* r = dsolve_method_builtin(res, dsolve_separable_try);
    if (!r) r = dsolve_method_builtin_implicit(res, dsolve_separable_implicit_try);
    return r;
}

void dsolve_separable_init(void) {
    symtab_add_builtin("DSolve`Separable", builtin_dsolve_separable);
    symtab_get_def("DSolve`Separable")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Separable",
        "DSolve`Separable[eqn, y, x] solves y'[x] == g(x) h(y) by separating "
        "variables: Integrate[1/h, y] == Integrate[g, x] + C[1], solved for y "
        "(or returned as the implicit first integral when it does not invert).");
}
