/* integrate_derivdivides.c
 *
 * Integration by substitution -- the "derivative-divides" method.  See
 * integrate_derivdivides.h for the high-level description.
 *
 * ---------------------------------------------------------------------
 * Algorithm
 * ---------------------------------------------------------------------
 * Candidate kernels u(x) are every distinct subexpression of f that depends
 * on x, except x itself and f itself (this reproduces
 * `Cases[Union[Level[f, {0, Infinity}]], e_ /; !FreeQ[e, x]]` from the user's
 * worked transcript, but is collected directly in C so we never round-trip
 * through the parser).  For each kernel u(x):
 *
 *   Direct quotient strategy
 *   ------------------------
 *     du = D[u(x), x];   if du == 0 skip
 *     q  = Cancel[Together[f / du]]
 *     qu = q /. u(x) -> u            (u a fresh symbol)
 *     if FreeQ[qu, x]:
 *         G  = Integrate[qu, u]      (recurse into the full integrator)
 *         if G closes (no leftover Integrate):
 *             r = G /. u -> u(x)
 *             if PossibleZeroQ[D[r, x] - f]:  return r
 *
 *   Eliminate / Solve strategy (explicit method only)
 *   -------------------------------------------------
 *     rel  = Eliminate[{Dt[y] == f Dt[x], u == u(x), Dt[u == u(x)]},
 *                      {x, Dt[x]}]
 *     sol  = Solve[rel, Dt[y]]
 *     gs   = Cancel[ PowerExpand[ Factor //@ (Dt[y] /. sol) ] / Dt[u] ]
 *     for each branch g in gs (gs is a List, one entry per Solve branch):
 *         if FreeQ[g, x] && FreeQ[g, Dt[x]] && FreeQ[g, Dt[u]]:
 *             G = Integrate[g, u];  if G closes:
 *                 r = G /. u -> u(x)
 *                 if PossibleZeroQ[D[r, x] - f]:  return r
 *
 * The PossibleZeroQ[D[r, x] - f] gate is unconditional and is what selects
 * the correct inverse-function branch among Solve's +- solutions and the
 * branch assumptions PowerExpand bakes in.  PossibleZeroQ's numeric sampling
 * settles these branches decisively; we deliberately avoid a follow-up
 * Simplify confirm, which on entangled radical-trig integrands (e.g.
 * Sqrt[Tan[x]]) cost ~1.1s of trigrat normalization for no behavioural gain.
 *
 * Memory: builtins take ownership of `res`; helpers below own every Expr they
 * construct and free intermediates explicitly.  `eval_take` consumes its
 * argument.  We never expr_free(res).
 */

/* clock_gettime / CLOCK_MONOTONIC (the TimeConstraint deadline clock) are POSIX,
 * not C99.  glibc hides them under -std=c99 unless a feature-test macro is
 * defined BEFORE the first include; below it, it has no effect (SPEC.md §10). */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "integrate_derivdivides.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include "internal.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "poly/eliminate.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* ---------------------------------------------------------------------- */
/* Small builders / evaluation helpers (mirrors integrate_unknown.c)      */
/* ---------------------------------------------------------------------- */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }

static Expr* mk_fn1(const char* name, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(expr_new_symbol(name), args, 1);
}

static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(name), args, 2);
}

/* Evaluate `call` to a fixed point, freeing `call`. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

static bool is_lit_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* True if `f` contains no subexpression structurally equal to `x`
 * (the structural meaning of FreeQ[f, x]). */
static bool expr_free_of(const Expr* f, const Expr* x) {
    if (expr_eq((Expr*)f, (Expr*)x)) return false;
    if (f->type == EXPR_FUNCTION) {
        if (!expr_free_of(f->data.function.head, x)) return false;
        for (size_t i = 0; i < f->data.function.arg_count; i++) {
            if (!expr_free_of(f->data.function.args[i], x)) return false;
        }
    }
    return true;
}

/* True if `e` contains any unevaluated Integrate[...] call (the integral did
 * not close).  Also catches a leftover Integrate`DerivativeDivides[...]. */
static bool contains_unintegrated(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Integrate)) return true;
    if (contains_unintegrated(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_unintegrated(e->data.function.args[i])) return true;
    }
    return false;
}

/* True if `f` contains a Power[b, e] whose base `b` involves x and whose
 * exponent `e` is not a concrete number -- a symbolic-exponent monomial such
 * as x^k, x^(k-1) or x^(a+b).
 *
 * DerivativeDivides forms the quotient f / u'(x) and normalises it with
 * Together/Cancel.  When x carries a symbolic exponent, that normalisation
 * feeds PolynomialGCD a pseudo-remainder sequence over the symbolic powers as
 * independent generators, which blows up (unbounded time) -- and such a
 * monomial is never a productive u-substitution kernel here anyway.  Bailing
 * early keeps the cascade responsive; the integral falls through to the
 * methods that own these forms (rational / Euler-Beta reductions, tables). */
static bool has_symbolic_power_of(const Expr* f, const Expr* x) {
    if (!f || f->type != EXPR_FUNCTION) return false;
    if (head_is((Expr*)f, SYM_Power) && f->data.function.arg_count == 2) {
        const Expr* b = f->data.function.args[0];
        const Expr* e = f->data.function.args[1];
        if (!expr_free_of(b, x) && !expr_is_numeric_like(e)) return true;
    }
    if (!expr_free_of(f->data.function.head, x)) return true;
    for (size_t i = 0; i < f->data.function.arg_count; i++)
        if (has_symbolic_power_of(f->data.function.args[i], x)) return true;
    return false;
}

/* D[expr, x]; borrows, returns owned (evaluated). */
static Expr* deriv_dx(const Expr* expr, const Expr* x) {
    return eval_take(mk_fn2("D", expr_copy((Expr*)expr), expr_copy((Expr*)x)));
}

/* Cancel[Together[e]], consuming `e`.  Canonicalises a rational expression so
 * a quotient collapses (e.g. f/u' with a shared factor). */
static Expr* cancel_together(Expr* e) {
    if (!e) return NULL;
    Expr* t = eval_take(internal_together((Expr*[]){ e }, 1));
    if (!t) return NULL;
    return eval_take(internal_cancel((Expr*[]){ t }, 1));
}

/* ReplaceAll[expr, from -> to]; consumes `expr`; returns owned (evaluated). */
static Expr* replace_one(Expr* expr, const Expr* from, const Expr* to) {
    Expr* rule = mk_fn2("Rule", expr_copy((Expr*)from), expr_copy((Expr*)to));
    Expr* call = internal_replace_all((Expr*[]){ expr, rule }, 2);
    return eval_take(call);
}

/* The differentiation-verification gate: true iff D[r, x] - f is zero.
 *
 * Acceptance rests solely on PossibleZeroQ[D[r, x] - f] === True.  PossibleZeroQ
 * uses numeric sampling, which settles these candidates (including the entangled
 * radical-trig forms) decisively.  We deliberately do NOT follow up with a
 * rigorous Simplify confirm: that confirm cost ~1.1s on
 * Integrate[Sqrt[Tan[x]], x] (a heavy trigrat normalization of the
 * back-substituted difference) for no behavioural gain over the sampler.
 * Borrows r and f. */
static bool differentiates_back(const Expr* r, const Expr* f, const Expr* x) {
    Expr* dr = deriv_dx(r, x);
    if (!dr) return false;
    Expr* diff = mk_fn2("Plus", dr,
                        mk_fn2("Times", mk_int(-1), expr_copy((Expr*)f)));

    Expr* pz = eval_take(mk_fn1("PossibleZeroQ", diff));
    bool zero = pz && pz->type == EXPR_SYMBOL && pz->data.symbol.name == SYM_True;
    if (pz) expr_free(pz);
    return zero;
}

/* Recurse into the full integrator: Integrate[g, u].  Returns the closed
 * antiderivative, or NULL if it does not close.  Borrows g and u. */
static Expr* integrate_in(const Expr* g, const Expr* u) {
    Expr* r = eval_take(mk_fn2("Integrate", expr_copy((Expr*)g),
                               expr_copy((Expr*)u)));
    if (!r) return NULL;
    if (contains_unintegrated(r)) { expr_free(r); return NULL; }
    return r;
}

/* ---------------------------------------------------------------------- */
/* Candidate-kernel collection                                            */
/* ---------------------------------------------------------------------- */

typedef struct { Expr** items; size_t n, cap; } ExprVec;

static void vec_init(ExprVec* v) { v->items = NULL; v->n = 0; v->cap = 0; }
static void vec_free(ExprVec* v) {
    for (size_t i = 0; i < v->n; i++) expr_free(v->items[i]);
    free(v->items); vec_init(v);
}
static void vec_add_unique(ExprVec* v, const Expr* e) {
    for (size_t i = 0; i < v->n; i++) if (expr_eq(v->items[i], (Expr*)e)) return;
    if (v->n == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 8;
        v->items = realloc(v->items, v->cap * sizeof(Expr*));
    }
    v->items[v->n++] = expr_copy((Expr*)e);
}

/* Number of nodes in e -- used to order kernels (simpler kernels first). */
static size_t leaf_count(const Expr* e) {
    if (!e) return 0;
    size_t n = 1;
    if (e->type == EXPR_FUNCTION) {
        n += leaf_count(e->data.function.head);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            n += leaf_count(e->data.function.args[i]);
    }
    return n;
}

/* Collect every distinct subexpression of `e` that depends on x, excluding
 * x itself and the whole integrand `f`. */
static void collect_kernels(const Expr* e, const Expr* x, const Expr* f,
                            ExprVec* out) {
    if (!e || expr_free_of(e, x)) return;          /* must depend on x */
    if (!expr_eq((Expr*)e, (Expr*)x) && !expr_eq((Expr*)e, (Expr*)f)) {
        vec_add_unique(out, e);
    }
    if (e->type == EXPR_FUNCTION) {
        collect_kernels(e->data.function.head, x, f, out);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            collect_kernels(e->data.function.args[i], x, f, out);
    }
}

/* Sort kernels by ascending leaf count (stable insertion sort -- the vector
 * is small).  Simpler kernels tend to give cleaner reduced integrands; the
 * differentiation gate guarantees correctness regardless of order. */
static void sort_kernels(ExprVec* v) {
    for (size_t i = 1; i < v->n; i++) {
        Expr* key = v->items[i];
        size_t kc = leaf_count(key);
        size_t j = i;
        while (j > 0 && leaf_count(v->items[j - 1]) > kc) {
            v->items[j] = v->items[j - 1];
            j--;
        }
        v->items[j] = key;
    }
}

/* ---------------------------------------------------------------------- */
/* Recursion guard (the reduced integral re-enters Integrate)             */
/* ---------------------------------------------------------------------- */

/* Every substitution trial reduces Integrate[f, x] to Integrate[g, u] and
 * recurses through the whole integrator via integrate_in().  Derivative-divides
 * (the cheap direct-quotient strategy) stays available on those sub-integrals;
 * we do NOT switch it off when nested.  Two guards keep the recursion finite
 * and cheap:
 *
 * 1. Integrand memo (loop guard).  Overlapping substitution branches routinely
 *    regenerate integrands seen earlier in the descent: Integrate[x Sin[x^2], x]
 *    spawns Integrate[Sin[u]/2, u], whose Eliminate/Solve branches reduce to a
 *    handful of radical forms that fold back into one another.  Without
 *    memoisation that fans out into an exponential tree -- an effective hang --
 *    even though the number of *distinct* integrands is tiny.  We therefore
 *    memoise every integrand attempted in the current top-level descent.  Each
 *    level mints a fresh substitution variable, so the recurrence is invisible
 *    to plain structural equality; we canonicalise by renaming the integration
 *    variable to a fixed sentinel before comparing.  A canonical integrand seen
 *    before short-circuits to NULL (its integrability does not depend on the
 *    path that reached it, so re-attempting it can only loop or repeat work).
 *    The memo is owned by the outermost dd_core frame and freed when it
 *    returns, so independent top-level integrals never share state.
 *
 * 2. Eliminate/Solve only at the outermost call (see dd_core).  The memo makes
 *    the recursion terminate, but the Eliminate + Solve + Factor + PowerExpand
 *    branch-search costs ~0.1-1s per node, so even a memo-bounded chain of a
 *    dozen radical nodes is a multi-second near-hang.  That heavyweight search
 *    only earns its keep on the *original* integrand (e.g. u = Sqrt[Tan[x]],
 *    whose reduction is rational and is finished by the rational stage without
 *    re-entering dd); reduced sub-integrals are better finished by the cheaper
 *    cascade stages.  Restricting it to the outermost frame takes the
 *    formerly-hanging cases to a few milliseconds while leaving every closed
 *    radical substitution intact.
 *
 * dd_depth remains as a hard backstop ceiling. */
#define DD_MAX_DEPTH 8
static int dd_depth = 0;

/* ------------------------------------------------------------------------
 * TimeConstraint deadline (seconds), set per top-level Integrate by the
 * cascade; HUGE_VAL means unbounded, which is what the explicit
 * `Method -> "DerivativeDivides"` surface gets.
 *
 * WHY THIS STAGE AND NOT ANOTHER. Measured from inside the cascade
 * (MATHILDA_INTEGRATE_PROFILE=1) over the integrands Integrate fails to close,
 * this stage is essentially the entire cost: on
 * `Sqrt[Sin[x]]/(Sin[x]^2 + 1)` it is 25.05 s of a 25.06 s cascade -- 100.0% --
 * while the ten stages ahead of it together cost 0.003 s. So the cascade is
 * NOT mis-ordered (the cheap exact methods are already first and are free);
 * one stage is unbounded.
 *
 * THE COST HAS TWO SOURCES, and the cascade pays both. try_derivdivides sends a
 * pseudo-elliptic integrand down the direct-only path but everything else to
 * integrate_derivdivides_full, so Pass 2's Eliminate/Solve search does run here
 * for most integrands. And Pass 1 is not cheap either: try_direct_kernel ->
 * integrate_in -> Integrate[qu, u] re-enters the WHOLE cascade once per candidate
 * kernel (the cascade switch runs at every depth), so a kernel list of n buys n
 * full descents, each able to fan out again. The gate in try_direct_kernel
 * addresses the second source structurally; this deadline bounds both, which is
 * why it is absolute -- set once at the outermost entry, so every nested
 * re-entry inherits the same cutoff. A per-call budget would be multiplied by
 * the fan-out, i.e. by the very thing being bounded.
 * ---------------------------------------------------------------------- */
static double dd_deadline = HUGE_VAL;

/* Monotonic seconds; CLOCK_MONOTONIC is immune to NTP steps, with the CPU clock
 * as fallback.  Mirrors simp_mono_seconds() / dt_wall_seconds(). */
static double dd_mono_seconds(void);
static bool   dd_expired(void);

/* Seconds left in the budget, for bounding a single sub-evaluation; 0 when
 * unbounded (so callers pass no TimeConstrained wrapper at all). */
static double dd_remaining(void) {
    if (dd_deadline >= HUGE_VAL) return 0.0;
    double left = dd_deadline - dd_mono_seconds();
    return left > 0.0 ? left : 0.001;      /* expired: ask for the minimum */
}

/* Evaluate `call` under what is left of the budget.  REQUIRED for the Eliminate/
 * Solve relation: it is ONE evaluation that can run for tens of seconds without
 * re-entering the cascade, so no between-stages check can bound it -- measured,
 * `TimeConstraint -> 3` let an integrand run 40 s because a single Eliminate call
 * consumed 37 of them. Consumes `call`; returns NULL when the budget is spent. */
static Expr* dd_eval_bounded(Expr* call) {
    double left = dd_remaining();
    if (left <= 0.0) return eval_take(call);          /* unbounded */
    Expr* guarded = expr_new_function(
        expr_new_symbol("TimeConstrained"),
        (Expr*[]){ call, expr_new_real(left), expr_new_symbol("$Aborted") }, 3);
    Expr* out = eval_take(guarded);                   /* consumes guarded */
    if (out && out->type == EXPR_SYMBOL
        && strcmp(out->data.symbol.name, "$Aborted") == 0) {
        expr_free(out);
        return NULL;
    }
    return out;
}

static double dd_mono_seconds(void) {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0)
        return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
    return (double)clock() / (double)CLOCKS_PER_SEC;
}

/* Armed ONCE per top-level evaluation, keyed on eval_toplevel_id().
 *
 * Not per call: the evaluator's fixed point re-enters builtin_integrate several
 * times for one user-level Integrate, so re-arming on entry multiplies the
 * budget by the number of passes -- measured, `TimeConstraint -> 3` let a case
 * run 30 s (ten passes) while `-> 1` stopped at 1.94 s (two). The epoch key is
 * the same one the fail-memo and DSolve`Separable's budget use. */
static uint64_t dd_budget_epoch = 0;
static bool     dd_budget_epoch_set = false;

void integrate_derivdivides_set_budget(double seconds) {
    uint64_t ep = eval_toplevel_id();
    if (dd_budget_epoch_set && ep == dd_budget_epoch) return;   /* same integral */
    dd_budget_epoch = ep;
    dd_budget_epoch_set = true;
    dd_deadline = (seconds > 0.0 && seconds < HUGE_VAL)
                      ? dd_mono_seconds() + seconds : HUGE_VAL;
}

static bool dd_expired(void) {
    return dd_deadline < HUGE_VAL && dd_mono_seconds() >= dd_deadline;
}

/* Canonical integrands already attempted in the current top-level descent.
 * Valid while dd_depth > 0; freed by the outermost frame on exit. */
static ExprVec dd_seen = { NULL, 0, 0 };

/* Canonical key for (f, x): f with the integration variable renamed to a
 * fixed sentinel, so integrands that differ only in their (gensym'd) variable
 * name compare equal.  Returns an owned Expr the caller must free. */
static Expr* dd_canon_key(const Expr* f, const Expr* x) {
    Expr* sentinel = expr_new_symbol("Integrate`DerivativeDivides`$var");
    Expr* key = replace_one(expr_copy((Expr*)f), x, sentinel); /* consumes copy */
    expr_free(sentinel);
    return key;
}

/* True iff `key` was already attempted in this descent. */
static bool dd_seen_contains(const Expr* key) {
    if (!key) return false;
    for (size_t i = 0; i < dd_seen.n; i++)
        if (expr_eq(dd_seen.items[i], (Expr*)key)) return true;
    return false;
}

/* Per-call fresh-symbol counter so nested substitutions never collide. */
static unsigned long dd_sym_counter = 0;

static void fresh_uy(Expr** u, Expr** y) {
    char nu[80], ny[80];
    snprintf(nu, sizeof(nu), "Integrate`DerivativeDivides`u$%lu", dd_sym_counter);
    snprintf(ny, sizeof(ny), "Integrate`DerivativeDivides`y$%lu", dd_sym_counter);
    dd_sym_counter++;
    *u = expr_new_symbol(nu);
    *y = expr_new_symbol(ny);
}

/* ---------------------------------------------------------------------- */
/* Strategy 1: direct quotient                                            */
/* ---------------------------------------------------------------------- */

/* If `w` is a pure unit root of the integration variable, w = x^(1/n)
 * (n >= 2), return n; else 0.  For such a kernel the substitution u = x^(1/n)
 * inverts cleanly as x = u^n, which lets the direct-quotient strategy replace
 * the residual x that a bare w -> u leaves behind (e.g. Sqrt[x + Sqrt[x]] with
 * kernel Sqrt[x] gives 2 u Sqrt[u + x]; closing x -> u^2 yields the clean
 * 2 u Sqrt[u + u^2], keeping the antiderivative in the integrand's own radical
 * Sqrt[x + Sqrt[x]] rather than introducing x^(1/4) generators). */
static int64_t kernel_root_degree(const Expr* w, const Expr* x) {
    if (!w || w->type != EXPR_FUNCTION) return 0;
    if (!head_is((Expr*)w, SYM_Power) || w->data.function.arg_count != 2) return 0;
    if (!expr_eq(w->data.function.args[0], (Expr*)x)) return 0;
    const Expr* e = w->data.function.args[1];
    if (e->type == EXPR_FUNCTION && head_is((Expr*)e, SYM_Rational)
        && e->data.function.arg_count == 2) {
        const Expr* a = e->data.function.args[0];
        const Expr* b = e->data.function.args[1];
        if (a->type == EXPR_INTEGER && b->type == EXPR_INTEGER
            && a->data.integer == 1 && b->data.integer >= 2)
            return b->data.integer;
    }
    return 0;
}

/* Does `e` mention the symbol `u`?  A plain structural test -- expr_free_of()
 * in this file is the derivative-based one, which is far more than needed here
 * and would evaluate D[] on every node of the walk. */
static bool dd_contains_sym(const Expr* e, const Expr* u) {
    if (!e || !u || u->type != EXPR_SYMBOL) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == u->data.symbol.name;
    if (e->type != EXPR_FUNCTION) return false;
    if (dd_contains_sym(e->data.function.head, u)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (dd_contains_sym(e->data.function.args[i], u)) return true;
    return false;
}

static bool dd_polynomial_in(const Expr* g, const Expr* u) {
    Expr* args[2] = { expr_copy((Expr*)g), expr_copy((Expr*)u) };
    Expr* call = internal_polynomialq(args, 2);
    Expr* val  = evaluate(call);
    expr_free(call);
    bool ok = val && val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_True;
    expr_free(val);
    return ok;
}

static bool dd_rational_in(const Expr* g, const Expr* u) {
    Expr* tog = eval_take(mk_fn1("Together", expr_copy((Expr*)g)));
    if (!tog) return false;
    Expr* num = eval_take(mk_fn1("Numerator", expr_copy(tog)));
    Expr* den = eval_take(mk_fn1("Denominator", tog));          /* consumes tog */
    bool ok = num && den && dd_polynomial_in(num, u) && dd_polynomial_in(den, u);
    expr_free(num); expr_free(den);
    return ok;
}

/* Every sub-structure containing u must be rational, or a rational power of a
 * polynomial of degree <= 2 in u.  Walks the tree; no evaluation except the
 * degree test, which only runs on a radicand. */
static bool dd_algebraic_rationalisable(const Expr* g, const Expr* u) {
    if (!g) return true;
    if (!dd_contains_sym(g, u)) return true;                 /* free of u: fine */
    if (g->type == EXPR_SYMBOL) return true;                 /* u itself */
    if (g->type != EXPR_FUNCTION) return true;
    const Expr* h = g->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    const char* hn = h->data.symbol.name;
    if (hn == SYM_Plus || hn == SYM_Times) {
        for (size_t i = 0; i < g->data.function.arg_count; i++)
            if (!dd_algebraic_rationalisable(g->data.function.args[i], u)) return false;
        return true;
    }
    /* Exp / Log of something already in the class: the reduced integrand is then
     * rational in the Risch monomial tower (u, E^p(u), Log p(u)), which the
     * transcendental Risch stage decides algorithmically -- the same standard as
     * "convertible to rational", one tower level up. Without this the gate loses
     * genuine closes: O079/O080 of the overintegration corpus,
     * E^ArcCsc[x] (... Log[...] ...), reduce to exactly this shape and were
     * closed by RischTranscendental before the gate existed. */
    if ((hn == SYM_Exp || hn == SYM_Log) && g->data.function.arg_count == 1)
        return dd_algebraic_rationalisable(g->data.function.args[0], u);
    if (hn == SYM_Power && g->data.function.arg_count == 2) {
        /* E^anything-in-the-class is the Exp case written as a Power. */
        const Expr* b0 = g->data.function.args[0];
        if (b0->type == EXPR_SYMBOL && b0->data.symbol.name == SYM_E)
            return dd_algebraic_rationalisable(g->data.function.args[1], u);
        const Expr* base = g->data.function.args[0];
        const Expr* ex   = g->data.function.args[1];
        if (ex->type == EXPR_INTEGER)                        /* rational so far */
            return dd_algebraic_rationalisable(base, u);
        /* fractional exponent: the radicand must be a polynomial of degree <= 2 */
        if (!dd_polynomial_in(base, u)) return false;
        Expr* deg = eval_take(mk_fn2("Exponent", expr_copy((Expr*)base),
                                     expr_copy((Expr*)u)));
        bool small = deg && deg->type == EXPR_INTEGER && deg->data.integer <= 2;
        expr_free(deg);
        return small;
    }
    return false;                                            /* transcendental in u */
}

MATHILDA_MAYBE_UNUSED static bool dd_decidable_class(const Expr* g, const Expr* u) {
    if (dd_polynomial_in(g, u)) return true;
    if (dd_rational_in(g, u))   return true;
    return dd_algebraic_rationalisable(g, u);
}

/* Try kernel `w` via  q = Cancel[Together[f / D[w, x]]]; q /. w -> u.
 * Returns a verified antiderivative or NULL.  Borrows everything. */
static Expr* try_direct_kernel(const Expr* f, const Expr* x, const Expr* w) {
    Expr* du = deriv_dx(w, x);
    if (!du || is_lit_zero(du)) { if (du) expr_free(du); return NULL; }

    /* q = Cancel[Together[f * du^-1]] */
    Expr* quot = mk_fn2("Times", expr_copy((Expr*)f),
                        mk_fn2("Power", du, mk_int(-1)));   /* consumes du */
    Expr* q = cancel_together(quot);
    if (!q) return NULL;

    Expr* u; Expr* y; fresh_uy(&u, &y);
    Expr* qu = replace_one(q, w, u);                        /* consumes q */

    /* A bare w -> u can leave the variable behind when w is a root of x (the
     * radicand still mentions x, e.g. Sqrt[x + Sqrt[x]] -> 2 u Sqrt[u + x]).
     * For a unit root w = x^(1/n) the inverse x = u^n closes the substitution
     * without introducing foreign radical generators; differentiate-back still
     * gates correctness. */
    if (qu && !expr_free_of(qu, x)) {
        int64_t n = kernel_root_degree(w, x);
        if (n >= 2) {
            Expr* upow = mk_fn2("Power", expr_copy(u), mk_int(n));
            qu = replace_one(qu, x, upow);                  /* consumes qu */
            expr_free(upow);
            /* PowerExpand collapses the (u^n)^(p/n) powers the x -> u^n
             * substitution creates (e.g. from x^(5/2) -> (u^2)^(5/2)) down to
             * u^p, under u = x^(1/n) >= 0 (the integration-positivity
             * convention).  It leaves genuine radicals of polynomials in u --
             * Sqrt[u^2 + u] etc. -- untouched, so the reduced integrand stays a
             * clean function of the integrand's own radical.  differentiate-back
             * against the original f re-verifies, so a branch-crossing expansion
             * is rejected rather than trusted. */
            if (qu) qu = eval_take(mk_fn1("PowerExpand", qu));
        }
    }

    Expr* result = NULL;
    if (qu && expr_free_of(qu, x)) {
        Expr* G = integrate_in(qu, u);
        if (G) {
            Expr* r = replace_one(G, u, w);                 /* consumes G */
            if (r && differentiates_back(r, f, x)) {
                result = r;
            } else if (r) {
                expr_free(r);
            }
        }
    }
    if (qu) expr_free(qu);
    expr_free(u); expr_free(y);
    return result;
}


/* ------------------------------------------------------------------------
 * Is `g` in a class the cascade's cheap deterministic stages can DECIDE?
 *
 * The recursion in Pass 1 is speculative: for every candidate kernel w it forms
 * the reduced integrand qu and calls Integrate[qu, u], which re-enters the whole
 * cascade. When qu is rational -- or one of the concrete algebraic shapes that
 * rationalise ALGORITHMICALLY -- that descent is cheap and decisive, and it is
 * the point of the substitution: the whole method exists to turn an integrand
 * into a rational one. When qu is anything else, the descent is a full cascade
 * search on an integrand no simpler than the one we started from, and the fan-out
 * (one descent per kernel, each able to fan out again) is what makes this stage
 * 47% of the cost of every integral the cascade fails to close.
 *
 * So: recurse when qu is polynomial or rational in u, or when its only
 * non-rational structure is a rational power of a polynomial of degree <= 2 in u
 * -- the linear-radical, quadratic-radical, linear-ratio-radical and
 * binomial-differential forms that stages 03-06 rationalise by construction.
 * A transcendental kernel in u (Log, Exp, a trig head) or a radical of a cubic or
 * higher polynomial is NOT decidable here and does not earn a descent.
 *
 * This is a gate on the RECURSION, not on the answer: the differentiate-back
 * check still guards every result, so a declined descent can only cost a kernel
 * that would have closed, never produce a wrong one.
 * ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Strategy 2: Eliminate / Solve with branch selection                    */
/* ---------------------------------------------------------------------- */

/* Try kernel `w` via the differential Eliminate/Solve relation.  Iterates the
 * Solve branches and returns the first that differentiates back to f, or NULL.
 * Borrows everything.  May emit Eliminate ::ifun / ::alg diagnostics. */
static Expr* try_eliminate_kernel(const Expr* f, const Expr* x, const Expr* w) {
    Expr* u; Expr* y; fresh_uy(&u, &y);
    Expr* dtx = mk_fn1("Dt", expr_copy((Expr*)x));
    Expr* dtu = mk_fn1("Dt", expr_copy(u));
    Expr* dty = mk_fn1("Dt", expr_copy(y));

    /* eqns = { Dt[y] == f Dt[x], u == w, Dt[u == w] } */
    Expr* eq1 = mk_fn2("Equal", expr_copy(dty),
                       mk_fn2("Times", expr_copy((Expr*)f), expr_copy(dtx)));
    Expr* eq2 = mk_fn2("Equal", expr_copy(u), expr_copy((Expr*)w));
    Expr* eq3 = mk_fn1("Dt", mk_fn2("Equal", expr_copy(u), expr_copy((Expr*)w)));
    Expr* eqns = expr_new_function(expr_new_symbol(SYM_List),
                                   (Expr*[]){ eq1, eq2, eq3 }, 3);
    /* vars = { x, Dt[x] } */
    Expr* vars = expr_new_function(expr_new_symbol(SYM_List),
                                   (Expr*[]){ expr_copy((Expr*)x), expr_copy(dtx) }, 2);

    /* Eliminate runs as a private sub-step here; mute its advisory ::ifun /
     * ::alg / ::nlin diagnostics so the user never sees branch caveats for an
     * elimination they did not request. */
    int saved_quiet = eliminate_suppress_messages;
    eliminate_suppress_messages = 1;
    /* Bounded: a single Eliminate on a radical relation can run for tens of
     * seconds without re-entering the cascade, so this is the one place a
     * between-stages deadline check cannot reach. */
    Expr* rel = dd_eval_bounded(mk_fn2("Eliminate", eqns, vars)); /* consumes eqns, vars */
    eliminate_suppress_messages = saved_quiet;

    Expr* result = NULL;
    /* Eliminate must have produced an equation, not bailed back unevaluated. */
    if (rel && head_is(rel, SYM_Equal)) {
        Expr* sol = dd_eval_bounded(mk_fn2("Solve", expr_copy(rel), expr_copy(dty)));
        if (sol) {
            /* gs = Cancel[ PowerExpand[ Factor //@ (Dt[y] /. sol) ] / Dt[u] ] */
            Expr* dvals = eval_take(internal_replace_all(
                (Expr*[]){ expr_copy(dty), expr_copy(sol) }, 2));
            Expr* fac = dvals ? eval_take(mk_fn2("MapAll",
                                expr_new_symbol(SYM_Factor), dvals)) : NULL;
            Expr* pe  = fac ? eval_take(mk_fn1("PowerExpand", fac)) : NULL;
            Expr* gs  = pe ? cancel_together(mk_fn2("Times", pe,
                                mk_fn2("Power", expr_copy(dtu), mk_int(-1)))) : NULL;

            if (gs) {
                /* gs is a List (one entry per branch) or a single expr. */
                size_t nb = (head_is(gs, SYM_List)) ? gs->data.function.arg_count : 1;
                for (size_t b = 0; b < nb && !result; b++) {
                    Expr* g = head_is(gs, SYM_List)
                            ? gs->data.function.args[b] : gs;
                    if (!expr_free_of(g, x))   continue;
                    if (!expr_free_of(g, dtx)) continue;
                    if (!expr_free_of(g, dtu)) continue;
                    if (dd_expired()) break;
                    Expr* G = integrate_in(g, u);
                    if (!G) continue;
                    Expr* r = replace_one(G, u, w);          /* consumes G */
                    if (r && differentiates_back(r, f, x)) result = r;
                    else if (r) expr_free(r);
                }
                expr_free(gs);
            }
        }
        if (sol) expr_free(sol);
    }
    if (rel) expr_free(rel);
    expr_free(dtx); expr_free(dtu); expr_free(dty);
    expr_free(u); expr_free(y);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Core driver                                                            */
/* ---------------------------------------------------------------------- */

/* `use_eliminate` selects whether the (diagnostic-emitting, more thorough)
 * Eliminate/Solve strategy is attempted in addition to the direct one.
 *
 * `forced_kernel`, when non-NULL, pins the substitution: instead of collecting
 * and trialing every x-dependent subexpression, only that single kernel u(x) is
 * attempted.  Used by the explicit
 * `Method -> {"DerivativeDivides", "Substitution" -> u}` option. */
static Expr* dd_core(Expr* f, Expr* x, bool use_eliminate,
                     const Expr* forced_kernel) {
    if (x->type != EXPR_SYMBOL) return NULL;
    if (expr_free_of(f, x))      return NULL;   /* nothing to integrate in x */
    if (dd_depth >= DD_MAX_DEPTH) return NULL;
    if (dd_expired())            return NULL;   /* TimeConstraint spent */
    /* Symbolic-exponent monomials (x^k, x^(k-1), ...) send the Together/Cancel
     * quotient normalisation into a PolynomialGCD pseudo-remainder blow-up and
     * are never a productive substitution kernel -- decline them up front. */
    if (has_symbolic_power_of(f, x)) return NULL;

    bool outermost = (dd_depth == 0);

    /* Loop guard: short-circuit any integrand (modulo the integration
     * variable) we have already attempted in this descent.  Breaks circular
     * substitution chains and collapses overlapping subproblems so the
     * recursion cannot fan out exponentially. */
    Expr* key = dd_canon_key(f, x);
    if (dd_seen_contains(key)) { expr_free(key); return NULL; }
    vec_add_unique(&dd_seen, key);              /* copies key into the memo */
    expr_free(key);

    ExprVec kernels; vec_init(&kernels);
    Expr* result = NULL;
    if (forced_kernel) {
        /* A user-pinned substitution must depend on x; otherwise it cannot
         * reduce the integral and the integral stays unevaluated (strict, no
         * fallback to the automatic kernel search). */
        if (expr_free_of(forced_kernel, x)) { vec_free(&kernels); goto done; }
        vec_add_unique(&kernels, forced_kernel);
    } else {
        collect_kernels(f, x, f, &kernels);
    }
    if (kernels.n == 0) { vec_free(&kernels); goto done; }
    sort_kernels(&kernels);

    dd_depth++;

    /* Pass 1: the cheap, quiet, branch-correct direct quotient. */
    for (size_t i = 0; i < kernels.n && !result && !dd_expired(); i++)
        result = try_direct_kernel(f, x, kernels.items[i]);

    /* Pass 2: the thorough Eliminate/Solve search.  Heavyweight (~0.1-1s per
     * kernel), so it runs only on the outermost integrand; reduced sub-integrals
     * are finished by the direct strategy above and the rest of the cascade.
     * See the recursion-guard note above for why. */
    if (!result && use_eliminate && outermost)
        for (size_t i = 0; i < kernels.n && !result && !dd_expired(); i++)
            result = try_eliminate_kernel(f, x, kernels.items[i]);

    dd_depth--;
    vec_free(&kernels);

done:
    /* The outermost frame owns the memo for the whole descent; reset it so the
     * next independent top-level integral starts clean. */
    if (outermost) vec_free(&dd_seen);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Public entry points                                                    */
/* ---------------------------------------------------------------------- */

Expr* integrate_derivdivides_try(Expr* f, Expr* x) {
    /* Automatic cascade: direct strategy only (fast, no diagnostics). */
    return dd_core(f, x, false, NULL);
}

Expr* integrate_derivdivides_full(Expr* f, Expr* x) {
    /* Explicit Method -> "DerivativeDivides": direct strategy, then the
     * thorough Eliminate/Solve branch-search. */
    return dd_core(f, x, true, NULL);
}

Expr* integrate_derivdivides_with_sub(Expr* f, Expr* x, Expr* sub) {
    /* Explicit Method -> {"DerivativeDivides", "Substitution" -> sub}: trial
     * only the user-pinned kernel `sub`, running both the direct and the
     * thorough Eliminate/Solve strategies. */
    return dd_core(f, x, true, sub);
}

Expr* builtin_integrate_derivdivides(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    /* Explicit method: direct strategy, then the thorough Eliminate search. */
    return dd_core(f, x, true, NULL);
}

void integrate_derivdivides_init(void) {
    symtab_add_builtin("Integrate`DerivativeDivides", builtin_integrate_derivdivides);
    symtab_get_def("Integrate`DerivativeDivides")->attributes |=
        ATTR_PROTECTED;
    symtab_set_docstring("Integrate`DerivativeDivides",
        "Integrate`DerivativeDivides[f, x] integrates by substitution: it finds a\n"
        "kernel u(x) whose derivative divides f, reduces to Integrate[h[u], u],\n"
        "back-substitutes u -> u(x), and verifies the result by differentiation.\n"
        "Tries a direct quotient first, then an Eliminate/Solve search that can\n"
        "select the correct branch for radical substitutions. Strict: returns\n"
        "unevaluated when no substitution closes the integral.");
}
