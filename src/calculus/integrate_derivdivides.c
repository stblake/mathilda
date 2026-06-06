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
 *             if Simplify[D[r, x] - f] === 0:  return r
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
 *                 if Simplify[D[r, x] - f] === 0:  return r
 *
 * The Simplify[D[r, x] - f] === 0 gate is unconditional and is what selects
 * the correct inverse-function branch among Solve's +- solutions and the
 * branch assumptions PowerExpand bakes in.
 *
 * Memory: builtins take ownership of `res`; helpers below own every Expr they
 * construct and free intermediates explicitly.  `eval_take` consumes its
 * argument.  We never expr_free(res).
 */

#include "integrate_derivdivides.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include "internal.h"
#include "sym_intern.h"
#include "sym_names.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
 * A cheap PossibleZeroQ pre-screen rejects the (common) clearly-non-zero
 * candidate before paying for Simplify: PossibleZeroQ uses numeric sampling
 * and so settles most wrong branches without any symbolic work.  Only when it
 * reports a likely zero do we confirm rigorously with Simplify, then a
 * Cancel[Together[Expand[.]]] fallback for the forms Simplify leaves
 * non-canonical.  Borrows r and f. */
static bool differentiates_back(const Expr* r, const Expr* f, const Expr* x) {
    Expr* dr = deriv_dx(r, x);
    if (!dr) return false;
    Expr* diff = mk_fn2("Plus", dr,
                        mk_fn2("Times", mk_int(-1), expr_copy((Expr*)f)));

    /* Cheap pre-screen: bail unless PossibleZeroQ[diff] === True. */
    Expr* pz = eval_take(mk_fn1("PossibleZeroQ", expr_copy(diff)));
    bool possible = pz && pz->type == EXPR_SYMBOL && pz->data.symbol == SYM_True;
    if (pz) expr_free(pz);
    if (!possible) { expr_free(diff); return false; }

    /* Confirm: Simplify[diff], then Cancel[Together[Expand[diff]]]. */
    Expr* s = eval_take(mk_fn1("Simplify", expr_copy(diff)));
    bool zero = is_lit_zero(s);
    if (s) expr_free(s);
    if (!zero) {
        Expr* ex = eval_take(internal_expand((Expr*[]){ expr_copy(diff) }, 1));
        Expr* ct = cancel_together(ex);
        zero = is_lit_zero(ct);
        if (ct) expr_free(ct);
    }
    expr_free(diff);
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

#define DD_MAX_DEPTH 8
static int dd_depth = 0;

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
    Expr* eqns = expr_new_function(expr_new_symbol("List"),
                                   (Expr*[]){ eq1, eq2, eq3 }, 3);
    /* vars = { x, Dt[x] } */
    Expr* vars = expr_new_function(expr_new_symbol("List"),
                                   (Expr*[]){ expr_copy((Expr*)x), expr_copy(dtx) }, 2);

    Expr* rel = eval_take(mk_fn2("Eliminate", eqns, vars));   /* consumes eqns, vars */

    Expr* result = NULL;
    /* Eliminate must have produced an equation, not bailed back unevaluated. */
    if (rel && head_is(rel, SYM_Equal)) {
        Expr* sol = eval_take(mk_fn2("Solve", expr_copy(rel), expr_copy(dty)));
        if (sol) {
            /* gs = Cancel[ PowerExpand[ Factor //@ (Dt[y] /. sol) ] / Dt[u] ] */
            Expr* dvals = eval_take(internal_replace_all(
                (Expr*[]){ expr_copy(dty), expr_copy(sol) }, 2));
            Expr* fac = dvals ? eval_take(mk_fn2("MapAll",
                                expr_new_symbol("Factor"), dvals)) : NULL;
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
 * Eliminate/Solve strategy is attempted in addition to the direct one. */
static Expr* dd_core(Expr* f, Expr* x, bool use_eliminate) {
    if (x->type != EXPR_SYMBOL) return NULL;
    if (expr_free_of(f, x))      return NULL;   /* nothing to integrate in x */
    if (dd_depth >= DD_MAX_DEPTH) return NULL;

    ExprVec kernels; vec_init(&kernels);
    collect_kernels(f, x, f, &kernels);
    if (kernels.n == 0) { vec_free(&kernels); return NULL; }
    sort_kernels(&kernels);

    dd_depth++;
    Expr* result = NULL;

    /* Pass 1: the cheap, quiet, branch-correct direct quotient. */
    for (size_t i = 0; i < kernels.n && !result; i++)
        result = try_direct_kernel(f, x, kernels.items[i]);

    /* Pass 2: the thorough Eliminate/Solve search (explicit method only). */
    if (!result && use_eliminate)
        for (size_t i = 0; i < kernels.n && !result; i++)
            result = try_eliminate_kernel(f, x, kernels.items[i]);

    dd_depth--;
    vec_free(&kernels);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Public entry points                                                    */
/* ---------------------------------------------------------------------- */

Expr* integrate_derivdivides_try(Expr* f, Expr* x) {
    /* Automatic cascade: direct strategy only (fast, no diagnostics). */
    return dd_core(f, x, false);
}

Expr* builtin_integrate_derivdivides(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    /* Explicit method: direct strategy, then the thorough Eliminate search. */
    return dd_core(f, x, true);
}

void integrate_derivdivides_init(void) {
    symtab_add_builtin("Integrate`DerivativeDivides", builtin_integrate_derivdivides);
    symtab_get_def("Integrate`DerivativeDivides")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`DerivativeDivides",
        "Integrate`DerivativeDivides[f, x] integrates by substitution: it finds a\n"
        "kernel u(x) whose derivative divides f, reduces to Integrate[h[u], u],\n"
        "back-substitutes u -> u(x), and verifies the result by differentiation.\n"
        "Tries a direct quotient first, then an Eliminate/Solve search that can\n"
        "select the correct branch for radical substitutions. Strict: returns\n"
        "unevaluated when no substitution closes the integral.");
}
