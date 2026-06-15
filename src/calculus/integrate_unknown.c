/* integrate_unknown.c
 *
 * Integration of expressions containing undefined (unknown) functions and
 * their derivatives.  Implements the "polynomial part" of Kelly Roach's
 * §1.7 algorithm ("Undefined Functions", Roach 1992) via an equivalent,
 * implementation-friendly integration-by-parts reformulation.  The
 * fraction part and resultant-based log part are layered on in later
 * phases.
 *
 * ---------------------------------------------------------------------
 * Algorithm (polynomial part)
 * ---------------------------------------------------------------------
 * The integrand I is viewed as a polynomial in the "generator atoms"
 *
 *     theta^(n) = u^(n)(g)      (u undefined, g an expression in x),
 *
 * with coefficients that are ordinary (rational / elementary) functions
 * of x.  Pick the atom A of highest derivative order appearing in I and
 * let B be its antiderivative atom (A = dB/dx).  If I is linear in A,
 * write I = C*A + D with C, D free of A.  Then by parts
 *
 *     Integrate[I, x] = W  -  Integrate[ dW/dx - C*A - D, x ]
 *
 * where W is the antiderivative of C with respect to B (a polynomial
 * integration in B).  dW/dx is computed with the native D builtin, which
 * already differentiates undefined functions via the chain rule, so the
 * new integrand is strictly free of A and the recursion terminates.
 *
 * To reuse the standard polynomial machinery (Coefficient,
 * IntegratePolynomial) the atoms are first substituted for fresh
 * symbols, the algebra is performed in that polynomial ring, and the
 * result is substituted back to real atoms before the D step.
 *
 * Memory: builtins take ownership of `res`; helpers below own every Expr
 * they construct and free intermediates explicitly.  `eval_take` and
 * `subst` consume their primary Expr argument.
 */

#include "integrate_unknown.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include "internal.h"
#include "sym_intern.h"
#include "sym_names.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------- */
/* Small builders / evaluation helpers                                    */
/* ---------------------------------------------------------------------- */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }

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

/* True if `f` contains no subexpression structurally equal to `x`. */
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

/* True iff `name` (an interned symbol string) currently has no builtin,
 * down-values, or own-values: i.e. it behaves as a free function head. */
static bool symbol_is_free(const char* name) {
    SymbolDef* d = symtab_lookup(name);
    if (!d) return true;
    if (d->builtin_func) return false;
    if (d->down_values) return false;
    if (d->own_values) return false;
    return true;
}

/* ---------------------------------------------------------------------- */
/* Undefined-function atom recognition                                    */
/* ---------------------------------------------------------------------- */

/* If `e` is an undefined-function atom that depends on `x` -- either
 *   f[g]                       (order 0), or
 *   Derivative[n][f][g]        (order n >= 1),
 * with f a free symbol and a single argument g (g depends on x) -- return
 * true and fill *order, *func (borrowed f), *arg (borrowed g). */
static bool classify_atom(const Expr* e, const Expr* x,
                          int* order, const Expr** func, const Expr** arg) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* head = e->data.function.head;

    /* Case (a): f[g]. */
    if (head->type == EXPR_SYMBOL) {
        if (e->data.function.arg_count != 1) return false;
        if (!symbol_is_free(head->data.symbol)) return false;
        const Expr* g = e->data.function.args[0];
        if (expr_free_of(g, x)) return false;
        *order = 0; *func = head; *arg = g;
        return true;
    }

    /* Case (b): Derivative[n][f][g].  head == Derivative[n][f]. */
    if (head->type == EXPR_FUNCTION) {
        if (head->data.function.arg_count != 1) return false;        /* [f] */
        const Expr* dop = head->data.function.head;                  /* Derivative[n] */
        if (dop->type != EXPR_FUNCTION) return false;
        if (!head_is(dop, SYM_Derivative)) return false;
        if (dop->data.function.arg_count != 1) return false;         /* single n */
        const Expr* nexp = dop->data.function.args[0];
        if (nexp->type != EXPR_INTEGER || nexp->data.integer < 1) return false;
        const Expr* f = head->data.function.args[0];
        if (f->type != EXPR_SYMBOL || !symbol_is_free(f->data.symbol)) return false;
        if (e->data.function.arg_count != 1) return false;
        const Expr* g = e->data.function.args[0];
        if (expr_free_of(g, x)) return false;
        *order = (int)nexp->data.integer; *func = f; *arg = g;
        return true;
    }
    return false;
}

/* Build the atom u^(order)(arg): f[arg] for order 0, else
 * Derivative[order][f][arg]. */
static Expr* make_atom(const Expr* func, int order, const Expr* arg) {
    if (order == 0) {
        Expr* a[1] = { expr_copy((Expr*)arg) };
        return expr_new_function(expr_copy((Expr*)func), a, 1);
    }
    Expr* nidx[1] = { mk_int(order) };
    Expr* dop = expr_new_function(expr_new_symbol(SYM_Derivative), nidx, 1);
    Expr* df[1] = { expr_copy((Expr*)func) };
    Expr* dop_f = expr_new_function(dop, df, 1);
    Expr* g[1] = { expr_copy((Expr*)arg) };
    return expr_new_function(dop_f, g, 1);
}

/* ---------------------------------------------------------------------- */
/* Atom map: unique atoms <-> fresh polynomial-ring symbols               */
/* ---------------------------------------------------------------------- */

typedef struct {
    Expr** atoms;   /* owned atom expressions                             */
    Expr** syms;    /* owned fresh symbols (one per atom)                 */
    int*   orders;  /* derivative order of each atom                      */
    size_t n;
    size_t cap;
} AtomMap;

static unsigned long iu_sym_counter = 0;

static void map_init(AtomMap* m) {
    m->atoms = NULL; m->syms = NULL; m->orders = NULL; m->n = 0; m->cap = 0;
}

static void map_free(AtomMap* m) {
    for (size_t i = 0; i < m->n; i++) {
        expr_free(m->atoms[i]);
        expr_free(m->syms[i]);
    }
    free(m->atoms); free(m->syms); free(m->orders);
    map_init(m);
}

/* Index of `atom` in the map, or -1. */
static long map_find(const AtomMap* m, const Expr* atom) {
    for (size_t i = 0; i < m->n; i++) {
        if (expr_eq(m->atoms[i], (Expr*)atom)) return (long)i;
    }
    return -1;
}

/* Add `atom` (copied) with a fresh symbol if not already present. */
static void map_add(AtomMap* m, const Expr* atom, int order) {
    if (map_find(m, atom) >= 0) return;
    if (m->n == m->cap) {
        m->cap = m->cap ? m->cap * 2 : 8;
        m->atoms  = realloc(m->atoms,  m->cap * sizeof(Expr*));
        m->syms   = realloc(m->syms,   m->cap * sizeof(Expr*));
        m->orders = realloc(m->orders, m->cap * sizeof(int));
    }
    char name[64];
    snprintf(name, sizeof(name), "Integrate`u$%lu", iu_sym_counter++);
    m->atoms[m->n]  = expr_copy((Expr*)atom);
    m->syms[m->n]   = expr_new_symbol(name);
    m->orders[m->n] = order;
    m->n++;
}

/* Walk `e`, recording every undefined-function atom that depends on x.
 * Atoms are opaque: we do not recurse into a recognised atom. */
static void collect_atoms(const Expr* e, const Expr* x, AtomMap* m) {
    int order; const Expr* func; const Expr* arg;
    if (classify_atom(e, x, &order, &func, &arg)) {
        map_add(m, e, order);
        return;
    }
    if (e->type == EXPR_FUNCTION) {
        collect_atoms(e->data.function.head, x, m);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            collect_atoms(e->data.function.args[i], x, m);
        }
    }
}

/* Ensure that for every atom of order n in the map, its whole
 * antiderivative chain (orders n-1, ..., 0 with the same func/arg) also
 * has a fresh symbol, so by-parts can always name the antiderivative
 * atom B even when B itself does not appear in the integrand. */
static void map_complete_chains(AtomMap* m, const Expr* x) {
    /* Snapshot the current atoms; new ones are appended as we go. */
    for (size_t i = 0; i < m->n; i++) {
        int order; const Expr* func; const Expr* arg;
        if (!classify_atom(m->atoms[i], x, &order, &func, &arg)) continue;
        for (int k = order - 1; k >= 0; k--) {
            Expr* lower = make_atom(func, k, arg);
            map_add(m, lower, k);
            expr_free(lower);
        }
    }
}

/* True iff `e` contains any undefined-function atom of x (order >= 0). */
static bool contains_unknown_atom(const Expr* e, const Expr* x) {
    AtomMap m; map_init(&m);
    collect_atoms(e, x, &m);
    bool has = (m.n > 0);
    map_free(&m);
    return has;
}

/* ---------------------------------------------------------------------- */
/* Transcendental (Log[eta]) generators                                   */
/* ---------------------------------------------------------------------- */

/* A simple growable vector of owned Expr*. */
typedef struct { Expr** items; size_t n, cap; } ExprVec;

static void evec_init(ExprVec* v) { v->items = NULL; v->n = 0; v->cap = 0; }
static void evec_free(ExprVec* v) {
    for (size_t i = 0; i < v->n; i++) expr_free(v->items[i]);
    free(v->items); evec_init(v);
}
static void evec_add_unique(ExprVec* v, const Expr* e) {
    for (size_t i = 0; i < v->n; i++) if (expr_eq(v->items[i], (Expr*)e)) return;
    if (v->n == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 4;
        v->items = realloc(v->items, v->cap * sizeof(Expr*));
    }
    v->items[v->n++] = expr_copy((Expr*)e);
}

/* Collect every `Log[eta]` subexpression whose argument depends on an
 * undefined function of x.  Such logarithms are transcendental over the
 * unknown-function field and act as additional generators (Roach §1.6/§1.7);
 * ordinary `Log[x]` (eta free of unknowns) is handled by the standard
 * integrators and is not collected.  Outer logs are added before inner ones. */
static void collect_log_gens(const Expr* e, const Expr* x, ExprVec* out) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (head_is(e, SYM_Log) && e->data.function.arg_count == 1 &&
        contains_unknown_atom(e->data.function.args[0], x)) {
        evec_add_unique(out, e);
    }
    collect_log_gens(e->data.function.head, x, out);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        collect_log_gens(e->data.function.args[i], x, out);
    }
}

static bool has_log_gen(const Expr* f, const Expr* x) {
    ExprVec v; evec_init(&v);
    collect_log_gens(f, x, &v);
    bool has = (v.n > 0);
    evec_free(&v);
    return has;
}

/* ---------------------------------------------------------------------- */
/* Substitution between atoms and fresh symbols                           */
/* ---------------------------------------------------------------------- */

/* ReplaceAll[expr, {from_i -> to_i}], consuming `expr`. */
static Expr* subst(Expr* expr, Expr** from, Expr** to, size_t n) {
    Expr** rules = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        rules[i] = mk_fn2("Rule", expr_copy(from[i]), expr_copy(to[i]));
    }
    Expr* rulelist = expr_new_function(expr_new_symbol(SYM_List), rules, n);
    free(rules);
    Expr* call = internal_replace_all((Expr*[]){ expr, rulelist }, 2);
    return eval_take(call);
}

static Expr* atoms_to_syms(Expr* expr, const AtomMap* m) {
    if (m->n == 0) return expr;
    return subst(expr, m->atoms, m->syms, m->n);
}

static Expr* syms_to_atoms(Expr* expr, const AtomMap* m) {
    if (m->n == 0) return expr;
    return subst(expr, m->syms, m->atoms, m->n);
}

/* ---------------------------------------------------------------------- */
/* Polynomial-ring helpers (operating on fresh symbols)                   */
/* ---------------------------------------------------------------------- */

/* Coefficient[expr, var, k].  Borrows; returns owned result. */
static Expr* coefficient_k(const Expr* expr, const Expr* var, int k) {
    Expr* call = internal_coefficient(
        (Expr*[]){ expr_copy((Expr*)expr), expr_copy((Expr*)var), mk_int(k) }, 3);
    return eval_take(call);
}

/* True if `e` contains any unevaluated Integrate[...] call. */
static bool contains_unintegrated(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Integrate)) return true;
    if (contains_unintegrated(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_unintegrated(e->data.function.args[i])) return true;
    }
    return false;
}

/* Antiderivative of `integrand` with respect to the ring symbol `var`,
 * computed with the full (rational + log + arctan) integrator.  Because
 * `var` is an ordinary symbol and `integrand` is rational in the
 * generator symbols, this dispatches to the rational integrator -- it
 * does NOT recurse back into the unknown-function engine.  Returns NULL
 * if the integral cannot be closed (any Integrate[...] left over).  This
 * single step subsumes Roach's polynomial, fraction, and log parts:
 *   1/tB     -> Log[tB],   1/tB^2 -> -1/tB,   -a/(a^2+tB^2) -> -ArcTan[tB/a]. */
static Expr* integrate_wrt(const Expr* integrand, const Expr* var) {
    Expr* call = mk_fn2("Integrate", expr_copy((Expr*)integrand),
                        expr_copy((Expr*)var));
    Expr* r = eval_take(call);
    if (!r) return NULL;
    if (contains_unintegrated(r)) { expr_free(r); return NULL; }
    return r;
}

/* Cancel[Together[Expand[e]]], consuming `e`.  Canonicalises a rational
 * expression so that a difference which is mathematically zero collapses
 * to the literal 0.  All three passes are needed:
 *   - Expand distributes products with sum factors, e.g. g(1+x^2) and
 *     (g + g x^2), so syntactically-different-but-equal terms combine;
 *   - Together puts a difference of fractions over a common denominator,
 *     e.g. 1/(1+(g'/f)^2) vs 1/(f^2+g'^2) (Expand alone cannot);
 *   - Cancel removes the resulting common factor. */
static Expr* canon(Expr* e) {
    if (!e) return NULL;
    Expr* x = eval_take(internal_expand((Expr*[]){ e }, 1));
    if (!x) return NULL;
    Expr* t = eval_take(internal_together((Expr*[]){ x }, 1));
    if (!t) return NULL;
    return eval_take(internal_cancel((Expr*[]){ t }, 1));
}

/* Cancel[Together[a - b]]; borrows a, b; returns owned. */
static Expr* canon_diff(const Expr* a, const Expr* b) {
    Expr* diff = mk_fn2("Plus", expr_copy((Expr*)a),
                        mk_fn2("Times", mk_int(-1), expr_copy((Expr*)b)));
    return canon(diff);
}

/* D[expr, x]; borrows, returns owned. */
static Expr* deriv_dx(const Expr* expr, const Expr* x) {
    Expr* call = expr_new_function(expr_new_symbol(SYM_D),
        (Expr*[]){ expr_copy((Expr*)expr), expr_copy((Expr*)x) }, 2);
    return eval_take(call);
}

/* ---------------------------------------------------------------------- */
/* Linearity-aware residual integral                                      */
/* ---------------------------------------------------------------------- */

/* Build a tidy representation of Integrate[integrand, x] by applying
 * linearity: split over Plus, pull x-free constant factors out of each
 * term, and integrate a bare constant as c*x.  Remaining genuinely
 * un-integrable pieces re-enter the global Integrate dispatcher (which
 * recurses back here for nested unknown-function structure).  Consumes
 * `integrand`; returns an owned, evaluated expression.
 *
 * This mirrors the linearity Mathematica's Integrate applies even to
 * unevaluated results, keeping outputs like  f[x] + Integrate[g[x], x]
 * instead of  f[x] - Integrate[-g[x], x]. */
static Expr* residual_integral(Expr* integrand, const Expr* x) {
    Expr* e = eval_take(internal_expand((Expr*[]){ integrand }, 1));
    if (!e) return NULL;

    /* IMPORTANT: do NOT split a Plus into term-by-term integrals here.
     * A sum such as f'[x]g'[x] + f[x]g''[x] is the exact derivative
     * (f[x]g'[x])' and integrates cleanly as a whole, while its individual
     * terms (e.g. f'[x]g'[x]) have no closed form and send integration-by-
     * parts into an infinite cycle.  The whole sum is handed back to the
     * dispatcher (default case), which re-enters this engine and reduces it
     * as a unit.  Only single-term constant factoring is applied for tidy
     * leftovers. */

    /* A constant (free of x) integrates to c*x. */
    if (!head_is(e, SYM_Plus) && expr_free_of(e, x)) {
        return eval_take(mk_fn2("Times", e, expr_copy((Expr*)x)));
    }

    /* Pull x-free factors out of a product. */
    if (head_is(e, SYM_Times)) {
        size_t n = e->data.function.arg_count;
        Expr** cf = malloc(n * sizeof(Expr*)); size_t cn = 0;
        Expr** rf = malloc(n * sizeof(Expr*)); size_t rn = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* a = e->data.function.args[i];
            if (expr_free_of(a, x)) cf[cn++] = expr_copy(a);
            else                    rf[rn++] = expr_copy(a);
        }
        if (cn > 0 && rn > 0) {
            Expr* c = (cn == 1) ? cf[0]
                    : expr_new_function(expr_new_symbol(SYM_Times), cf, cn);
            Expr* r = (rn == 1) ? rf[0]
                    : expr_new_function(expr_new_symbol(SYM_Times), rf, rn);
            free(cf); free(rf); expr_free(e);
            Expr* integ = mk_fn2("Integrate", r, expr_copy((Expr*)x));
            return eval_take(mk_fn2("Times", c, integ));
        }
        for (size_t i = 0; i < cn; i++) expr_free(cf[i]);
        for (size_t i = 0; i < rn; i++) expr_free(rf[i]);
        free(cf); free(rf);
    }

    /* Default: hand back to the global integrator. */
    return eval_take(mk_fn2("Integrate", e, expr_copy((Expr*)x)));
}

/* Complete a by-parts step  Integrate[Iorig, x] = bdy - Integrate[newI, x].
 * Consumes `bdy` and `newI`; borrows `Iorig`.
 *
 * Self-referential resolution: if the residual integrand is a constant
 * multiple of the original, newI = k*Iorig with k free of x, then
 *   Integrate[Iorig,x] = bdy - k Integrate[Iorig,x]
 *   => Integrate[Iorig,x] = bdy / (1 + k)        (when 1 + k != 0).
 * This closes perfect-power cases such as Integrate[Log[f] f'/f, x] =
 * Log[f]^2/2 (k = 1) and also the ordinary newI == 0 case (k = 0).
 * Otherwise fall back to bdy - Integrate[newI, x]. */
static Expr* finish_byparts(Expr* bdy, Expr* newI, const Expr* Iorig,
                            const Expr* x) {
    Expr* ratio = canon(mk_fn2("Times", expr_copy(newI),
                  mk_fn2("Power", expr_copy((Expr*)Iorig), mk_int(-1))));
    if (ratio && expr_free_of(ratio, x)) {
        Expr* denom = eval_take(mk_fn2("Plus", mk_int(1), expr_copy(ratio)));
        if (denom && !is_lit_zero(denom)) {
            expr_free(ratio); expr_free(newI);
            return eval_take(mk_fn2("Times", bdy,
                             mk_fn2("Power", denom, mk_int(-1))));
        }
        if (denom) expr_free(denom);
    }
    if (ratio) expr_free(ratio);
    Expr* rec = residual_integral(newI, x);   /* consumes newI */
    return eval_take(mk_fn2("Plus", bdy, mk_fn2("Times", mk_int(-1), rec)));
}

/* ---------------------------------------------------------------------- */
/* By-parts reduction for one chosen top atom                             */
/* ---------------------------------------------------------------------- */

/* Attempt the polynomial-part by-parts step using map atom index `ai`
 * (which must have order >= 1) as the top generator A.  `Isub` is the
 * integrand already substituted into the fresh-symbol ring.  On success
 * returns the closed/recursive form  W - Integrate[newI, x]  (owned);
 * returns NULL when this atom does not yield a valid linear reduction. */
static Expr* try_byparts(const AtomMap* m, size_t ai, const Expr* Isub,
                         const Expr* x, const Expr* f_orig) {
    Expr* tA = m->syms[ai];

    /* Linearity: Isub == C1*tA + C0 with C1, C0 free of tA. */
    Expr* C1 = coefficient_k(Isub, tA, 1);
    Expr* C0 = coefficient_k(Isub, tA, 0);
    if (!C1 || !C0) { if (C1) expr_free(C1); if (C0) expr_free(C0); return NULL; }

    Expr* recon = mk_fn2("Plus", mk_fn2("Times", expr_copy(C1), expr_copy(tA)),
                         expr_copy(C0));
    Expr* lin_check = canon_diff(Isub, recon);
    expr_free(recon);
    bool linear = is_lit_zero(lin_check);
    if (lin_check) expr_free(lin_check);
    /* Also require C1 itself to be free of tA (defensive). */
    if (!linear || !expr_free_of(C1, tA)) { expr_free(C1); expr_free(C0); return NULL; }

    /* Antiderivative atom B (order-1 lower) and its ring symbol tB. */
    int order; const Expr* func; const Expr* arg;
    if (!classify_atom(m->atoms[ai], x, &order, &func, &arg) || order < 1) {
        expr_free(C1); expr_free(C0); return NULL;
    }
    Expr* B = make_atom(func, order - 1, arg);
    long bi = map_find(m, B);
    expr_free(B);
    if (bi < 0) { expr_free(C1); expr_free(C0); return NULL; }
    Expr* tB = m->syms[bi];

    /* Chain factor for a composite argument g = arg:  dB/dx = A * g'.
     * The by-parts antiderivative integrates C/g' (not C) with respect to
     * B, since Integrate[C*A, x] = Integrate[(C/g')*dB/dx, x].  For g = x
     * we have g' = 1 and this is a no-op.  Nested unknown functions in the
     * argument (f[g[x]]) are deferred. */
    if (contains_unknown_atom(arg, x)) { expr_free(C1); expr_free(C0); return NULL; }
    Expr* gp = deriv_dx(arg, x);                       /* g' (real, in x) */
    Expr* gp_sub = atoms_to_syms(gp, m);               /* consumes gp */
    Expr* Cp = canon(mk_fn2("Times", expr_copy(C1),
                            mk_fn2("Power", gp_sub, mk_int(-1))));   /* C / g' */
    if (!Cp) { expr_free(C1); expr_free(C0); return NULL; }

    /* W = antiderivative of C/g' with respect to B.  Uses the full
     * integrator over tB, subsuming the polynomial / fraction / log
     * parts (e.g. 1/tB -> Log, 1/tB^2 -> -1/tB, a/(a^2+tB^2) -> ArcTan). */
    Expr* W_sub = integrate_wrt(Cp, tB);
    expr_free(Cp);
    if (!W_sub) { expr_free(C1); expr_free(C0); return NULL; }

    /* Back-substitute to real atoms. */
    Expr* W      = syms_to_atoms(expr_copy(W_sub), m);
    Expr* C1real = syms_to_atoms(expr_copy(C1), m);
    Expr* C0real = syms_to_atoms(expr_copy(C0), m);
    Expr* A      = expr_copy(m->atoms[ai]);
    expr_free(W_sub); expr_free(C1); expr_free(C0);

    /* newI = Cancel[Together[ dW/dx - C1*A - C0 ]]  (must be free of A).
     * Together/Cancel are essential: dW/dx from an ArcTan/Log antiderivative
     * carries a denominator in a different surface form than C1, so plain
     * Expand would not cancel the A-terms and the free-of-A guard would
     * wrongly fail. */
    Expr* dWdx = deriv_dx(W, x);
    Expr* CA   = mk_fn2("Times", C1real, A);          /* consumes C1real, A */
    Expr* resid = mk_fn2("Plus", dWdx,
                         mk_fn2("Times", mk_int(-1),
                                mk_fn2("Plus", CA, C0real)));
    Expr* newI = canon(resid);                        /* consumes resid */

    if (!expr_free_of(newI, m->atoms[ai])) {  /* safety: A must be gone */
        expr_free(W); expr_free(newI);
        return NULL;
    }

    /* Result: W - Integrate[newI, x], resolving the self-referential case
     * and recursing on the residual for nested structure. */
    return finish_byparts(W, newI, f_orig, x);
}

/* ---------------------------------------------------------------------- */
/* By-parts reduction for a Log[eta] generator                            */
/* ---------------------------------------------------------------------- */

/* Reduce an integrand linear in the transcendental generator L = Log[eta]:
 *
 *     Integrate[C*L + D, x] = G*L - Integrate[G*L' - D, x]
 *
 * where G = Integrate[C, x] (the coefficient's antiderivative) and
 * L' = D[eta, x] / eta.  This lowers the integrand's degree in L by one;
 * the residual is free of L and re-enters the global integrator.  Returns
 * NULL if the integrand is not linear in L or if G has no closed form. */
static Expr* try_byparts_log(const Expr* L, const Expr* I, const Expr* x) {
    const Expr* eta = L->data.function.args[0];

    char name[64];
    snprintf(name, sizeof(name), "Integrate`L$%lu", iu_sym_counter++);
    Expr* Lsym = expr_new_symbol(name);

    Expr* from[1] = { (Expr*)L };
    Expr* to[1]   = { Lsym };
    Expr* Isub = subst(expr_copy((Expr*)I), from, to, 1);
    if (!Isub) { expr_free(Lsym); return NULL; }

    Expr* C  = coefficient_k(Isub, Lsym, 1);
    Expr* D0 = coefficient_k(Isub, Lsym, 0);
    if (!C || !D0) { if (C) expr_free(C); if (D0) expr_free(D0);
                     expr_free(Isub); expr_free(Lsym); return NULL; }

    /* Linearity: Isub == C*Lsym + D0, with C, D0 free of Lsym. */
    Expr* recon = mk_fn2("Plus", mk_fn2("Times", expr_copy(C), expr_copy(Lsym)),
                         expr_copy(D0));
    Expr* lin = canon_diff(Isub, recon);
    expr_free(recon); expr_free(Isub);
    bool ok = is_lit_zero(lin) && expr_free_of(C, Lsym) && expr_free_of(D0, Lsym);
    if (lin) expr_free(lin);
    expr_free(Lsym);
    if (!ok) { expr_free(C); expr_free(D0); return NULL; }

    /* G = Integrate[C, x] must close (no leftover Integrate). */
    Expr* G = integrate_wrt(C, x);
    expr_free(C);
    if (!G) { expr_free(D0); return NULL; }

    /* L' = D[eta, x] / eta. */
    Expr* eta_p = deriv_dx(eta, x);
    Expr* Lp = canon(mk_fn2("Times", eta_p,
                            mk_fn2("Power", expr_copy((Expr*)eta), mk_int(-1))));
    if (!Lp) { expr_free(G); expr_free(D0); return NULL; }

    /* newI = G*L' - D0  (free of L). */
    Expr* newI = canon(mk_fn2("Plus", mk_fn2("Times", expr_copy(G), Lp),
                              mk_fn2("Times", mk_int(-1), D0)));
    if (!newI) { expr_free(G); return NULL; }

    /* result = G*L - Integrate[newI, x], resolving the self-referential
     * case (e.g. Integrate[Log[f] f'/f, x] = Log[f]^2/2). */
    Expr* bdy = mk_fn2("Times", G, expr_copy((Expr*)L));
    return finish_byparts(bdy, newI, I, x);
}

/* ---------------------------------------------------------------------- */
/* Core integrator                                                        */
/* ---------------------------------------------------------------------- */

/* True iff `f` contains at least one undefined-function derivative
 * (order >= 1) of `x`.  Cheap gate run before the expensive canonical
 * cycle-detection key is built. */
static bool has_deriv_atom(const Expr* f, const Expr* x) {
    AtomMap m; map_init(&m);
    collect_atoms(f, x, &m);
    bool has = false;
    for (size_t i = 0; i < m.n; i++) if (m.orders[i] >= 1) { has = true; break; }
    map_free(&m);
    return has;
}

static Expr* iu_integrate_core(Expr* f, Expr* x) {
    if (x->type != EXPR_SYMBOL) return NULL;

    /* Transcendental generators (Log[eta] over the unknown-function field)
     * sit above the derivative tower, so reduce them first by parts.  The
     * residual is free of the chosen log and re-enters the integrator. */
    ExprVec logs; evec_init(&logs);
    collect_log_gens(f, x, &logs);
    for (size_t i = 0; i < logs.n; i++) {
        Expr* r = try_byparts_log(logs.items[i], f, x);
        if (r) { evec_free(&logs); return r; }
    }
    bool had_logs = (logs.n > 0);
    evec_free(&logs);

    AtomMap m; map_init(&m);
    collect_atoms(f, x, &m);

    /* Gate: need at least one undefined-function *derivative* (order >= 1).
     * (A log-only integrand with no closed coefficient falls through here.) */
    bool has_deriv = false;
    for (size_t i = 0; i < m.n; i++) if (m.orders[i] >= 1) { has_deriv = true; break; }
    if (!has_deriv) { map_free(&m); (void)had_logs; return NULL; }

    map_complete_chains(&m, x);

    Expr* Isub = atoms_to_syms(expr_copy(f), &m);

    /* Try by-parts on present top atoms, highest order first. */
    Expr* result = NULL;
    int max_order = 0;
    for (size_t i = 0; i < m.n; i++) if (m.orders[i] > max_order) max_order = m.orders[i];

    for (int ord = max_order; ord >= 1 && !result; ord--) {
        for (size_t i = 0; i < m.n && !result; i++) {
            if (m.orders[i] != ord) continue;
            /* Skip atoms not actually present in Isub (chain-only fillers). */
            if (expr_free_of(Isub, m.syms[i])) continue;
            result = try_byparts(&m, i, Isub, x, f);
        }
    }

    expr_free(Isub);
    map_free(&m);
    return result;
}

/* Cycle guard.  Integration by parts can transform an integrand back into
 * one already being processed further up the stack (e.g. the genuinely
 * non-elementary Integrate[f'[x] g'[x], x] -> f g' - Integrate[f g'', x]
 * -> f g' - (f g' - Integrate[f' g', x]) -> ...).  We record the canonical
 * form of every integrand currently in flight and refuse to recurse into a
 * repeat, returning NULL so the dispatcher leaves it unevaluated. */
#define IU_MAX_DEPTH 64
static Expr* iu_stack[IU_MAX_DEPTH];
static int   iu_depth = 0;

static Expr* iu_integrate(Expr* f, Expr* x) {
    if (x->type != EXPR_SYMBOL) return NULL;
    if (!has_deriv_atom(f, x) && !has_log_gen(f, x)) return NULL;

    Expr* key = canon(expr_copy(f));
    if (!key) return NULL;
    for (int i = 0; i < iu_depth; i++) {
        if (expr_eq(iu_stack[i], key)) { expr_free(key); return NULL; }
    }
    if (iu_depth >= IU_MAX_DEPTH) { expr_free(key); return NULL; }
    iu_stack[iu_depth++] = key;

    Expr* r = iu_integrate_core(f, x);

    iu_depth--;
    expr_free(iu_stack[iu_depth]);
    return r;
}

/* ---------------------------------------------------------------------- */
/* Public entry points                                                    */
/* ---------------------------------------------------------------------- */

Expr* integrate_unknown_try(Expr* f, Expr* x) {
    return iu_integrate(f, x);
}

Expr* builtin_integrate_unknown(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    return iu_integrate(f, x);
}

void integrate_unknown_init(void) {
    symtab_add_builtin("Integrate`Undefined", builtin_integrate_unknown);
    symtab_get_def("Integrate`Undefined")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`Undefined",
        "Integrate`Undefined[f, x] integrates expressions that are rational\n"
        "in undefined functions u[x] and their derivatives (Roach 1992, §1.7)\n"
        "by recognising total-derivative and integration-by-parts structure.\n"
        "Returns unevaluated when no closed antiderivative in the generators\n"
        "exists.");
}
