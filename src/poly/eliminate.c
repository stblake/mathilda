/* eliminate.c
 *
 * `Eliminate[eqns, vars]` -- the user-facing variable-elimination front
 * door.  Takes a list / conjunction of `lhs == rhs` equations and a list
 * of variables to eliminate, then drives the lex-order Buchberger engine
 * in groebner.c with the elimination-block layout that GroebnerBasis's
 * 3-arg form already uses internally.  Surviving basis polynomials are
 * re-presented as balanced `Equal[posPart, negPart]` equations and
 * combined with `And` for multiples.
 *
 * The pre-pass tries to handle simple transcendental equations of the
 * shape `f[poly] == const` (or `const == f[poly]`) for invertible
 * elementaries f -- one principal-branch rewrite per layer, with an
 * `Eliminate::ifun` diagnostic to flag that solutions may be missed.
 *
 * Memory contract: the evaluator owns `res`.  Every early-return path
 * frees all temporaries we own (wrapper Lists, all_vars arrays,
 * partially-built GBPoly* arrays); we never `expr_free(res)`.  See
 * SPEC.md §4 for the builtin ownership rule.
 */

#include "eliminate.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arithmetic.h"
#include "attr.h"
#include "eval.h"
#include "groebner.h"
#include "internal.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ */
/*  Small shape predicates                                             */
/* ------------------------------------------------------------------ */

static bool head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == sym;
}

static bool is_list(const Expr* e)  { return head_is(e, SYM_List); }
static bool is_and(const Expr* e)   { return head_is(e, SYM_And); }
static bool is_equal(const Expr* e) { return head_is(e, SYM_Equal); }
static bool is_plus(const Expr* e)  { return head_is(e, SYM_Plus); }
static bool is_times(const Expr* e) { return head_is(e, SYM_Times); }

/* "Polynomial-recursive" heads: structural heads we walk INTO when
 * looking for variable atoms (the args are still candidate atoms).
 * Power is handled separately because we only recurse into its base
 * when the exponent is an integer. */
static bool is_poly_recursive_head(const char* hsym) {
    return hsym == SYM_Plus  || hsym == SYM_Times
        || hsym == SYM_List  || hsym == SYM_And
        || hsym == SYM_Equal || hsym == SYM_Or
        || hsym == SYM_Not;
}

/* Anything that should never be treated as a polynomial variable atom:
 * pure numeric literals and the boolean / infinity sentinels.
 *
 * Mathematical constants like Pi, E, I, EulerGamma are deliberately
 * NOT filtered here — they flow through the polynomial pipeline as
 * parameter symbols (this matches Mathematica's behaviour: `Eliminate[
 * y == Pi/2, y]` returns True via `Pi` as a free parameter). */
static bool is_const_atom(const Expr* e) {
    if (!e) return true;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_REAL:
        case EXPR_BIGINT:
            return true;
        default: break;
    }
    if (e->type == EXPR_MPFR) return true;
    if (e->type == EXPR_SYMBOL) {
        const char* s = e->data.symbol;
        return s == SYM_True || s == SYM_False
            || s == SYM_Infinity || s == SYM_ComplexInfinity
            || s == SYM_Indeterminate;
    }
    if (head_is(e, SYM_Rational)) return true;
    if (head_is(e, SYM_Complex))  return true;
    return false;
}

/* True iff a symbol named `name` (string compare, NOT pointer compare —
 * the callers may pass a non-interned probe) appears anywhere inside e. */
static bool walk_has_symbol_name(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        return strcmp(e->data.symbol, name) == 0;
    }
    if (e->type != EXPR_FUNCTION) return false;
    if (walk_has_symbol_name(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (walk_has_symbol_name(e->data.function.args[i], name)) return true;
    }
    return false;
}

/* ------------------------------------------------------------------ */
/*  Variable-set membership                                            */
/* ------------------------------------------------------------------ */

/* True iff `e` matches any element of `vars[0..n-1]` under structural
 * equality.  Variables may be arbitrary expressions (Mathematica allows
 * `Eliminate[..., {a[1], a[2]}]`), so we use `expr_eq` rather than a
 * pointer or symbol comparison. */
static bool var_in_list(const Expr* e, Expr** vars, size_t n) {
    for (size_t i = 0; i < n; i++) {
        if (expr_eq(e, vars[i])) return true;
    }
    return false;
}

/* True iff `e` contains any sub-expression structurally equal to any of
 * the elim variables.  Conservative: matches at every level.  Used to
 * decide whether `f[u] == v` is safe to invert (we require `v` to be
 * elim-free). */
static bool contains_any_var(const Expr* e, Expr** vars, size_t n) {
    if (!e) return false;
    if (var_in_list(e, vars, n)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_any_var(e->data.function.head, vars, n)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_any_var(e->data.function.args[i], vars, n)) return true;
    }
    return false;
}

/* ------------------------------------------------------------------ */
/*  Symbol collection (variable discovery)                             */
/* ------------------------------------------------------------------ */

typedef struct {
    Expr** items;     /* borrowed pointers, deduped by expr_eq */
    size_t n;
    size_t cap;
} VarSet;

static void varset_init(VarSet* v) {
    v->items = NULL; v->n = 0; v->cap = 0;
}

static void varset_free(VarSet* v) {
    free(v->items);
    v->items = NULL; v->n = 0; v->cap = 0;
}

static void varset_push(VarSet* v, Expr* e) {
    for (size_t i = 0; i < v->n; i++) {
        if (expr_eq(v->items[i], e)) return;
    }
    if (v->n == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 8;
        v->items = (Expr**)realloc(v->items, sizeof(Expr*) * v->cap);
    }
    v->items[v->n++] = e;
}

/* Walk `e` collecting every variable atom that is not in `elim`.
 *
 * An "atom" here is anything we want `gb_from_expr` to treat as a single
 * polynomial indeterminate:
 *   - any non-constant symbol (`x`, `y`, ...);
 *   - any function-shaped expression whose head is NOT one of the
 *     structural operators we recurse through (Plus/Times/Power/List/
 *     And/Equal/...).  This is the rule that turns `Dt[y]`, `f[a,b]`,
 *     `Sin[x]` (when not algebraised) into single variables instead of
 *     mis-decomposing them into `Dt * y`.
 *
 * Heads handled specially:
 *   - Plus, Times, List, And, Equal, Or, Not: recurse into args.
 *   - Power[base, exp]:
 *       * exp is an integer  -> recurse into base only (`x^3` is
 *         polynomial in `x`).
 *       * otherwise          -> treat the whole `Power[]` as an atom
 *         (post-algebraisation this should not occur; safety net).
 *   - any other head         -> treat the whole `Function[...]` as an
 *     atom (do NOT recurse into args — we don't want to add `y` from
 *     `Dt[y]` to the main-vars set).
 *
 * Elim-list membership and constant-atom checks are by `expr_eq` so an
 * elim entry like `Dt[x]` correctly suppresses the entire function. */
static void collect_main_vars(const Expr* e, Expr** elim, size_t n_elim,
                              VarSet* out) {
    if (!e || is_const_atom(e)) return;
    if (var_in_list(e, elim, n_elim)) return;
    if (e->type == EXPR_SYMBOL) {
        varset_push(out, (Expr*)e);
        return;
    }
    if (e->type != EXPR_FUNCTION) return;

    const char* hsym = (e->data.function.head
                        && e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol : NULL;

    if (hsym && is_poly_recursive_head(hsym)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            collect_main_vars(e->data.function.args[i], elim, n_elim, out);
        }
        return;
    }
    if (hsym == SYM_Power && e->data.function.arg_count == 2
        && e->data.function.args[1]->type == EXPR_INTEGER) {
        collect_main_vars(e->data.function.args[0], elim, n_elim, out);
        return;
    }
    /* Function-shaped atom: register the whole expression. */
    varset_push(out, (Expr*)e);
}

/* ------------------------------------------------------------------ */
/*  Transcendental pre-pass                                            */
/* ------------------------------------------------------------------ */

/* For invertible-elementary heads, return the principal-branch inverse
 * symbol (interned).  Returns NULL if `head_sym` is not in the table. */
static const char* invertible_head(const char* head_sym) {
    if (head_sym == SYM_Sin)  return SYM_ArcSin;
    if (head_sym == SYM_Cos)  return SYM_ArcCos;
    if (head_sym == SYM_Tan)  return SYM_ArcTan;
    if (head_sym == SYM_Sinh) return SYM_ArcSinh;
    if (head_sym == SYM_Cosh) return SYM_ArcCosh;
    if (head_sym == SYM_Tanh) return SYM_ArcTanh;
    if (head_sym == SYM_Exp)  return SYM_Log;
    if (head_sym == SYM_Log)  return SYM_Exp;
    if (head_sym == SYM_ArcSin)  return SYM_Sin;
    if (head_sym == SYM_ArcCos)  return SYM_Cos;
    if (head_sym == SYM_ArcTan)  return SYM_Tan;
    if (head_sym == SYM_ArcSinh) return SYM_Sinh;
    if (head_sym == SYM_ArcCosh) return SYM_Cosh;
    if (head_sym == SYM_ArcTanh) return SYM_Tanh;
    return NULL;
}

/* If `side` is `f[u]` with `f` invertible, return its inverse symbol
 * string and store `u` into `*inner_out` (borrowed, not copied).  Else
 * return NULL. */
static const char* peel_invertible(Expr* side, Expr** inner_out) {
    if (!side || side->type != EXPR_FUNCTION) return NULL;
    if (!side->data.function.head
     || side->data.function.head->type != EXPR_SYMBOL) return NULL;
    if (side->data.function.arg_count != 1) return NULL;
    const char* inv = invertible_head(side->data.function.head->data.symbol);
    if (!inv) return NULL;
    *inner_out = side->data.function.args[0];
    return inv;
}

/* Attempt to rewrite a single `Equal[lhs, rhs]` so the elim variables
 * appear only in polynomial position.  Returns a fresh `Equal[]` Expr
 * the caller owns.  On no rewrite, returns a fresh copy of `eq`.
 * Sets `*fired = true` if any inverse-function rewrite applied. */
static Expr* try_inverse_rewrite(const Expr* eq, Expr** elim, size_t n_elim,
                                 bool* fired) {
    if (!is_equal(eq) || eq->data.function.arg_count != 2) {
        return expr_copy((Expr*)eq);
    }
    Expr* lhs = eq->data.function.args[0];
    Expr* rhs = eq->data.function.args[1];

    /* Iterate while we can peel a layer.  Cap at a few rounds so a
     * pathological alternation can't loop forever. */
    Expr* cur_lhs = expr_copy(lhs);
    Expr* cur_rhs = expr_copy(rhs);
    for (int round = 0; round < 6; round++) {
        bool peeled = false;
        /* Case 1: f[u] == v with v elim-free and u elim-relevant. */
        Expr* inner = NULL;
        const char* inv = peel_invertible(cur_lhs, &inner);
        if (inv && !contains_any_var(cur_rhs, elim, n_elim)
                && contains_any_var(inner, elim, n_elim)) {
            Expr* new_rhs_arg = expr_copy(cur_rhs);
            Expr* new_rhs = expr_new_function(expr_new_symbol(inv),
                (Expr*[]){ new_rhs_arg }, 1);
            Expr* new_lhs = expr_copy(inner);
            expr_free(cur_lhs); expr_free(cur_rhs);
            cur_lhs = new_lhs; cur_rhs = new_rhs;
            peeled = true; *fired = true;
        } else {
            /* Case 2: mirror -- v == f[u]. */
            inv = peel_invertible(cur_rhs, &inner);
            if (inv && !contains_any_var(cur_lhs, elim, n_elim)
                    && contains_any_var(inner, elim, n_elim)) {
                Expr* new_lhs_arg = expr_copy(cur_lhs);
                Expr* new_lhs = expr_new_function(expr_new_symbol(inv),
                    (Expr*[]){ new_lhs_arg }, 1);
                Expr* new_rhs = expr_copy(inner);
                expr_free(cur_lhs); expr_free(cur_rhs);
                cur_lhs = new_lhs; cur_rhs = new_rhs;
                peeled = true; *fired = true;
            }
        }
        if (!peeled) break;
    }
    /* Evaluate the inverse-function applications so simple cases like
     * `ArcSin[1]` fold to closed form.  evaluate() returns a fresh tree,
     * so we still own the inputs and must free them. */
    Expr* ev_lhs = evaluate(cur_lhs);
    Expr* ev_rhs = evaluate(cur_rhs);
    expr_free(cur_lhs);
    expr_free(cur_rhs);
    Expr* out = expr_new_function(expr_new_symbol("Equal"),
        (Expr*[]){ ev_lhs, ev_rhs }, 2);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Algebraisation pre-pass                                            */
/* ------------------------------------------------------------------ */

/* Replace every `Power[base, p/q]` whose `base` mentions an elim
 * variable with `Power[aux, p*L/q]` where `aux` is a fresh symbol and
 * `L = lcm({q's seen for this base})`.  In tandem we add an algebraic
 * constraint `Power[aux, L] == base` to the system and push `aux` onto
 * the elim list.
 *
 * Effect: the system becomes polynomial in the (extended) variable
 * set, and any output mentioning `aux` is dropped by the elimination
 * filter just like an original elim variable would be.  The
 * generically-true cross-multiplied form is what Mathematica's
 * `Eliminate` returns for radical inputs — see the headline example
 * in `tests/test_eliminate.c::test_radical_dt_y`. */

typedef struct {
    Expr**   bases;   /* owned deep copies of the unique algebraic bases */
    int64_t* lcms;    /* per-base LCM of all rational-exponent denoms    */
    Expr**   auxs;    /* owned EXPR_SYMBOL nodes, one per base           */
    size_t   n;
    size_t   cap;
} AlgState;

static void alg_init(AlgState* a) {
    a->bases = NULL; a->lcms = NULL; a->auxs = NULL;
    a->n = 0; a->cap = 0;
}

static void alg_free(AlgState* a) {
    if (a->bases) {
        for (size_t i = 0; i < a->n; i++) expr_free(a->bases[i]);
    }
    if (a->auxs) {
        for (size_t i = 0; i < a->n; i++) expr_free(a->auxs[i]);
    }
    free(a->bases);
    free(a->lcms);
    free(a->auxs);
    alg_init(a);
}

/* Register a (base, denom) pair, deduping by `expr_eq` on the base.
 * `q` is the denominator of the rational exponent and is normalised to
 * positive here so `Sqrt[u]` and `1/Sqrt[u]` collapse onto the same
 * base. */
static void alg_add(AlgState* a, const Expr* base, int64_t q) {
    if (q < 0) q = -q;
    if (q <= 0) return;
    for (size_t i = 0; i < a->n; i++) {
        if (expr_eq(a->bases[i], base)) {
            a->lcms[i] = lcm(a->lcms[i], q);
            return;
        }
    }
    if (a->n == a->cap) {
        a->cap = a->cap ? a->cap * 2 : 4;
        a->bases = (Expr**)realloc(a->bases, sizeof(Expr*) * a->cap);
        a->lcms  = (int64_t*)realloc(a->lcms,  sizeof(int64_t) * a->cap);
    }
    a->bases[a->n] = expr_copy((Expr*)base);
    a->lcms[a->n]  = q;
    a->n++;
}

/* Walk `e`, registering every `Power[base, Rational[p, q]]` with q != 1
 * provided `base` mentions at least one elim variable.  Bases without
 * any elim variable are left alone: they're parameter constants like
 * `Sqrt[2]` that we cannot eliminate without losing information, and
 * the user can pre-rationalise such inputs manually if needed. */
static void alg_collect(const Expr* e, AlgState* a,
                        Expr** elim, size_t n_elim) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        int64_t p, q;
        if (is_rational(exp, &p, &q) && q != 1
            && contains_any_var(base, elim, n_elim)) {
            alg_add(a, base, q);
        }
    }
    if (e->data.function.head) {
        alg_collect(e->data.function.head, a, elim, n_elim);
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        alg_collect(e->data.function.args[i], a, elim, n_elim);
    }
}

/* Build a fresh tree replacing each registered `Power[base, p/q]` with
 * `Power[aux, p*L/q]` (collapsed to `1`, `aux`, or `Power[aux, k]` as
 * the exponent simplifies).  Returns a freshly-allocated Expr* the
 * caller owns. */
static Expr* alg_substitute(const Expr* e, const AlgState* a) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        int64_t p, q;
        if (is_rational(exp, &p, &q) && q != 1) {
            for (size_t i = 0; i < a->n; i++) {
                if (!expr_eq(a->bases[i], base)) continue;
                int64_t L = a->lcms[i];
                /* L is the LCM of all q's seen for this base, so q | L
                 * and (p * L / q) is exactly an integer. */
                int64_t new_exp = (p * L) / q;
                if (new_exp == 0) return expr_new_integer(1);
                if (new_exp == 1) return expr_copy(a->auxs[i]);
                Expr** pa = (Expr**)malloc(sizeof(Expr*) * 2);
                pa[0] = expr_copy(a->auxs[i]);
                pa[1] = expr_new_integer(new_exp);
                Expr* r = expr_new_function(expr_new_symbol("Power"), pa, 2);
                free(pa);
                return r;
            }
        }
    }

    Expr* new_head = alg_substitute(e->data.function.head, a);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = alg_substitute(e->data.function.args[i], a);
    }
    Expr* r = expr_new_function(new_head, new_args, n);
    free(new_args);
    return r;
}

/* `Power[aux_i, L_i] == base_i` (the unsigned algebraic constraint).
 *
 * The base is first run through `alg_substitute` so any nested
 * algebraic atoms inside it (e.g. `Sqrt[x]` inside the base of
 * `Sqrt[x + Sqrt[x]]`) are rewritten in terms of their own aux vars.
 * Without this, nested-radical inputs leave a `Power[u, 1/2]` inside
 * the constraint and the polynomial pipeline correctly rejects them. */
static Expr* alg_make_constraint(const AlgState* a, size_t i) {
    Expr** pa = (Expr**)malloc(sizeof(Expr*) * 2);
    pa[0] = expr_copy(a->auxs[i]);
    pa[1] = expr_new_integer(a->lcms[i]);
    Expr* pow_e = expr_new_function(expr_new_symbol("Power"), pa, 2);
    free(pa);

    Expr** ea = (Expr**)malloc(sizeof(Expr*) * 2);
    ea[0] = pow_e;
    ea[1] = alg_substitute(a->bases[i], a);
    Expr* eq = expr_new_function(expr_new_symbol("Equal"), ea, 2);
    free(ea);
    return eq;
}

/* Pick a symbol name not appearing anywhere in `eqs`, `elim`, the
 * registered bases, or the auxes assigned so far.  The chosen name is
 * interned and returned (the interner is global; the same string maps
 * to the same `const char*` across the program). */
static const char* alg_fresh_name(Expr** eqs, size_t n_eq,
                                  Expr** elim, size_t n_elim,
                                  const AlgState* a, size_t i_so_far,
                                  int* next_idx) {
    char buf[64];
    for (int tries = 0; tries < 10000; tries++) {
        snprintf(buf, sizeof(buf), "$el%d$", *next_idx);
        (*next_idx)++;
        bool collide = false;
        for (size_t j = 0; j < n_eq && !collide; j++) {
            if (walk_has_symbol_name(eqs[j], buf)) collide = true;
        }
        for (size_t j = 0; j < n_elim && !collide; j++) {
            if (walk_has_symbol_name(elim[j], buf)) collide = true;
        }
        for (size_t j = 0; j < a->n && !collide; j++) {
            if (walk_has_symbol_name(a->bases[j], buf)) collide = true;
        }
        for (size_t j = 0; j < i_so_far && !collide; j++) {
            if (walk_has_symbol_name(a->auxs[j], buf)) collide = true;
        }
        if (!collide) return intern_symbol(buf);
    }
    /* Vanishingly unlikely fallback. */
    return intern_symbol("$el$overflow$");
}

/* ------------------------------------------------------------------ */
/*  Equation -> polynomial                                             */
/* ------------------------------------------------------------------ */

/* Reduce `Equal[a, b]` to the polynomial `a - b`.
 *
 * Pipeline:  Together  ->  Numerator  ->  Expand  ->  evaluate.
 *
 * The Together+Numerator step clears any negative powers
 * (`Power[t, -1]`) introduced by the algebraisation pre-pass — the
 * downstream `gb_from_expr` only accepts non-negative integer
 * exponents, so we have to pull every `t` out of the denominator
 * BEFORE building the polynomial.  For a bare polynomial input this is
 * idempotent: Together returns the polynomial unchanged and Numerator
 * is the identity.  We still Expand at the end to canonicalise to
 * `Plus[Times[...], ...]` form, which `expr_term` expects. */
static Expr* normalise_equation(const Expr* eq) {
    Expr* base;
    if (is_equal(eq) && eq->data.function.arg_count == 2) {
        Expr* a = expr_copy(eq->data.function.args[0]);
        Expr* b = expr_copy(eq->data.function.args[1]);
        base = internal_subtract((Expr*[]){ a, b }, 2);
    } else {
        base = expr_copy((Expr*)eq);
    }
    Expr* together = internal_together((Expr*[]){ base }, 1);
    Expr* numer    = internal_numerator((Expr*[]){ together }, 1);
    Expr* expanded = internal_expand((Expr*[]){ numer }, 1);
    Expr* normalised = evaluate(expanded);
    expr_free(expanded);
    return normalised;
}

/* ------------------------------------------------------------------ */
/*  Term-by-sign classification + balanced Equal building              */
/* ------------------------------------------------------------------ */

/* +1 / -1 / 0 ("don't know -- treat as positive").  Reads the leading
 * coefficient of a single term emitted by `gb_to_expr`. */
static int term_sign(const Expr* t) {
    if (!t) return 1;
    if (t->type == EXPR_INTEGER) {
        return (t->data.integer < 0) ? -1 : 1;
    }
    if (t->type == EXPR_BIGINT) {
        return (mpz_sgn(t->data.bigint) < 0) ? -1 : 1;
    }
    if (head_is(t, SYM_Rational) && t->data.function.arg_count == 2) {
        return term_sign(t->data.function.args[0]);
    }
    if (is_times(t) && t->data.function.arg_count >= 1) {
        return term_sign(t->data.function.args[0]);
    }
    return 1;
}

/* Build a fresh copy of `t` with its leading coefficient negated. */
static Expr* term_negate(const Expr* t) {
    if (!t) return expr_new_integer(0);
    if (t->type == EXPR_INTEGER) {
        return expr_new_integer(-t->data.integer);
    }
    if (t->type == EXPR_BIGINT) {
        mpz_t neg; mpz_init(neg);
        mpz_neg(neg, t->data.bigint);
        Expr* r = expr_new_bigint_from_mpz(neg);
        mpz_clear(neg);
        return r;
    }
    if (head_is(t, SYM_Rational) && t->data.function.arg_count == 2) {
        Expr* num_neg = term_negate(t->data.function.args[0]);
        Expr* den = expr_copy(t->data.function.args[1]);
        return expr_new_function(expr_new_symbol("Rational"),
            (Expr*[]){ num_neg, den }, 2);
    }
    if (is_times(t) && t->data.function.arg_count >= 1) {
        size_t n = t->data.function.arg_count;
        Expr* lead = t->data.function.args[0];
        if (lead->type == EXPR_INTEGER && lead->data.integer == -1) {
            /* Strip the leading -1.  If only one other factor remains,
             * return it directly; otherwise rebuild Times[rest...]. */
            if (n == 2) return expr_copy(t->data.function.args[1]);
            Expr** args = (Expr**)malloc(sizeof(Expr*) * (n - 1));
            for (size_t i = 0; i < n - 1; i++) {
                args[i] = expr_copy(t->data.function.args[i + 1]);
            }
            Expr* r = expr_new_function(expr_new_symbol("Times"), args, n - 1);
            free(args);
            return r;
        }
        if (lead->type == EXPR_INTEGER
         || lead->type == EXPR_BIGINT
         || head_is(lead, SYM_Rational)) {
            /* Negate the leading coefficient in-place in a fresh copy. */
            Expr** args = (Expr**)malloc(sizeof(Expr*) * n);
            args[0] = term_negate(lead);
            for (size_t i = 1; i < n; i++) {
                args[i] = expr_copy(t->data.function.args[i]);
            }
            Expr* r = expr_new_function(expr_new_symbol("Times"), args, n);
            free(args);
            return r;
        }
        /* No coefficient -- prepend -1. */
        Expr** args = (Expr**)malloc(sizeof(Expr*) * (n + 1));
        args[0] = expr_new_integer(-1);
        for (size_t i = 0; i < n; i++) {
            args[i + 1] = expr_copy(t->data.function.args[i]);
        }
        Expr* r = expr_new_function(expr_new_symbol("Times"), args, n + 1);
        free(args);
        return r;
    }
    /* Symbol, Power, or anything else: wrap as Times[-1, t]. */
    Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
    args[0] = expr_new_integer(-1);
    args[1] = expr_copy((Expr*)t);
    Expr* r = expr_new_function(expr_new_symbol("Times"), args, 2);
    free(args);
    return r;
}

/* Build `Equal[posPart, -negPart]` from `p = lhs - rhs`.  Calls evaluate
 * at the end so Plus[x] collapses to x and constants fold. */
static Expr* balance_polynomial(const GBPoly* g, Expr** all_vars) {
    Expr* poly = gb_to_expr(g, all_vars);

    /* Normalise: a single term becomes a 1-arg Plus for uniform splitting. */
    Expr** terms = NULL;
    size_t n_terms = 0;
    bool plus_wrapped = is_plus(poly);
    if (plus_wrapped) {
        n_terms = poly->data.function.arg_count;
        terms = (Expr**)malloc(sizeof(Expr*) * (n_terms > 0 ? n_terms : 1));
        for (size_t i = 0; i < n_terms; i++) {
            terms[i] = poly->data.function.args[i];
        }
    } else {
        n_terms = 1;
        terms = (Expr**)malloc(sizeof(Expr*) * 1);
        terms[0] = poly;
    }

    /* Split by sign. */
    Expr** pos = (Expr**)malloc(sizeof(Expr*) * n_terms);
    Expr** neg = (Expr**)malloc(sizeof(Expr*) * n_terms);
    size_t n_pos = 0, n_neg = 0;
    for (size_t i = 0; i < n_terms; i++) {
        int s = term_sign(terms[i]);
        if (s < 0) {
            neg[n_neg++] = term_negate(terms[i]);
        } else {
            pos[n_pos++] = expr_copy(terms[i]);
        }
    }
    free(terms);
    expr_free(poly);

    /* Pack each side into Integer 0 / single Expr / Plus[...]. */
    Expr* lhs;
    if (n_pos == 0) {
        lhs = expr_new_integer(0);
    } else if (n_pos == 1) {
        lhs = pos[0];
    } else {
        lhs = expr_new_function(expr_new_symbol("Plus"), pos, n_pos);
    }
    free(pos);

    Expr* rhs;
    if (n_neg == 0) {
        rhs = expr_new_integer(0);
    } else if (n_neg == 1) {
        rhs = neg[0];
    } else {
        rhs = expr_new_function(expr_new_symbol("Plus"), neg, n_neg);
    }
    free(neg);

    Expr* eq = expr_new_function(expr_new_symbol("Equal"),
        (Expr*[]){ lhs, rhs }, 2);
    Expr* evd = evaluate(eq);
    expr_free(eq);
    return evd;
}

/* Divide `p` by the largest monomial X^(m_1, ..., m_n) that divides
 * every term — equivalently, for each variable v subtract
 * `m_v = min_t exp[t][v]` from every term's exponent for v.
 *
 * After algebraisation+Buchberger we often pick up an extraneous
 * factor of a main variable from cross-multiplication (e.g. a final
 * `u^5 Dt[y] - ...` that is just `u` times the primitive form
 * `u^4 Dt[y] - ...`).  Stripping the monomial factor makes our output
 * match Mathematica's reduced-form `Eliminate` result and shortens
 * downstream `Factor` / `Solve` work.
 *
 * Generically valid: dividing both sides of `P == 0` by `u^m` is
 * equivalent to the original equation as long as `u != 0`.  Eliminate
 * already returns the generic consequence, so this is consistent with
 * the rest of the pipeline. */
static void gb_poly_strip_monomial(GBPoly* p) {
    if (p->n_terms == 0 || p->n_vars == 0) return;
    int* mins = (int*)malloc(sizeof(int) * (size_t)p->n_vars);
    for (int v = 0; v < p->n_vars; v++) {
        mins[v] = p->exps[v];  /* first term */
    }
    for (size_t t = 1; t < p->n_terms; t++) {
        const int* row = &p->exps[t * (size_t)p->n_vars];
        for (int v = 0; v < p->n_vars; v++) {
            if (row[v] < mins[v]) mins[v] = row[v];
        }
    }
    bool any = false;
    for (int v = 0; v < p->n_vars; v++) {
        if (mins[v] > 0) { any = true; break; }
    }
    if (any) {
        for (size_t t = 0; t < p->n_terms; t++) {
            int* row = &p->exps[t * (size_t)p->n_vars];
            for (int v = 0; v < p->n_vars; v++) row[v] -= mins[v];
        }
    }
    free(mins);
}

/* ------------------------------------------------------------------ */
/*  Diagnostics                                                        */
/* ------------------------------------------------------------------ */

static void emit_argt(size_t argc) {
    fprintf(stderr, "Eliminate::argt: Eliminate called with %zu argument(s); "
                    "2 expected.\n", argc);
}

static void emit_eqf(void) {
    fprintf(stderr, "Eliminate::eqf: equations must be given as Equal[lhs, rhs] "
                    "(==) or a list/And of such.\n");
}

static void emit_nlin(void) {
    fprintf(stderr, "Eliminate::nlin: system is not polynomial in the given "
                    "variables.\n");
}

static void emit_ifun(void) {
    fprintf(stderr, "Eliminate::ifun: Inverse functions are being used by "
                    "Eliminate, so some solutions may not be found; use Reduce "
                    "for complete solution information.\n");
}

/* Emit at most once per distinct input expression.  `RepeatedTiming`,
 * `Table[...]`, etc. evaluate the same `Eliminate[...]` form hundreds
 * of times; without the hash guard each iteration re-prints the warning
 * and floods the terminal.  Mirrors `solve.c:warn_bad_option`. */
static void emit_alg(const Expr* res) {
    static uint64_t last_warned_hash = 0;
    uint64_t h = res ? expr_hash((Expr*)res) : 0;
    if (h == last_warned_hash) return;
    last_warned_hash = h;
    fprintf(stderr, "Eliminate::alg: Radical (rational-power) subexpressions "
                    "were replaced by auxiliary variables; the returned "
                    "polynomial relation is the cross-multiplied generic "
                    "consequence (sign / branch information may be lost).\n");
}

/* ------------------------------------------------------------------ */
/*  Builtin entry                                                      */
/* ------------------------------------------------------------------ */

Expr* builtin_eliminate(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) {
        emit_argt(argc);
        return NULL;
    }

    Expr* eqns_arg = res->data.function.args[0];
    Expr* vars_arg = res->data.function.args[1];

    /* ----- Sentinel arguments produced by evaluator pre-evaluation -----
     * `Eliminate[1==2, x]` reaches us as `Eliminate[False, x]` because
     * Eliminate is not HoldAll; same for `True`.  Fold these into the
     * obvious answer without emitting `eqf`. */
    if (eqns_arg && eqns_arg->type == EXPR_SYMBOL) {
        if (eqns_arg->data.symbol == SYM_True)  return expr_new_symbol("True");
        if (eqns_arg->data.symbol == SYM_False) return expr_new_symbol("False");
    }

    /* ----- Unpack equation list (filtering True/False sentinels) ----- */
    Expr** raw_in = NULL;
    size_t raw_n = 0;
    if (is_list(eqns_arg) || is_and(eqns_arg)) {
        raw_n = eqns_arg->data.function.arg_count;
        raw_in = eqns_arg->data.function.args;  /* borrowed */
    } else {
        raw_n = 1;
        raw_in = &eqns_arg;  /* borrowed single-element slot */
    }
    if (raw_n == 0) {
        /* Empty system: no constraint at all. */
        return expr_new_symbol("True");
    }
    /* Pre-filter for True/False sentinels mixed inside a list/And. */
    Expr** eq_in = (Expr**)malloc(sizeof(Expr*) * raw_n);
    size_t n_eq = 0;
    for (size_t i = 0; i < raw_n; i++) {
        Expr* x = raw_in[i];
        if (x && x->type == EXPR_SYMBOL) {
            if (x->data.symbol == SYM_True) continue;  /* drop */
            if (x->data.symbol == SYM_False) {
                /* Any False element makes the whole conjunction False. */
                free(eq_in);
                return expr_new_symbol("False");
            }
        }
        eq_in[n_eq++] = x;
    }
    if (n_eq == 0) {
        free(eq_in);
        return expr_new_symbol("True");
    }

    /* ----- Unpack elim-variable list (allow bare symbol/expr) ----- */
    Expr* wrap_vars = NULL;
    Expr** elim_items = NULL;
    size_t n_elim = 0;
    if (is_list(vars_arg)) {
        n_elim = vars_arg->data.function.arg_count;
        elim_items = vars_arg->data.function.args;
    } else {
        Expr** wrapped = (Expr**)malloc(sizeof(Expr*));
        wrapped[0] = expr_copy(vars_arg);
        wrap_vars = expr_new_function(expr_new_symbol("List"), wrapped, 1);
        free(wrapped);
        n_elim = 1;
        elim_items = wrap_vars->data.function.args;
    }
    if (n_elim == 0) {
        /* `Eliminate[eqns, {}]` -- nothing to eliminate; just return the
         * original equations as `And` (or single `Equal`). */
        if (wrap_vars) expr_free(wrap_vars);
        if (n_eq == 1) return expr_copy(eq_in[0]);
        Expr** copies = (Expr**)malloc(sizeof(Expr*) * n_eq);
        for (size_t i = 0; i < n_eq; i++) copies[i] = expr_copy(eq_in[i]);
        Expr* and_e = expr_new_function(expr_new_symbol("And"), copies, n_eq);
        free(copies);
        Expr* ev = evaluate(and_e);
        expr_free(and_e);
        return ev;
    }

    /* ----- Transcendental pre-pass ----- */
    bool fired_ifun = false;
    Expr** rewritten = (Expr**)malloc(sizeof(Expr*) * n_eq);
    for (size_t i = 0; i < n_eq; i++) {
        rewritten[i] = try_inverse_rewrite(eq_in[i], elim_items, n_elim,
                                           &fired_ifun);
    }
    /* `eq_in` is just a working array of borrowed pointers; the source
     * (raw_in / arg list) still owns the real expressions.  Free it now. */
    free(eq_in);
    if (fired_ifun) emit_ifun();

    /* Validate every rewritten equation is `Equal[lhs, rhs]`. */
    for (size_t i = 0; i < n_eq; i++) {
        if (!is_equal(rewritten[i])
         || rewritten[i]->data.function.arg_count != 2) {
            emit_eqf();
            for (size_t k = 0; k < n_eq; k++) expr_free(rewritten[k]);
            free(rewritten);
            if (wrap_vars) expr_free(wrap_vars);
            return NULL;
        }
    }

    /* ----- Algebraisation pre-pass -----
     * Identifies `Power[base, p/q]` (q != 1) subexpressions whose base
     * mentions an elim variable, replaces each with `Power[aux, p*L/q]`
     * for fresh `aux`, and appends `aux^L == base` to the system + `aux`
     * to the elim list.  Drops out cheaply (alg.n == 0) when the input
     * is already polynomial in the elim vars. */
    AlgState alg; alg_init(&alg);
    for (size_t i = 0; i < n_eq; i++) {
        alg_collect(rewritten[i], &alg, elim_items, n_elim);
    }
    Expr** alg_elim_array = NULL;
    if (alg.n > 0) {
        /* Assign fresh aux symbols. */
        alg.auxs = (Expr**)calloc(alg.n, sizeof(Expr*));
        int next_idx = 1;
        for (size_t i = 0; i < alg.n; i++) {
            const char* nm = alg_fresh_name(rewritten, n_eq,
                                            elim_items, n_elim,
                                            &alg, i, &next_idx);
            alg.auxs[i] = expr_new_symbol(nm);
        }
        /* Substitute every equation in-place: build new tree, free old. */
        Expr** combined = (Expr**)malloc(sizeof(Expr*) * (n_eq + alg.n));
        for (size_t i = 0; i < n_eq; i++) {
            combined[i] = alg_substitute(rewritten[i], &alg);
            expr_free(rewritten[i]);
        }
        for (size_t i = 0; i < alg.n; i++) {
            combined[n_eq + i] = alg_make_constraint(&alg, i);
        }
        free(rewritten);
        rewritten = combined;
        n_eq += alg.n;
        /* Extend elim list with aux symbols.  Borrows pointers; we own
         * the array itself and must free it at exit. */
        alg_elim_array = (Expr**)malloc(sizeof(Expr*) * (n_elim + alg.n));
        for (size_t i = 0; i < n_elim; i++) {
            alg_elim_array[i] = elim_items[i];
        }
        for (size_t i = 0; i < alg.n; i++) {
            alg_elim_array[n_elim + i] = alg.auxs[i];
        }
        elim_items = alg_elim_array;
        n_elim    += alg.n;
        emit_alg(res);
    }

    /* ----- Normalise each equation to polynomial form ----- */
    Expr** polys = (Expr**)malloc(sizeof(Expr*) * n_eq);
    for (size_t i = 0; i < n_eq; i++) {
        polys[i] = normalise_equation(rewritten[i]);
    }
    for (size_t i = 0; i < n_eq; i++) expr_free(rewritten[i]);
    free(rewritten);

    /* ----- Variable discovery: collect main variables ----- */
    VarSet mains; varset_init(&mains);
    for (size_t i = 0; i < n_eq; i++) {
        collect_main_vars(polys[i], elim_items, n_elim, &mains);
    }
    /* main_vars borrows pointers into `polys` -- we must not free `polys`
     * before we're done with `all_vars` below. */
    /* If any function-shaped main atom (e.g. `Sin[x*y]`) STILL contains
     * an elim variable as a subexpression, treating it as a single
     * polynomial variable would silently lose the relationship between
     * that atom and the elim var.  Bail out with nlin instead so the
     * caller sees the input returned unevaluated rather than a wrong
     * answer. */
    {
        bool nlin_in_atom = false;
        for (size_t i = 0; i < mains.n; i++) {
            Expr* m = mains.items[i];
            if (m && m->type == EXPR_FUNCTION
                && contains_any_var(m, elim_items, n_elim)) {
                nlin_in_atom = true;
                break;
            }
        }
        if (nlin_in_atom) {
            emit_nlin();
            for (size_t i = 0; i < n_eq; i++) expr_free(polys[i]);
            free(polys);
            varset_free(&mains);
            if (wrap_vars) expr_free(wrap_vars);
            free(alg_elim_array);
            alg_free(&alg);
            return NULL;
        }
    }
    size_t n_main = mains.n;
    size_t n_vars = n_elim + n_main;
    if (n_vars == 0) {
        /* Nothing variable in the system at all: every equation must be
         * a numeric identity.  Fold into True/False by evaluating each. */
        bool any_false = false;
        for (size_t i = 0; i < n_eq; i++) {
            Expr* zero_eq = polys[i];
            if (zero_eq && zero_eq->type == EXPR_INTEGER
                && zero_eq->data.integer != 0) {
                any_false = true; break;
            }
            if (zero_eq && zero_eq->type == EXPR_BIGINT
                && mpz_sgn(zero_eq->data.bigint) != 0) {
                any_false = true; break;
            }
        }
        for (size_t i = 0; i < n_eq; i++) expr_free(polys[i]);
        free(polys);
        varset_free(&mains);
        if (wrap_vars) expr_free(wrap_vars);
        free(alg_elim_array);
        alg_free(&alg);
        return expr_new_symbol(any_false ? "False" : "True");
    }

    /* ----- Build joint variable array (elim block first) -----
     * Within the main block, place main variables in *reverse* order of
     * first appearance.  Under lex this makes later-mentioned variables
     * smaller, so they tend to end up as the "free parameter" of the
     * triangular basis -- matching the convention Mathematica's
     * Eliminate uses (e.g. `Eliminate[{f==x^5+y^5, a==x+y, b==x*y},
     * {x,y}]` yields `f == a^5 - 5 a^3 b + 5 a b^2` with `f` free). */
    Expr** all_vars = (Expr**)malloc(sizeof(Expr*) * n_vars);
    for (size_t i = 0; i < n_elim; i++) all_vars[i] = elim_items[i];
    for (size_t i = 0; i < n_main; i++) {
        all_vars[n_elim + i] = mains.items[n_main - 1 - i];
    }

    /* ----- Convert each polynomial to GBPoly ----- */
    GBPoly** F = (GBPoly**)malloc(sizeof(GBPoly*) * n_eq);
    size_t nF = 0;
    bool nlin = false;
    bool inconsistent_const = false;
    for (size_t i = 0; i < n_eq; i++) {
        GBPoly* p = gb_from_expr(polys[i], all_vars, (int)n_vars,
                                 GB_ORDER_LEX, (int)n_elim, NULL);
        if (!p) { nlin = true; break; }
        if (p->n_terms == 0) { gb_poly_free(p); continue; }
        if (gb_poly_is_constant(p)) {
            /* Equation `c == 0` with `c != 0`: contradiction. */
            inconsistent_const = true;
            gb_poly_free(p);
            break;
        }
        F[nF++] = p;
    }

    /* Cleanup that's safe to do now: input polynomials and main-vars
     * borrow array.  `all_vars[n_elim..]` borrows into `polys[i]`, so
     * keep those alive until we're done with gb_to_expr too. */

    if (nlin) {
        emit_nlin();
        for (size_t i = 0; i < nF; i++) gb_poly_free(F[i]);
        free(F);
        free(all_vars);
        for (size_t i = 0; i < n_eq; i++) expr_free(polys[i]);
        free(polys);
        varset_free(&mains);
        if (wrap_vars) expr_free(wrap_vars);
        free(alg_elim_array);
        alg_free(&alg);
        return NULL;
    }
    if (inconsistent_const) {
        for (size_t i = 0; i < nF; i++) gb_poly_free(F[i]);
        free(F);
        free(all_vars);
        for (size_t i = 0; i < n_eq; i++) expr_free(polys[i]);
        free(polys);
        varset_free(&mains);
        if (wrap_vars) expr_free(wrap_vars);
        free(alg_elim_array);
        alg_free(&alg);
        return expr_new_symbol("False");
    }
    if (nF == 0) {
        /* All inputs collapsed to 0 -> no constraint remains. */
        free(F);
        free(all_vars);
        for (size_t i = 0; i < n_eq; i++) expr_free(polys[i]);
        free(polys);
        varset_free(&mains);
        if (wrap_vars) expr_free(wrap_vars);
        free(alg_elim_array);
        alg_free(&alg);
        return expr_new_symbol("True");
    }

    /* ----- Buchberger ----- */
    size_t out_n = 0;
    GBPoly** G = gb_buchberger(F, nF, &out_n);
    for (size_t i = 0; i < nF; i++) gb_poly_free(F[i]);
    free(F);

    /* ----- Elimination filter: drop polys that still mention any elim var ----- */
    int* elim_idx = (int*)malloc(sizeof(int) * (n_elim > 0 ? n_elim : 1));
    for (size_t i = 0; i < n_elim; i++) elim_idx[i] = (int)i;
    size_t k = 0;
    for (size_t i = 0; i < out_n; i++) {
        if (gb_poly_free_of_vars(G[i], elim_idx, (int)n_elim)) {
            G[k++] = G[i];
        } else {
            gb_poly_free(G[i]);
        }
    }
    out_n = k;
    free(elim_idx);

    /* ----- Strip extraneous monomial factors from each surviving poly -----
     * After cross-multiplying through the algebraisation denominators
     * the elimination ideal often contains a factor of a main variable
     * (e.g. `u * (...)` where the `(...)` is the primitive form
     * Mathematica reports).  Divide it out for cosmetic / downstream
     * agreement.  Cheap: O(sum of term * n_vars). */
    for (size_t i = 0; i < out_n; i++) {
        gb_poly_strip_monomial(G[i]);
    }

    /* ----- Sentinels: empty -> True; contains constant -> False ----- */
    if (out_n == 0) {
        gb_basis_free(G, 0);
        free(all_vars);
        for (size_t i = 0; i < n_eq; i++) expr_free(polys[i]);
        free(polys);
        varset_free(&mains);
        if (wrap_vars) expr_free(wrap_vars);
        free(alg_elim_array);
        alg_free(&alg);
        return expr_new_symbol("True");
    }
    for (size_t i = 0; i < out_n; i++) {
        if (gb_poly_is_constant(G[i]) && !gb_poly_is_zero(G[i])) {
            gb_basis_free(G, out_n);
            free(all_vars);
            for (size_t i2 = 0; i2 < n_eq; i2++) expr_free(polys[i2]);
            free(polys);
            varset_free(&mains);
            if (wrap_vars) expr_free(wrap_vars);
            free(alg_elim_array);
            alg_free(&alg);
            return expr_new_symbol("False");
        }
    }

    /* ----- Convert each surviving basis poly to a balanced Equal[] ----- */
    Expr** eqs = (Expr**)malloc(sizeof(Expr*) * out_n);
    for (size_t i = 0; i < out_n; i++) {
        eqs[i] = balance_polynomial(G[i], all_vars);
    }
    gb_basis_free(G, out_n);
    free(all_vars);
    for (size_t i = 0; i < n_eq; i++) expr_free(polys[i]);
    free(polys);
    varset_free(&mains);
    if (wrap_vars) expr_free(wrap_vars);
    free(alg_elim_array);
    alg_free(&alg);

    /* ----- Combine ----- */
    if (out_n == 1) {
        Expr* out = eqs[0];
        free(eqs);
        return out;
    }
    Expr* and_e = expr_new_function(expr_new_symbol("And"), eqs, out_n);
    free(eqs);
    Expr* evd = evaluate(and_e);
    expr_free(and_e);
    return evd;
}

/* ------------------------------------------------------------------ */
/*  Initialisation                                                     */
/* ------------------------------------------------------------------ */

void eliminate_init(void) {
    symtab_add_builtin("Eliminate", builtin_eliminate);
    SymbolDef* def = symtab_get_def("Eliminate");
    if (def) def->attributes |= ATTR_PROTECTED;
}
