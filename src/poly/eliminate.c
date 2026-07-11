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
#include "trigsimp.h"
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
    Expr* out = expr_new_function(expr_new_symbol(SYM_Equal),
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
                Expr* r = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
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
    Expr* pow_e = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
    free(pa);

    Expr** ea = (Expr**)malloc(sizeof(Expr*) * 2);
    ea[0] = pow_e;
    ea[1] = alg_substitute(a->bases[i], a);
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal), ea, 2);
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
/*  Transcendental algebraisation pre-pass                             */
/* ------------------------------------------------------------------ */

/* The radical pass above handles `Power[base, p/q]`.  This pass is its
 * transcendental sibling: when an elim variable sits inside a circular /
 * hyperbolic trig function, an exponential (`b^x`), or a logarithm, we
 * replace each kernel with fresh aux symbols and push them onto the elim
 * list.  Circular/hyperbolic kernels contribute a `Sin`/`Cos` aux pair
 * plus the Pythagorean constraint `s^2 + c^2 == 1` (circular) /
 * `c^2 - s^2 == 1` (hyperbolic); Exp and Log kernels contribute a single
 * algebraically-free aux.  After substitution the system is polynomial in
 * the (extended) variable set, exactly as the radical pass leaves it, so
 * the rest of the pipeline is unchanged.
 *
 * An expansion step makes this fully general rather than a special case:
 *   - trig: reciprocal heads -> Sin/Cos, then `TrigExpand` reduces every
 *     multiple/sum angle (`Sin[2x]`, `Sin[x+y]`) to a polynomial in
 *     Sin/Cos of *atomic* angles;
 *   - exp:  `b^(p+q) -> b^p b^q`, `b^(k m) -> (b^m)^k`;
 *   - log:  `Log[a b] -> Log[a]+Log[b]`, `Log[a^n] -> n Log[a]`.
 * So `Sin[x]` and `Sin[3x]`, or `Exp[x]` and `Exp[2x]`, collapse onto the
 * same atomic kernel and stay correctly related through one shared aux.
 *
 * Soundness is preserved by two conservative gates (kernel_gate_a and the
 * Gate-B check in builtin_eliminate); when either fires we bail with
 * `nlin` and return the input unevaluated rather than emit an unsound
 * relation. */

typedef enum {
    TK_NONE = 0, TK_SIN, TK_COS, TK_TAN, TK_COT, TK_SEC, TK_CSC
} TrigWhich;

/* Classify a head symbol as a circular/hyperbolic trig function.
 * `*hyper` is set true for the hyperbolic family. */
static TrigWhich trig_which(const char* h, bool* hyper) {
    *hyper = false;
    if (h == SYM_Sin) return TK_SIN;
    if (h == SYM_Cos) return TK_COS;
    if (h == SYM_Tan) return TK_TAN;
    if (h == SYM_Cot) return TK_COT;
    if (h == SYM_Sec) return TK_SEC;
    if (h == SYM_Csc) return TK_CSC;
    *hyper = true;
    if (h == SYM_Sinh) return TK_SIN;
    if (h == SYM_Cosh) return TK_COS;
    if (h == SYM_Tanh) return TK_TAN;
    if (h == SYM_Coth) return TK_COT;
    if (h == SYM_Sech) return TK_SEC;
    if (h == SYM_Csch) return TK_CSC;
    *hyper = false;
    return TK_NONE;
}

/* If `e` is a single-argument trig call, return its classification and
 * store the argument in `*arg_out` (borrowed).  Else TK_NONE. */
static TrigWhich trig_call(const Expr* e, bool* hyper, Expr** arg_out) {
    if (!e || e->type != EXPR_FUNCTION) return TK_NONE;
    if (!e->data.function.head
     || e->data.function.head->type != EXPR_SYMBOL) return TK_NONE;
    if (e->data.function.arg_count != 1) return TK_NONE;
    TrigWhich w = trig_which(e->data.function.head->data.symbol, hyper);
    if (w != TK_NONE) *arg_out = e->data.function.args[0];
    return w;
}

/* Exponential kernel: `Power[base, exp]` with `base` elim-free and `exp`
 * mentioning an elim variable (covers `E^x`, `2^x`, `a^x`).  This is
 * disjoint from the radical pass, which fires only when the *base*
 * contains an elim variable and the exponent is rational. */
static bool exp_kernel(const Expr* e, Expr** elim, size_t n_elim,
                       Expr** base_out, Expr** exp_out) {
    if (!head_is(e, SYM_Power) || e->data.function.arg_count != 2) return false;
    Expr* base = e->data.function.args[0];
    Expr* exp  = e->data.function.args[1];
    if (contains_any_var(base, elim, n_elim)) return false;
    if (!contains_any_var(exp, elim, n_elim)) return false;
    if (base_out) *base_out = base;
    if (exp_out)  *exp_out  = exp;
    return true;
}

/* Logarithmic kernel: `Log[arg]` with `arg` mentioning an elim variable. */
static bool log_kernel(const Expr* e, Expr** elim, size_t n_elim,
                       Expr** arg_out) {
    if (!head_is(e, SYM_Log) || e->data.function.arg_count != 1) return false;
    Expr* arg = e->data.function.args[0];
    if (!contains_any_var(arg, elim, n_elim)) return false;
    if (arg_out) *arg_out = arg;
    return true;
}

/* True iff any transcendental kernel (circular/hyperbolic trig, Exp, or
 * Log) over an elim-containing argument appears in `e`.  This is the
 * cheap gate that keeps polynomial / radical systems on exactly the same
 * code path as before (no transcendental => no overhead). */
static bool has_transcendental_with_elim(const Expr* e,
                                         Expr** elim, size_t n_elim) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    bool hyper;
    Expr* arg = NULL;
    if (trig_call(e, &hyper, &arg) != TK_NONE
        && contains_any_var(arg, elim, n_elim)) {
        return true;
    }
    if (exp_kernel(e, elim, n_elim, NULL, NULL)) return true;
    if (log_kernel(e, elim, n_elim, NULL)) return true;
    if (has_transcendental_with_elim(e->data.function.head, elim, n_elim)) {
        return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_transcendental_with_elim(e->data.function.args[i],
                                         elim, n_elim)) {
            return true;
        }
    }
    return false;
}

/* Build `Sin[theta]` / `Cos[theta]` (or the hyperbolic variants). */
static Expr* make_sin(Expr* arg, bool hyper) {
    return expr_new_function(expr_new_symbol(hyper ? "Sinh" : "Sin"),
                             (Expr*[]){ arg }, 1);
}
static Expr* make_cos(Expr* arg, bool hyper) {
    return expr_new_function(expr_new_symbol(hyper ? "Cosh" : "Cos"),
                             (Expr*[]){ arg }, 1);
}
/* `Power[base, -1]` (a reciprocal). */
static Expr* make_recip(Expr* base) {
    return expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ base, expr_new_integer(-1) }, 2);
}

/* Rewrite Tan/Cot/Sec/Csc (and hyperbolic kin) of an elim-containing
 * argument into Sin/Cos, so the whole system is expressed through the
 * Sin/Cos algebra before TrigExpand and substitution.  Returns a fresh
 * tree the caller owns.  Non-matching nodes are copied verbatim. */
static Expr* to_sincos_rewrite(const Expr* e, Expr** elim, size_t n_elim) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    bool hyper;
    Expr* arg = NULL;
    TrigWhich w = trig_call(e, &hyper, &arg);
    if ((w == TK_TAN || w == TK_COT || w == TK_SEC || w == TK_CSC)
        && contains_any_var(arg, elim, n_elim)) {
        Expr* inner = to_sincos_rewrite(arg, elim, n_elim);
        switch (w) {
            case TK_TAN:  /* Sin/Cos */
                return expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ make_sin(expr_copy(inner), hyper),
                               make_recip(make_cos(inner, hyper)) }, 2);
            case TK_COT:  /* Cos/Sin */
                return expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ make_cos(expr_copy(inner), hyper),
                               make_recip(make_sin(inner, hyper)) }, 2);
            case TK_SEC:  /* 1/Cos */
                return make_recip(make_cos(inner, hyper));
            case TK_CSC:  /* 1/Sin */
                return make_recip(make_sin(inner, hyper));
            default: break;  /* unreachable */
        }
    }

    Expr* new_head = to_sincos_rewrite(e->data.function.head, elim, n_elim);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = to_sincos_rewrite(e->data.function.args[i], elim, n_elim);
    }
    Expr* r = expr_new_function(new_head, new_args, n);
    free(new_args);
    return r;
}

/* Run a single equation through `TrigExpand`, reducing multiple- and
 * sum-angle trig to polynomials in Sin/Cos of atomic angles.  TrigExpand
 * threads over Equal, so the result stays an `Equal[...]`.  Returns a
 * fresh tree the caller owns. */
static Expr* trig_expand_eq(const Expr* eq) {
    Expr* call = expr_new_function(expr_new_symbol(SYM_TrigExpand),
        (Expr*[]){ expr_copy((Expr*)eq) }, 1);
    Expr* out = evaluate(call);
    expr_free(call);
    return out;
}

/* `Power[base, exp]` from owned parts. */
static Expr* make_power(Expr* base, Expr* exp) {
    return expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ base, exp }, 2);
}

/* Split `Times[k, rest...]` into its integer literal factor `*k_out` and
 * the product of the remaining factors `*rest_out` (a fresh owned tree;
 * `1` if none remain).  Returns false if there is no integer factor. */
static bool times_split_int(const Expr* e, int64_t* k_out, Expr** rest_out) {
    if (!is_times(e) || e->data.function.arg_count < 2) return false;
    size_t n = e->data.function.arg_count;
    int64_t k = 0;
    bool found = false;
    Expr** rest = (Expr**)malloc(sizeof(Expr*) * n);
    size_t nr = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* f = e->data.function.args[i];
        if (!found && f->type == EXPR_INTEGER) {
            k = f->data.integer; found = true;
        } else {
            rest[nr++] = expr_copy(f);
        }
    }
    if (!found) {
        for (size_t i = 0; i < nr; i++) expr_free(rest[i]);
        free(rest);
        return false;
    }
    Expr* rest_e;
    if (nr == 0)      { rest_e = expr_new_integer(1); }
    else if (nr == 1) { rest_e = rest[0]; }
    else { rest_e = expr_new_function(expr_new_symbol(SYM_Times), rest, nr); }
    free(rest);  /* frees the array; rest[0] is kept as rest_e when nr==1 */
    *k_out = k;
    *rest_out = rest_e;
    return true;
}

/* Algebraically expand exponential kernels so multiple/sum exponents land
 * on a common atomic exponent:  b^(p+q) -> b^p b^q,  b^(k m) -> (b^m)^k
 * (k an integer).  Only fires on `Power[base, exp]` with elim-free base
 * and elim-bearing exponent; everything else is copied.  No evaluator
 * call is made (which would recombine the split powers), so the caller
 * must substitute before any re-evaluation. */
static Expr* exp_expand(const Expr* e, Expr** elim, size_t n_elim) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    Expr* base = NULL;
    Expr* exp  = NULL;
    if (exp_kernel(e, elim, n_elim, &base, &exp)) {
        if (is_plus(exp) && exp->data.function.arg_count >= 2) {
            size_t n = exp->data.function.arg_count;
            Expr** factors = (Expr**)malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) {
                Expr* term = make_power(expr_copy(base),
                                        expr_copy(exp->data.function.args[i]));
                factors[i] = exp_expand(term, elim, n_elim);
                expr_free(term);
            }
            Expr* r = expr_new_function(expr_new_symbol(SYM_Times), factors, n);
            free(factors);
            return r;
        }
        int64_t k;
        Expr* rest = NULL;
        if (times_split_int(exp, &k, &rest) && k != 1
            && contains_any_var(rest, elim, n_elim)) {
            Expr* inner = make_power(expr_copy(base), rest);  /* owns rest */
            Expr* inner_x = exp_expand(inner, elim, n_elim);
            expr_free(inner);
            Expr* r = make_power(inner_x, expr_new_integer(k));
            return r;
        }
        if (rest) expr_free(rest);
        /* Atomic exponent: leave as Power[base, exp]. */
        return expr_copy((Expr*)e);
    }

    Expr* new_head = exp_expand(e->data.function.head, elim, n_elim);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = exp_expand(e->data.function.args[i], elim, n_elim);
    }
    Expr* r = expr_new_function(new_head, new_args, n);
    free(new_args);
    return r;
}

/* Algebraically expand logarithmic kernels onto atomic arguments:
 *   Log[a b] -> Log[a] + Log[b],   Log[a^n] -> n Log[a]   (n an integer).
 * Same contract as exp_expand (no evaluator call). */
static Expr* log_expand(const Expr* e, Expr** elim, size_t n_elim) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    Expr* arg = NULL;
    if (log_kernel(e, elim, n_elim, &arg)) {
        if (is_times(arg) && arg->data.function.arg_count >= 2) {
            size_t n = arg->data.function.arg_count;
            Expr** terms = (Expr**)malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) {
                Expr* lg = expr_new_function(expr_new_symbol(SYM_Log),
                    (Expr*[]){ expr_copy(arg->data.function.args[i]) }, 1);
                terms[i] = log_expand(lg, elim, n_elim);
                expr_free(lg);
            }
            Expr* r = expr_new_function(expr_new_symbol(SYM_Plus), terms, n);
            free(terms);
            return r;
        }
        if (head_is(arg, SYM_Power) && arg->data.function.arg_count == 2
            && arg->data.function.args[1]->type == EXPR_INTEGER) {
            Expr* lg = expr_new_function(expr_new_symbol(SYM_Log),
                (Expr*[]){ expr_copy(arg->data.function.args[0]) }, 1);
            Expr* inner = log_expand(lg, elim, n_elim);
            expr_free(lg);
            Expr* r = expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_copy(arg->data.function.args[1]), inner }, 2);
            return r;
        }
        /* Atomic argument: leave as Log[arg]. */
        return expr_copy((Expr*)e);
    }

    Expr* new_head = log_expand(e->data.function.head, elim, n_elim);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = log_expand(e->data.function.args[i], elim, n_elim);
    }
    Expr* r = expr_new_function(new_head, new_args, n);
    free(new_args);
    return r;
}

/* ------------------------------------------------------------------ */
/*  Commensurate-exponent folding                                      */
/* ------------------------------------------------------------------ */

/* `exp_expand` above splits sum exponents and factors *integer* multiples
 * onto a common atomic kernel (`b^(k m) -> (b^m)^k`).  This pass is its
 * fractional sibling: exponentials whose exponents differ only by a
 * *rational* factor of a shared monomial are commensurate and must collapse
 * onto one atomic kernel, e.g. `E^(x/3)` and `E^(x/2)` are both powers of
 * `E^(x/6)` (`E^(x/3) = (E^(x/6))^2`, `E^(x/2) = (E^(x/6))^3`).  Without
 * this the three kernels register as algebraically-independent auxes and
 * Gate A rejects the system as `nlin`.
 *
 * The atomic exponent for a `(base, monomial)` group is `(1/L) * monomial`
 * where `L` is the LCM of the denominators of every rational coefficient
 * seen for that group; each `base^((p/q) m)` then rewrites to
 * `(base^((1/L) m))^(p L / q)` with an integer outer power (since q | L).
 * This subsumes the integer case (q = 1 => L unaffected), so it composes
 * cleanly with whatever `exp_expand` already produced. */

/* Decompose an exponent into a rational scalar coefficient `num/den`
 * (den > 0) and the remaining monomial (fresh owned tree).  A leading
 * Integer or Rational factor of a `Times` is peeled as the coefficient;
 * anything else yields coefficient 1 and monomial = copy(exp).  Generalises
 * `times_split_int` to fractional coefficients. */
static void exp_split_rational(const Expr* exp, int64_t* num, int64_t* den,
                               Expr** mono_out) {
    *num = 1; *den = 1;
    if (is_times(exp) && exp->data.function.arg_count >= 2) {
        const Expr* lead = exp->data.function.args[0];
        int64_t p, q;
        if (is_rational(lead, &p, &q)) {
            if (q < 0) { q = -q; p = -p; }
            *num = p; *den = q;
            size_t n = exp->data.function.arg_count;
            if (n == 2) { *mono_out = expr_copy(exp->data.function.args[1]); return; }
            Expr** rest = (Expr**)malloc(sizeof(Expr*) * (n - 1));
            for (size_t i = 1; i < n; i++) {
                rest[i - 1] = expr_copy(exp->data.function.args[i]);
            }
            *mono_out = expr_new_function(expr_new_symbol(SYM_Times), rest, n - 1);
            free(rest);
            return;
        }
    }
    *mono_out = expr_copy((Expr*)exp);
}

/* Build the atomic exponent `(1/L) * monomial` as a flat `Times`, or just
 * `monomial` when L == 1.  Deterministic so every rewrite of the same
 * (base, monomial) group yields a structurally identical inner exponent
 * (the invariant the downstream kernel dedup relies on). */
static Expr* make_scaled_monomial(int64_t L, const Expr* mono) {
    if (L == 1) return expr_copy((Expr*)mono);
    Expr* coeff = expr_new_function(expr_new_symbol(SYM_Rational),
        (Expr*[]){ expr_new_integer(1), expr_new_integer(L) }, 2);
    if (is_times(mono)) {
        size_t n = mono->data.function.arg_count;
        Expr** args = (Expr**)malloc(sizeof(Expr*) * (n + 1));
        args[0] = coeff;
        for (size_t i = 0; i < n; i++) args[i + 1] = expr_copy(mono->data.function.args[i]);
        Expr* r = expr_new_function(expr_new_symbol(SYM_Times), args, n + 1);
        free(args);
        return r;
    }
    return expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ coeff, expr_copy((Expr*)mono) }, 2);
}

/* Per-(base, monomial) LCM of exponent-coefficient denominators. */
typedef struct {
    Expr**   base;
    Expr**   mono;
    int64_t* lcms;
    size_t   n;
    size_t   cap;
} ExpLcmState;

static void explcm_init(ExpLcmState* s) {
    s->base = NULL; s->mono = NULL; s->lcms = NULL; s->n = 0; s->cap = 0;
}

static void explcm_free(ExpLcmState* s) {
    if (s->base) for (size_t i = 0; i < s->n; i++) expr_free(s->base[i]);
    if (s->mono) for (size_t i = 0; i < s->n; i++) expr_free(s->mono[i]);
    free(s->base); free(s->mono); free(s->lcms);
    explcm_init(s);
}

static void explcm_add(ExpLcmState* s, const Expr* base, const Expr* mono,
                       int64_t den) {
    if (den < 0) den = -den;
    if (den == 0) den = 1;
    for (size_t i = 0; i < s->n; i++) {
        if (expr_eq(s->base[i], base) && expr_eq(s->mono[i], mono)) {
            s->lcms[i] = lcm(s->lcms[i], den);
            return;
        }
    }
    if (s->n == s->cap) {
        s->cap = s->cap ? s->cap * 2 : 4;
        s->base = (Expr**)realloc(s->base, sizeof(Expr*) * s->cap);
        s->mono = (Expr**)realloc(s->mono, sizeof(Expr*) * s->cap);
        s->lcms = (int64_t*)realloc(s->lcms, sizeof(int64_t) * s->cap);
    }
    s->base[s->n] = expr_copy((Expr*)base);
    s->mono[s->n] = expr_copy((Expr*)mono);
    s->lcms[s->n] = den;
    s->n++;
}

/* Walk `e`, registering the coefficient denominator of every exponential
 * kernel `base^(coeff*monomial)` under its (base, monomial) group. */
static void explcm_collect(const Expr* e, ExpLcmState* s,
                           Expr** elim, size_t n_elim) {
    if (!e || e->type != EXPR_FUNCTION) return;
    Expr* base = NULL;
    Expr* exp  = NULL;
    if (exp_kernel(e, elim, n_elim, &base, &exp)) {
        int64_t num, den;
        Expr* mono = NULL;
        exp_split_rational(exp, &num, &den, &mono);
        explcm_add(s, base, mono, den);
        expr_free(mono);
    }
    if (e->data.function.head) explcm_collect(e->data.function.head, s, elim, n_elim);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        explcm_collect(e->data.function.args[i], s, elim, n_elim);
    }
}

/* Rewrite each exponential kernel `base^((p/q) m)` to `(base^((1/L) m))^k`
 * with `k = p L / q` an integer.  Returns a fresh tree the caller owns. */
static Expr* exp_commensurate(const Expr* e, const ExpLcmState* s,
                              Expr** elim, size_t n_elim) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    Expr* base = NULL;
    Expr* exp  = NULL;
    if (exp_kernel(e, elim, n_elim, &base, &exp)) {
        int64_t num, den;
        Expr* mono = NULL;
        exp_split_rational(exp, &num, &den, &mono);
        int64_t L = den;  /* den itself divides the group LCM */
        for (size_t i = 0; i < s->n; i++) {
            if (expr_eq(s->base[i], base) && expr_eq(s->mono[i], mono)) {
                L = s->lcms[i];
                break;
            }
        }
        Expr* inner_exp = make_scaled_monomial(L, mono);
        expr_free(mono);
        Expr* inner = make_power(expr_copy(base), inner_exp);
        int64_t k = (num * L) / den;  /* den | L => exact */
        if (k == 1) return inner;
        return make_power(inner, expr_new_integer(k));
    }

    Expr* new_head = exp_commensurate(e->data.function.head, s, elim, n_elim);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = exp_commensurate(e->data.function.args[i], s, elim, n_elim);
    }
    Expr* r = expr_new_function(new_head, new_args, n);
    free(new_args);
    return r;
}

/* Kernel families handled by the transcendental algebraisation pass.
 * CIRC / HYP carry an aux pair (Sin/Cos) plus a Pythagorean constraint;
 * EXP / LOG carry a single aux and no constraint. */
typedef enum { K_CIRC = 0, K_HYP, K_EXP, K_LOG } KernelKind;

/* One distinct transcendental kernel group. */
typedef struct {
    int*     kind;   /* KernelKind                                       */
    Expr**   key;    /* CIRC/HYP: angle theta; EXP: base; LOG: log-arg   */
    Expr**   key2;   /* EXP: atomic exponent monomial; else NULL         */
    Expr**   aux1;   /* owned primary aux (Sin/Sinh, Exp-aux, Log-aux)   */
    Expr**   aux2;   /* owned secondary aux (Cos/Cosh); NULL for EXP/LOG */
    size_t   n;
    size_t   cap;
} KernelState;

static void kernel_init(KernelState* t) {
    t->kind = NULL; t->key = NULL; t->key2 = NULL;
    t->aux1 = NULL; t->aux2 = NULL;
    t->n = 0; t->cap = 0;
}

static void kernel_free(KernelState* t) {
    if (t->key)  for (size_t i = 0; i < t->n; i++) expr_free(t->key[i]);
    if (t->key2) for (size_t i = 0; i < t->n; i++) expr_free(t->key2[i]);
    if (t->aux1) for (size_t i = 0; i < t->n; i++) expr_free(t->aux1[i]);
    if (t->aux2) for (size_t i = 0; i < t->n; i++) expr_free(t->aux2[i]);
    free(t->kind); free(t->key); free(t->key2); free(t->aux1); free(t->aux2);
    kernel_init(t);
}

/* The elim-bearing part of group `g` (what the soundness gates test):
 * the exponent for EXP kernels, the key (angle / log-arg) otherwise. */
static Expr* kernel_elim_part(const KernelState* t, size_t g) {
    return t->key2[g] ? t->key2[g] : t->key[g];
}

/* Register a kernel group, deduping by (kind, key, key2) under structural
 * equality.  `key2` may be NULL.  Aux symbols are filled in later. */
static void kernel_add(KernelState* t, KernelKind kind,
                       const Expr* key, const Expr* key2) {
    for (size_t i = 0; i < t->n; i++) {
        if (t->kind[i] != (int)kind) continue;
        if (!expr_eq(t->key[i], key)) continue;
        bool k2eq = (!t->key2[i] && !key2)
                 || (t->key2[i] && key2 && expr_eq(t->key2[i], key2));
        if (k2eq) return;
    }
    if (t->n == t->cap) {
        t->cap = t->cap ? t->cap * 2 : 4;
        t->kind = (int*)  realloc(t->kind, sizeof(int)   * t->cap);
        t->key  = (Expr**)realloc(t->key,  sizeof(Expr*) * t->cap);
        t->key2 = (Expr**)realloc(t->key2, sizeof(Expr*) * t->cap);
        t->aux1 = (Expr**)realloc(t->aux1, sizeof(Expr*) * t->cap);
        t->aux2 = (Expr**)realloc(t->aux2, sizeof(Expr*) * t->cap);
    }
    t->kind[t->n] = (int)kind;
    t->key[t->n]  = expr_copy((Expr*)key);
    t->key2[t->n] = key2 ? expr_copy((Expr*)key2) : NULL;
    t->aux1[t->n] = NULL;
    t->aux2[t->n] = NULL;
    t->n++;
}

/* Walk `e`, registering every transcendental kernel whose elim-bearing
 * part mentions an elim variable.  Runs AFTER the expansion steps, so
 * trig is over atomic angles, exp over atomic exponents, and log over
 * atomic arguments (any surviving compound form shares an elim var
 * across distinct groups and is rejected by Gate A). */
static void kernel_collect(const Expr* e, KernelState* t,
                           Expr** elim, size_t n_elim) {
    if (!e || e->type != EXPR_FUNCTION) return;
    bool hyper;
    Expr* arg = NULL;
    Expr* base = NULL;
    Expr* exp = NULL;
    TrigWhich w = trig_call(e, &hyper, &arg);
    if (w != TK_NONE && contains_any_var(arg, elim, n_elim)) {
        kernel_add(t, hyper ? K_HYP : K_CIRC, arg, NULL);
    } else if (exp_kernel(e, elim, n_elim, &base, &exp)) {
        kernel_add(t, K_EXP, base, exp);
    } else if (log_kernel(e, elim, n_elim, &arg)) {
        kernel_add(t, K_LOG, arg, NULL);
    }
    kernel_collect(e->data.function.head, t, elim, n_elim);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        kernel_collect(e->data.function.args[i], t, elim, n_elim);
    }
}

/* Gate A: an elim variable that appears inside two *distinct* kernel
 * groups cannot be captured (e.g. `Sin[x]` and `Sin[x*y]` share `x`, but
 * a product argument has no multiple-angle relation to `x`; likewise
 * `Sin[x]` and `Sinh[x]`, or `Sin[x]` and `Exp[x]`, are algebraically
 * independent).  Returns true if the algebraisation would be unsound. */
static bool kernel_gate_a(const KernelState* t, Expr** elim, size_t n_elim) {
    for (size_t e = 0; e < n_elim; e++) {
        int cnt = 0;
        for (size_t g = 0; g < t->n; g++) {
            if (contains_any_var(kernel_elim_part(t, g), &elim[e], 1)) cnt++;
            if (cnt > 1) return true;
        }
    }
    return false;
}

/* True iff `target` occurs as a genuine polynomial atom in `e` -- i.e.
 * reachable from the root through only the structural poly-recursive
 * heads (Plus/Times/...) and integer Power bases, NOT buried inside a
 * function-shaped atom (which the Groebner engine treats as a single
 * opaque variable).  Mirrors `collect_main_vars`' traversal so the two
 * agree on what counts as a live variable.
 *
 * Gate B uses this to distinguish a real free occurrence of a trig-bound
 * variable `x` (e.g. bare `x` in `x + y == 2`, which severs the lost
 * `x <-> Sin[x]` link) from a harmless one inside another elim atom
 * (e.g. `Dt[x]`, which the engine eliminates as a unit). */
static bool poly_atom_occurs(const Expr* e, const Expr* target) {
    if (!e) return false;
    if (expr_eq(e, target)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    const char* hsym = (e->data.function.head
                        && e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol : NULL;
    if (hsym && is_poly_recursive_head(hsym)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (poly_atom_occurs(e->data.function.args[i], target)) return true;
        }
        return false;
    }
    if (hsym == SYM_Power && e->data.function.arg_count == 2
        && e->data.function.args[1]->type == EXPR_INTEGER) {
        return poly_atom_occurs(e->data.function.args[0], target);
    }
    /* Function-shaped atom not equal to target: opaque, do not recurse. */
    return false;
}

/* Substitute every registered trig call (matching head family + atomic
 * angle) with the corresponding rational expression in the aux symbols:
 *   Sin->s   Cos->c   Tan->s/c   Cot->c/s   Sec->1/c   Csc->1/s
 * The 1/c, 1/s denominators are cleared downstream by the
 * Together->Numerator step in `normalise_equation`.  Handling every
 * family here (not just Sin/Cos) is essential: `trig_expand_eq` runs the
 * evaluator, which re-canonicalises `Sin/Cos` back into `Tan`, so a Tan
 * kernel can reappear after the to_sincos rewrite.  Returns a fresh
 * tree. */
static Expr* kernel_substitute(const Expr* e, const KernelState* t) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    /* Circular / hyperbolic trig: match head family + atomic angle. */
    bool hyper;
    Expr* arg = NULL;
    TrigWhich w = trig_call(e, &hyper, &arg);
    if (w != TK_NONE) {
        int want = hyper ? K_HYP : K_CIRC;
        for (size_t g = 0; g < t->n; g++) {
            if (t->kind[g] != want) continue;
            if (!expr_eq(t->key[g], arg)) continue;
            Expr* s = t->aux1[g];
            Expr* c = t->aux2[g];
            switch (w) {
                case TK_SIN: return expr_copy(s);
                case TK_COS: return expr_copy(c);
                case TK_TAN: return expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_copy(s), make_recip(expr_copy(c)) }, 2);
                case TK_COT: return expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_copy(c), make_recip(expr_copy(s)) }, 2);
                case TK_SEC: return make_recip(expr_copy(c));
                case TK_CSC: return make_recip(expr_copy(s));
                default: break;
            }
        }
    }

    /* Exponential: match Power[base, exponent] structurally. */
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* b = e->data.function.args[0];
        Expr* x = e->data.function.args[1];
        for (size_t g = 0; g < t->n; g++) {
            if (t->kind[g] != K_EXP) continue;
            if (expr_eq(t->key[g], b) && expr_eq(t->key2[g], x)) {
                return expr_copy(t->aux1[g]);
            }
        }
    }

    /* Logarithm: match Log[arg]. */
    if (head_is(e, SYM_Log) && e->data.function.arg_count == 1) {
        Expr* a = e->data.function.args[0];
        for (size_t g = 0; g < t->n; g++) {
            if (t->kind[g] != K_LOG) continue;
            if (expr_eq(t->key[g], a)) return expr_copy(t->aux1[g]);
        }
    }

    Expr* new_head = kernel_substitute(e->data.function.head, t);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = kernel_substitute(e->data.function.args[i], t);
    }
    Expr* r = expr_new_function(new_head, new_args, n);
    free(new_args);
    return r;
}

/* The algebraic constraint for group `i`, or NULL when the kernel needs
 * none (Exp/Log auxes are algebraically free):
 *   circular:   s^2 + c^2 == 1
 *   hyperbolic: c^2 - s^2 == 1   (cosh^2 - sinh^2 == 1) */
static Expr* kernel_make_constraint(const KernelState* t, size_t i) {
    if (t->kind[i] != K_CIRC && t->kind[i] != K_HYP) return NULL;
    Expr* s2 = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_copy(t->aux1[i]), expr_new_integer(2) }, 2);
    Expr* c2 = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_copy(t->aux2[i]), expr_new_integer(2) }, 2);
    Expr* sum;
    if (t->kind[i] == K_HYP) {
        Expr* neg_s2 = expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_new_integer(-1), s2 }, 2);
        sum = expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ c2, neg_s2 }, 2);
    } else {
        sum = expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ s2, c2 }, 2);
    }
    return expr_new_function(expr_new_symbol(SYM_Equal),
        (Expr*[]){ sum, expr_new_integer(1) }, 2);
}

/* Pick an aux symbol name `<prefix><idx>$` not colliding with anything in
 * `eqs` or the already-assigned auxes. */
static const char* kernel_fresh_name(const char* prefix, Expr** eqs,
                                     size_t n_eq, const KernelState* t,
                                     size_t i_so_far, int* next_idx) {
    char buf[64];
    for (int tries = 0; tries < 10000; tries++) {
        snprintf(buf, sizeof(buf), "%s%d$", prefix, *next_idx);
        (*next_idx)++;
        bool collide = false;
        for (size_t j = 0; j < n_eq && !collide; j++) {
            if (walk_has_symbol_name(eqs[j], buf)) collide = true;
        }
        for (size_t j = 0; j < i_so_far && !collide; j++) {
            if (t->aux1[j] && walk_has_symbol_name(t->aux1[j], buf)) collide = true;
            if (t->aux2[j] && walk_has_symbol_name(t->aux2[j], buf)) collide = true;
        }
        if (!collide) return intern_symbol(buf);
    }
    return intern_symbol("$t$overflow$");
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
        return expr_new_function(expr_new_symbol(SYM_Rational),
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
            Expr* r = expr_new_function(expr_new_symbol(SYM_Times), args, n - 1);
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
            Expr* r = expr_new_function(expr_new_symbol(SYM_Times), args, n);
            free(args);
            return r;
        }
        /* No coefficient -- prepend -1. */
        Expr** args = (Expr**)malloc(sizeof(Expr*) * (n + 1));
        args[0] = expr_new_integer(-1);
        for (size_t i = 0; i < n; i++) {
            args[i + 1] = expr_copy(t->data.function.args[i]);
        }
        Expr* r = expr_new_function(expr_new_symbol(SYM_Times), args, n + 1);
        free(args);
        return r;
    }
    /* Symbol, Power, or anything else: wrap as Times[-1, t]. */
    Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
    args[0] = expr_new_integer(-1);
    args[1] = expr_copy((Expr*)t);
    Expr* r = expr_new_function(expr_new_symbol(SYM_Times), args, 2);
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
        lhs = expr_new_function(expr_new_symbol(SYM_Plus), pos, n_pos);
    }
    free(pos);

    Expr* rhs;
    if (n_neg == 0) {
        rhs = expr_new_integer(0);
    } else if (n_neg == 1) {
        rhs = neg[0];
    } else {
        rhs = expr_new_function(expr_new_symbol(SYM_Plus), neg, n_neg);
    }
    free(neg);

    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
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

/* See eliminate.h: internal drivers (DerivativeDivides) raise this to mute
 * the advisory diagnostics for an elimination the user never requested. */
int eliminate_suppress_messages = 0;

static void emit_argt(size_t argc) {
    if (eliminate_suppress_messages) return;
    fprintf(stderr, "Eliminate::argt: Eliminate called with %zu argument(s); "
                    "2 expected.\n", argc);
}

static void emit_eqf(void) {
    if (eliminate_suppress_messages) return;
    fprintf(stderr, "Eliminate::eqf: equations must be given as Equal[lhs, rhs] "
                    "(==) or a list/And of such.\n");
}

static void emit_nlin(void) {
    if (eliminate_suppress_messages) return;
    fprintf(stderr, "Eliminate::nlin: system is not polynomial in the given "
                    "variables.\n");
}

static void emit_ifun(void) {
    if (eliminate_suppress_messages) return;
    fprintf(stderr, "Eliminate::ifun: Inverse functions are being used by "
                    "Eliminate, so some solutions may not be found; use Reduce "
                    "for complete solution information.\n");
}

/* Emit at most once per distinct input expression.  `RepeatedTiming`,
 * `Table[...]`, etc. evaluate the same `Eliminate[...]` form hundreds
 * of times; without the hash guard each iteration re-prints the warning
 * and floods the terminal.  Mirrors `solve.c:warn_bad_option`. */
static void emit_alg(const Expr* res) {
    if (eliminate_suppress_messages) return;
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
/*  Inverse-function substitution pre-pass                             */
/* ------------------------------------------------------------------ */

/* The forward transcendental pass (Sin/Cos/Exp/Log of an *elim* variable)
 * and the top-level inverse-rewrite (`f[u]==v` -> `u==f^-1[v]`) do not
 * cover the u-substitution shape that arises when integrating by parts /
 * by substitution, e.g.
 *
 *   Eliminate[{Dt[y] == x ArcSin[x]/Sqrt[1-x^2] Dt[x],
 *              u == ArcSin[x],
 *              Dt[u] == Dt[x]/Sqrt[1-x^2]}, {x, Dt[x]}]
 *
 * Here `ArcSin[x]` sits *inside a product* in the first equation, so it is
 * an opaque function-atom still mentioning the elim variable `x` -> the
 * pipeline bails `nlin`.  Yet the second equation *defines* the inverse
 * function: `u == ArcSin[x]`.
 *
 * This pass exploits such a defining equation.  Given `M == ArcF[x]` with
 * `M` elim-free and `x` a single elim symbol, it rewrites the whole system
 * so the elim variable appears only polynomially:
 *
 *   (1) replace every `ArcF[x]` by the main expression `M` (from the
 *       defining equation), everywhere;
 *   (2) replace the companion radical `(1 -/+ x^2)^(p/2)` (or `(x^2-1)^(p/2)`)
 *       by `coFn[M]^(co_sign*p)`, where `coFn` is the trigonometric
 *       co-function of `ArcF` -- e.g. `Sqrt[1-x^2] -> Cos[M]` for ArcSin.
 *       Because `coFn[M]` is a trig call of the *main* variable `M`, it is
 *       an ordinary (non-eliminated) polynomial atom, so any factor of it
 *       left over after elimination is divided out by
 *       `gb_poly_strip_monomial`.  (Algebraising the radical as an
 *       independent aux instead would leave a squared / branch-ambiguous
 *       factor and a messy answer.)
 *   (3) replace the defining equation with `x == F[M]` (F the forward
 *       function), tying the elim variable to a main-variable trig atom.
 *
 * The result is exactly the hand-substituted system a human would write
 * before running the Groebner elimination, and it reduces to Mathematica's
 * `Eliminate::ifun`-flagged answer.  Principal branches are assumed (hence
 * the `ifun` advisory), matching the forward pass's soundness contract. */

/* One inverse-function family. `arc`/`fwd`/`co` are interned head symbols;
 * `has_comp` gates the companion-radical rewrite; `const_term`/`x2_sign`
 * build the companion base `const_term + x2_sign*x^2`. */
typedef struct {
    const char* fwd;      /* forward function head (Sin/Cos/.../Exp)       */
    const char* co;       /* co-function head for the companion radical    */
    bool        has_comp; /* false for Log (no companion radical)          */
    int         const_term; /* +1 / -1  : the constant of the comp base    */
    int         x2_sign;    /* +1 / -1  : sign of x^2 in the comp base     */
    int         co_sign;    /* +1 / -1  : radical == co[M]^(co_sign*p)      */
} InvEntry;

/* Map an inverse-function head symbol to its family record. Uses the
 * runtime-interned SYM_* pointers (valid once core_init() has run), so it
 * cannot be a compile-time table. Returns false for non-inverse heads. */
static bool inv_lookup(const char* h, InvEntry* out) {
    if (h == SYM_ArcSin)  { *out = (InvEntry){ SYM_Sin,  SYM_Cos,  true,  1, -1, +1 }; return true; }
    if (h == SYM_ArcCos)  { *out = (InvEntry){ SYM_Cos,  SYM_Sin,  true,  1, -1, +1 }; return true; }
    if (h == SYM_ArcTan)  { *out = (InvEntry){ SYM_Tan,  SYM_Cos,  true,  1, +1, -1 }; return true; }
    if (h == SYM_ArcSinh) { *out = (InvEntry){ SYM_Sinh, SYM_Cosh, true,  1, +1, +1 }; return true; }
    if (h == SYM_ArcCosh) { *out = (InvEntry){ SYM_Cosh, SYM_Sinh, true, -1, +1, +1 }; return true; }
    if (h == SYM_ArcTanh) { *out = (InvEntry){ SYM_Tanh, SYM_Cosh, true,  1, -1, -1 }; return true; }
    /* Log: forward function is Exp, no companion radical.  Enabled only as a
     * *fallback* — `inv_detect` restricts a `M == Log[x]` defining equation
     * to the case where `x` also occurs polynomially (so the forward Log
     * algebraisation alone would trip Gate B).  When `x` appears solely
     * inside logs, the forward pass gives a tidier answer, so we defer to it.
     * Pinning `x == E^M` (rather than substituting `x -> E^M`) keeps `E^M`
     * as a single main-variable atom, sidestepping the exponential-blowup the
     * naive substitution would cause. */
    if (h == SYM_Log)     { *out = (InvEntry){ SYM_Exp,  NULL,     false, 0,  0,  0 }; return true; }
    return false;
}

/* Build the (evaluated, canonical) companion base `const_term + x2_sign*x^2`
 * for the elim symbol `xsym`, ready for `expr_eq` matching against radical
 * bases in the tree. */
static Expr* inv_build_comp_base(const char* xsym, const InvEntry* e) {
    Expr* x2 = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_new_symbol(xsym), expr_new_integer(2) }, 2);
    Expr* x2term = (e->x2_sign < 0)
        ? expr_new_function(expr_new_symbol(SYM_Times),
              (Expr*[]){ expr_new_integer(-1), x2 }, 2)
        : x2;
    Expr* base = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_new_integer(e->const_term), x2term }, 2);
    Expr* ev = evaluate(base);
    expr_free(base);
    return ev;
}

/* Recursively rewrite `e`: `ArcF[xsym] -> M`, and (when comp_base != NULL)
 * `Power[comp_base, Rational[p,2]] -> co[M]^(co_sign*p)`.  Returns a fresh
 * tree the caller owns.  `arc` is the inverse-function head being matched. */
static Expr* inv_rewrite(const Expr* e, const char* xsym, const Expr* M,
                         const char* arc, const Expr* comp_base,
                         const char* co, int co_sign) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    /* (1) ArcF[x] -> M */
    if (head_is(e, arc) && e->data.function.arg_count == 1) {
        Expr* a = e->data.function.args[0];
        if (a && a->type == EXPR_SYMBOL && a->data.symbol == xsym) {
            return expr_copy((Expr*)M);
        }
    }
    /* (2) companion power (comp_base)^(p/q) -> co[M]^(2 co_sign p / q).
     *
     * The companion base is exactly `co[M]^(2 co_sign)` (e.g. for ArcTan,
     * x = Tan[M] gives 1 + x^2 = Sec[M]^2 = Cos[M]^-2, so co = Cos,
     * co_sign = -1).  A half-integer exponent (q = 2) is the radical case
     * `Sqrt[1-x^2] -> Cos[M]` that ArcSin/ArcSinh/... need; an integer
     * exponent (q = 1) is the *rational-denominator* case `1/(1+x^2) ->
     * Cos[M]^2` that ArcTan/ArcTanh need — without it the `1+x^2` /
     * `1-x^2` factor survives elimination as a spurious `Sec/Sech`
     * multiple.  Handling any q for which `2 co_sign p` divides evenly
     * covers both (and is branch-safe: the identity is a squared one). */
    if (comp_base && co && head_is(e, SYM_Power)
        && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        int64_t p, q;
        if (is_rational(exp, &p, &q) && expr_eq(base, comp_base)) {
            int64_t numer = 2 * (int64_t)co_sign * p;
            if (q != 0 && numer % q == 0) {
                int64_t np = numer / q;
                if (np == 0) return expr_new_integer(1);
                Expr* com = expr_new_function(expr_new_symbol(co),
                    (Expr*[]){ expr_copy((Expr*)M) }, 1);
                if (np == 1) return com;
                return expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ com, expr_new_integer(np) }, 2);
            }
        }
    }

    Expr* new_head = inv_rewrite(e->data.function.head, xsym, M, arc,
                                 comp_base, co, co_sign);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = inv_rewrite(e->data.function.args[i], xsym, M, arc,
                                  comp_base, co, co_sign);
    }
    Expr* r = expr_new_function(new_head, new_args, n);
    free(new_args);
    return r;
}

/* Scan `eqs` for a defining equation `M == ArcF[x]` (either side): ArcF an
 * invertible transcendental, `x` a single elim symbol, `M` elim-free.  On a
 * hit, fill `*fe`, `*arc` (the interned arc head), `*xsym`, and an *owned*
 * copy `*M_out` (the tree it points into is about to be rewritten), then
 * return the equation index.  Returns -1 when none applies. */
static int inv_detect(Expr** eqs, size_t n_eq, Expr** elim, size_t n_elim,
                      InvEntry* fe, const char** arc, const char** xsym,
                      Expr** M_out) {
    for (size_t i = 0; i < n_eq; i++) {
        Expr* eq = eqs[i];
        if (!is_equal(eq) || eq->data.function.arg_count != 2) continue;
        for (int s = 0; s < 2; s++) {
            Expr* side  = eq->data.function.args[s];
            Expr* other = eq->data.function.args[1 - s];
            if (!side || side->type != EXPR_FUNCTION
                || side->data.function.arg_count != 1) continue;
            Expr* h = side->data.function.head;
            if (!h || h->type != EXPR_SYMBOL) continue;
            if (!inv_lookup(h->data.symbol, fe)) continue;
            Expr* a = side->data.function.args[0];
            if (!a || a->type != EXPR_SYMBOL) continue;
            if (!var_in_list(a, elim, n_elim)) continue;
            if (contains_any_var(other, elim, n_elim)) continue;
            /* Log fallback gate: only take the inverse-substitution route
             * when `x` also appears as a genuine polynomial atom somewhere
             * (e.g. `1/x`), which is exactly when the forward Log pass would
             * fail.  Otherwise leave `M == Log[x]` for the forward pass. */
            if (h->data.symbol == SYM_Log) {
                bool poly = false;
                for (size_t j = 0; j < n_eq && !poly; j++) {
                    if (poly_atom_occurs(eqs[j], a)) poly = true;
                }
                if (!poly) continue;
            }
            *arc    = h->data.symbol;
            *xsym   = a->data.symbol;
            *M_out  = expr_copy(other);
            return (int)i;
        }
    }
    return -1;
}

/* Drive the inverse-function substitution to a fixed point over `eqs`
 * (owned array of `n_eq` owned `Equal[]` trees; rewritten in place).  Each
 * round handles one defining equation.  Returns true if any rewrite fired
 * (so the caller emits `ifun`).  `n_eq` is unchanged: each round overwrites
 * exactly the defining-equation slot with `x == F[M]`. */
static bool inv_substitute_all(Expr** eqs, size_t n_eq,
                               Expr** elim, size_t n_elim) {
    bool fired = false;
    /* Bound the loop generously: each round removes one inverse-function
     * kernel from play; the elim count caps how many distinct ones matter. */
    for (size_t round = 0; round < n_elim + n_eq + 4; round++) {
        InvEntry fe;
        const char* arc = NULL;
        const char* xsym = NULL;
        Expr* M = NULL;
        int di = inv_detect(eqs, n_eq, elim, n_elim, &fe, &arc, &xsym, &M);
        if (di < 0) break;
        fired = true;

        Expr* comp = fe.has_comp ? inv_build_comp_base(xsym, &fe) : NULL;

        for (size_t i = 0; i < n_eq; i++) {
            Expr* r = inv_rewrite(eqs[i], xsym, M, arc, comp, fe.co, fe.co_sign);
            expr_free(eqs[i]);
            eqs[i] = r;
        }

        /* Overwrite the defining slot with `x == F[M]` (F the forward fn). */
        Expr* fwd = expr_new_function(expr_new_symbol(fe.fwd),
            (Expr*[]){ expr_copy(M) }, 1);
        Expr* neweq = expr_new_function(expr_new_symbol(SYM_Equal),
            (Expr*[]){ expr_new_symbol(xsym), fwd }, 2);
        expr_free(eqs[di]);
        eqs[di] = neweq;

        if (comp) expr_free(comp);
        expr_free(M);
    }
    return fired;
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
        if (eqns_arg->data.symbol == SYM_True)  return expr_new_symbol(SYM_True);
        if (eqns_arg->data.symbol == SYM_False) return expr_new_symbol(SYM_False);
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
        return expr_new_symbol(SYM_True);
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
                return expr_new_symbol(SYM_False);
            }
        }
        eq_in[n_eq++] = x;
    }
    if (n_eq == 0) {
        free(eq_in);
        return expr_new_symbol(SYM_True);
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
        wrap_vars = expr_new_function(expr_new_symbol(SYM_List), wrapped, 1);
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
        Expr* and_e = expr_new_function(expr_new_symbol(SYM_And), copies, n_eq);
        free(copies);
        Expr* ev = evaluate(and_e);
        expr_free(and_e);
        return ev;
    }

    /* ----- Owned working copies of the equations -----
     * `eq_in` is a borrowed array; the pre-passes below build new trees
     * and may append constraints, so take ownership now. */
    Expr** eq_work = (Expr**)malloc(sizeof(Expr*) * n_eq);
    for (size_t i = 0; i < n_eq; i++) eq_work[i] = expr_copy(eq_in[i]);
    free(eq_in);

    bool fired_ifun = false;
    KernelState ker; kernel_init(&ker);
    Expr** ker_elim_array = NULL;

    /* ----- Inverse-function substitution pre-pass -----
     * When a defining equation `M == ArcF[x]` is present (M elim-free, x a
     * single elim symbol), propagate `ArcF[x] -> M` and the companion
     * radical `Sqrt[1-x^2] -> Cos[M]` (etc.) through the whole system, then
     * pin `x == F[M]`.  This turns u-substitution shapes (with the inverse
     * buried inside a product) into a system that is polynomial in the elim
     * variables, which the forward/radical passes and Buchberger then solve.
     * See the section header above for the full rationale + identity table. */
    if (inv_substitute_all(eq_work, n_eq, elim_items, n_elim)) {
        fired_ifun = true;
    }

    /* ----- Transcendental algebraisation pre-pass -----
     * When an elim variable sits inside a circular/hyperbolic trig, an
     * exponential (`b^x`), or a logarithm, expand multiple/sum/product
     * arguments down to atomic kernels, replace each kernel with fresh
     * aux symbols (Sin/Cos pairs carry the Pythagorean constraint; Exp
     * and Log auxes are algebraically free), and push the auxes onto the
     * elim list.  After this the radical pass and Buchberger see a
     * polynomial system.  The cheap `has_transcendental_with_elim` gate
     * keeps purely polynomial / radical inputs on the identical code path
     * as before (no transcendental => no overhead). */
    bool trans_seen = false;
    for (size_t i = 0; i < n_eq && !trans_seen; i++) {
        if (has_transcendental_with_elim(eq_work[i], elim_items, n_elim)) {
            trans_seen = true;
        }
    }
    if (trans_seen) {
        /* (a) trig: reciprocal heads -> Sin/Cos, then TrigExpand (handles
         *     circular and hyperbolic multiple/sum angles). */
        for (size_t i = 0; i < n_eq; i++) {
            Expr* r = to_sincos_rewrite(eq_work[i], elim_items, n_elim);
            expr_free(eq_work[i]); eq_work[i] = r;
        }
        for (size_t i = 0; i < n_eq; i++) {
            Expr* r = trig_expand_eq(eq_work[i]);
            expr_free(eq_work[i]); eq_work[i] = r;
        }
        /* (b) exp/log: split sum exponents / product logs onto atomic
         *     kernels.  No evaluator runs after this until substitution
         *     (which would recombine the split powers/logs). */
        for (size_t i = 0; i < n_eq; i++) {
            Expr* r = exp_expand(eq_work[i], elim_items, n_elim);
            expr_free(eq_work[i]); eq_work[i] = r;
        }
        for (size_t i = 0; i < n_eq; i++) {
            Expr* r = log_expand(eq_work[i], elim_items, n_elim);
            expr_free(eq_work[i]); eq_work[i] = r;
        }
        /* (b') fold commensurate rational exponents (E^(x/3), E^(x/2))
         *      onto one atomic kernel E^(x/L) so they share a single aux
         *      instead of tripping Gate A as independent groups. */
        {
            ExpLcmState el; explcm_init(&el);
            for (size_t i = 0; i < n_eq; i++) {
                explcm_collect(eq_work[i], &el, elim_items, n_elim);
            }
            if (el.n > 0) {
                for (size_t i = 0; i < n_eq; i++) {
                    Expr* r = exp_commensurate(eq_work[i], &el,
                                               elim_items, n_elim);
                    expr_free(eq_work[i]); eq_work[i] = r;
                }
            }
            explcm_free(&el);
        }
        /* (c) collect distinct kernel groups. */
        for (size_t i = 0; i < n_eq; i++) {
            kernel_collect(eq_work[i], &ker, elim_items, n_elim);
        }
        if (ker.n > 0) {
            /* Gate A: an elim var shared across two distinct kernel groups
             * cannot be captured -> bail to nlin (input unevaluated). */
            if (kernel_gate_a(&ker, elim_items, n_elim)) {
                emit_nlin();
                for (size_t i = 0; i < n_eq; i++) expr_free(eq_work[i]);
                free(eq_work);
                kernel_free(&ker);
                if (wrap_vars) expr_free(wrap_vars);
                return NULL;
            }
            /* Assign fresh aux symbols: one per Exp/Log kernel, a pair per
             * Sin/Cos (circular/hyperbolic) kernel. */
            int idx1 = 1, idx2 = 1;
            size_t n_aux = 0, n_constr = 0;
            for (size_t g = 0; g < ker.n; g++) {
                bool pair = (ker.kind[g] == K_CIRC || ker.kind[g] == K_HYP);
                const char* p1 = (ker.kind[g] == K_EXP) ? "$te"
                               : (ker.kind[g] == K_LOG) ? "$tl" : "$ts";
                ker.aux1[g] = expr_new_symbol(
                    kernel_fresh_name(p1, eq_work, n_eq, &ker, g, &idx1));
                if (pair) {
                    ker.aux2[g] = expr_new_symbol(
                        kernel_fresh_name("$tc", eq_work, n_eq, &ker, g, &idx2));
                    n_aux += 2; n_constr += 1;
                } else {
                    n_aux += 1;
                }
            }
            /* Substitute every kernel -> aux expression everywhere. */
            for (size_t i = 0; i < n_eq; i++) {
                Expr* r = kernel_substitute(eq_work[i], &ker);
                expr_free(eq_work[i]); eq_work[i] = r;
            }
            /* Gate B: a kernel-bound elim var still present as a genuine
             * polynomial atom after substitution (e.g. bare `x` alongside
             * `Sin[x]`) is not captured by the relations -> bail to nlin. */
            bool gate_b = false;
            for (size_t e = 0; e < n_elim && !gate_b; e++) {
                bool bound = false;
                for (size_t g = 0; g < ker.n && !bound; g++) {
                    if (contains_any_var(kernel_elim_part(&ker, g),
                                         &elim_items[e], 1)) {
                        bound = true;
                    }
                }
                if (!bound) continue;
                for (size_t i = 0; i < n_eq && !gate_b; i++) {
                    if (poly_atom_occurs(eq_work[i], elim_items[e])) {
                        gate_b = true;
                    }
                }
            }
            if (gate_b) {
                emit_nlin();
                for (size_t i = 0; i < n_eq; i++) expr_free(eq_work[i]);
                free(eq_work);
                kernel_free(&ker);
                if (wrap_vars) expr_free(wrap_vars);
                return NULL;
            }
            /* Append the constraints (only Sin/Cos kernels have one). */
            eq_work = (Expr**)realloc(eq_work,
                                      sizeof(Expr*) * (n_eq + n_constr));
            size_t ci = 0;
            for (size_t g = 0; g < ker.n; g++) {
                Expr* c = kernel_make_constraint(&ker, g);
                if (c) eq_work[n_eq + ci++] = c;
            }
            /* Extend the elim list with the aux symbols (borrows aux
             * pointers; we own the array and free it at every exit). */
            ker_elim_array =
                (Expr**)malloc(sizeof(Expr*) * (n_elim + n_aux));
            for (size_t i = 0; i < n_elim; i++) {
                ker_elim_array[i] = elim_items[i];
            }
            size_t ai = n_elim;
            for (size_t g = 0; g < ker.n; g++) {
                ker_elim_array[ai++] = ker.aux1[g];
                if (ker.aux2[g]) ker_elim_array[ai++] = ker.aux2[g];
            }
            elim_items = ker_elim_array;
            n_elim    += n_aux;
            n_eq      += n_constr;
            fired_ifun = true;  /* branch/sign info lost: report `ifun` */
        }
    }

    /* ----- Transcendental (inverse-function) pre-pass ----- */
    Expr** rewritten = (Expr**)malloc(sizeof(Expr*) * n_eq);
    for (size_t i = 0; i < n_eq; i++) {
        rewritten[i] = try_inverse_rewrite(eq_work[i], elim_items, n_elim,
                                           &fired_ifun);
    }
    for (size_t i = 0; i < n_eq; i++) expr_free(eq_work[i]);
    free(eq_work);
    if (fired_ifun) emit_ifun();

    /* Validate every rewritten equation is `Equal[lhs, rhs]`. */
    for (size_t i = 0; i < n_eq; i++) {
        if (!is_equal(rewritten[i])
         || rewritten[i]->data.function.arg_count != 2) {
            emit_eqf();
            for (size_t k = 0; k < n_eq; k++) expr_free(rewritten[k]);
            free(rewritten);
            free(ker_elim_array);
            kernel_free(&ker);
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
            free(ker_elim_array);
            kernel_free(&ker);
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
        free(ker_elim_array);
        kernel_free(&ker);
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
        free(ker_elim_array);
        kernel_free(&ker);
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
        free(ker_elim_array);
        kernel_free(&ker);
        return expr_new_symbol(SYM_False);
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
        free(ker_elim_array);
        kernel_free(&ker);
        return expr_new_symbol(SYM_True);
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
        free(ker_elim_array);
        kernel_free(&ker);
        return expr_new_symbol(SYM_True);
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
            free(ker_elim_array);
            kernel_free(&ker);
            return expr_new_symbol(SYM_False);
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
    free(ker_elim_array);
    kernel_free(&ker);

    /* ----- Combine ----- */
    if (out_n == 1) {
        Expr* out = eqs[0];
        free(eqs);
        return out;
    }
    Expr* and_e = expr_new_function(expr_new_symbol(SYM_And), eqs, out_n);
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
