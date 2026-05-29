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

#include "attr.h"
#include "eval.h"
#include "groebner.h"
#include "internal.h"
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

/* Booleans / sentinels that should never appear inside a polynomial.
 * Mathematical constants like Pi, E, EulerGamma are *not* filtered: the
 * inverse-function pre-pass routinely produces them (e.g. `ArcSin[1] ->
 * Pi/2`), and treating them as parameter symbols in the polynomial ring
 * is the natural way for them to flow through to the final equation. */
static bool is_system_constant_sym(const char* s) {
    return s == SYM_True || s == SYM_False || s == SYM_Infinity;
}

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

/* Walk `e` collecting every non-constant symbol that isn't in `elim`. */
static void collect_main_vars(const Expr* e, Expr** elim, size_t n_elim,
                              VarSet* out) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (!is_system_constant_sym(e->data.symbol)) {
            if (!var_in_list(e, elim, n_elim)) {
                varset_push(out, (Expr*)e);
            }
        }
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    /* The head itself may be a symbol like `f` in `f[x]`; if `f` is a
     * generic head with no elim binding, treat it as a parameter symbol
     * (Mathematica's `Eliminate[{f == x + y, x == 1}, x]` keeps `f` as a
     * parameter).  Builtin heads will already be elim-irrelevant so this
     * is safe. */
    if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* hsym = e->data.function.head->data.symbol;
        /* Skip "structural" heads we know are operators, not parameters. */
        if (hsym != SYM_Plus && hsym != SYM_Times && hsym != SYM_Power
         && hsym != SYM_List && hsym != SYM_And  && hsym != SYM_Equal
         && hsym != SYM_Sin  && hsym != SYM_Cos  && hsym != SYM_Tan
         && hsym != SYM_Sinh && hsym != SYM_Cosh && hsym != SYM_Tanh
         && hsym != SYM_Exp  && hsym != SYM_Log
         && hsym != SYM_ArcSin && hsym != SYM_ArcCos && hsym != SYM_ArcTan
         && hsym != SYM_ArcSinh && hsym != SYM_ArcCosh && hsym != SYM_ArcTanh
         && !is_system_constant_sym(hsym)) {
            /* Treat `f` in `f[args]` as a free parameter unless an arg
             * is itself elim-relevant. */
            if (!var_in_list(e->data.function.head, elim, n_elim)) {
                varset_push(out, e->data.function.head);
            }
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        collect_main_vars(e->data.function.args[i], elim, n_elim, out);
    }
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
/*  Equation -> polynomial                                             */
/* ------------------------------------------------------------------ */

/* Reduce `Equal[a, b]` to the polynomial `a - b` (expanded and
 * evaluated to fix-point).  A bare polynomial input passes through
 * untouched (after Expand).  Caller owns the return. */
static Expr* normalise_equation(const Expr* eq) {
    Expr* base;
    if (is_equal(eq) && eq->data.function.arg_count == 2) {
        Expr* a = expr_copy(eq->data.function.args[0]);
        Expr* b = expr_copy(eq->data.function.args[1]);
        base = internal_subtract((Expr*[]){ a, b }, 2);
    } else {
        base = expr_copy((Expr*)eq);
    }
    Expr* expanded = internal_expand((Expr*[]){ base }, 1);
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
                                 GB_ORDER_LEX, (int)n_elim);
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

    /* ----- Sentinels: empty -> True; contains constant -> False ----- */
    if (out_n == 0) {
        gb_basis_free(G, 0);
        free(all_vars);
        for (size_t i = 0; i < n_eq; i++) expr_free(polys[i]);
        free(polys);
        varset_free(&mains);
        if (wrap_vars) expr_free(wrap_vars);
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
