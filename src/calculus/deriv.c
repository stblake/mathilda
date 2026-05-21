/*
 * deriv.c -- Native C implementation of Mathematica-style differentiation.
 *
 * This module replaces the fragile rule-based bootstrap in
 * src/internal/deriv.m with a direct, dispatch-driven implementation.
 *
 * Overview
 * --------
 * The two key entry points are the builtins D (partial derivative) and
 * Dt (total derivative). Both ultimately funnel through a single
 * recursive core, ``compute_deriv``, parameterised by an optional
 * differentiation variable. When the variable is non-NULL we compute a
 * partial derivative treating everything else as constant (using a fast
 * FreeQ-style walk to short-circuit constant sub-trees). When the
 * variable is NULL we compute a total derivative -- unknown symbols
 * then participate as ``Dt[sym]`` terms.
 *
 * Why this is faster than the rule-based implementation
 * -----------------------------------------------------
 * The old deriv.m relied on ~60 DownValues. Each call to D[f, x] would:
 *   * scan the DownValues list for D linearly,
 *   * attempt pattern matching against every rule head (Plus, Times,
 *     Power, every elementary function, ...),
 *   * run ``/;`` side-conditions such as FreeQ,
 *   * perform attempt-evaluate/backtrack cycles in the matcher,
 *   * recursively re-evaluate the result through the full rule engine.
 *
 * In contrast, this module performs a single head-symbol strcmp dispatch
 * per call, constructs the derivative expression directly, and lets the
 * outer evaluator simplify arithmetic. Crucially, the constant-detection
 * step uses a tailored structural traversal (expr_free_of) that avoids
 * calling out to the generic FreeQ builtin.
 *
 * Returned expressions
 * --------------------
 * Every builder below produces plain un-reduced expression trees (e.g.
 * Plus[0, x] or Times[1, x]). The outer Mathilda evaluator runs a full
 * fixed-point reduction on the value we return, so Plus[0, ...],
 * Times[1, ...], and all subsequent chain-rule simplifications fold
 * automatically. This keeps the code readable and avoids duplicating
 * the arithmetic simplifier.
 *
 * Memory ownership
 * ----------------
 * Every helper that returns an ``Expr*`` returns a freshly allocated
 * tree owned by the caller. Input expressions are never mutated;
 * sub-expressions that need to be reused are always deep-copied.
 */

#include "deriv.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "print.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------- */
/* Tiny expression builders                                                */
/* ---------------------------------------------------------------------- */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }

/* Build a Function expression with the given head (already owned) and
 * a stack-provided argument array. Arguments are moved (ownership
 * transferred) into the new node. */
static Expr* mk_fn_take(Expr* head, Expr** args, size_t n) {
    /* expr_new_function memcpys the pointer slice, so a local stack
     * buffer is enough. */
    return expr_new_function(head, args, n);
}

static Expr* mk_fn1(const char* name, Expr* a) {
    Expr* args[1] = { a };
    return mk_fn_take(mk_sym(name), args, 1);
}

static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return mk_fn_take(mk_sym(name), args, 2);
}

/* Take ownership of `items` pointers and wrap them as name[items...].
 * The `items` array itself is freed here. */
static Expr* mk_fnN_adopt(const char* name, Expr** items, size_t n) {
    Expr* r = expr_new_function(mk_sym(name), items, n);
    free(items);
    return r;
}

/* Like mk_fnN_adopt but the head is an arbitrary pre-built expression
 * (e.g. Derivative[k][f] when producing Derivative[k][f][g]). The
 * `items` pointer must be a heap-allocated array owned by the caller;
 * this helper frees it. */
static Expr* mk_fn_head_adopt(Expr* head, Expr** items, size_t n) {
    Expr* r = expr_new_function(head, items, n);
    free(items);
    return r;
}

/* Wrap a single expression `arg` under an arbitrary pre-built head,
 * producing `head[arg]`. Used to build things like Derivative[1][f]
 * and Derivative[1][f][g]. */
static Expr* mk_fn_head1(Expr* head, Expr* arg) {
    Expr* items[1] = { arg };
    return expr_new_function(head, items, 1);
}

/* Negation helper: build Times[-1, a]. */
static Expr* mk_neg(Expr* a) {
    return mk_fn2("Times", mk_int(-1), a);
}

/* ---------------------------------------------------------------------- */
/* Predicates                                                              */
/* ---------------------------------------------------------------------- */

static bool is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol, name) == 0;
}

static bool is_lit_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* True if ``f`` contains no subexpression structurally equal to ``x``.
 *
 * This is a hot path: it is called once per recursive deriv step to
 * short-circuit constant sub-trees. We avoid expr_hash and the generic
 * FreeQ builtin to keep the inner loop allocation-free. */
static bool expr_free_of(const Expr* f, const Expr* x) {
    if (expr_eq(f, x)) return false;
    if (f->type == EXPR_FUNCTION) {
        if (!expr_free_of(f->data.function.head, x)) return false;
        for (size_t i = 0; i < f->data.function.arg_count; i++) {
            if (!expr_free_of(f->data.function.args[i], x)) return false;
        }
    }
    return true;
}

/* Recognise the symbols that Dt treats as absolute constants. These
 * mirror the rule set in the original deriv.m plus a few additional
 * Mathematica-compatible constants. */
static bool is_dt_constant_symbol(const Expr* e) {
    if (!e || e->type != EXPR_SYMBOL) return false;
    const char* s = e->data.symbol;
    return s == SYM_Pi || s == SYM_E || s == SYM_I ||
           s == SYM_Infinity || s == SYM_ComplexInfinity ||
           s == SYM_EulerGamma || s == SYM_Catalan ||
           s == SYM_GoldenRatio || s == SYM_Degree ||
           s == SYM_True || s == SYM_False;
}

/* ---------------------------------------------------------------------- */
/* Derivative of an elementary unary function                              */
/* ---------------------------------------------------------------------- */

/* Given the head name of a known elementary single-argument function F
 * and its argument g, produce F'(g). Returns NULL if F is not
 * recognised.
 *
 * The argument `g` is NOT consumed -- we copy it as needed.
 */
static Expr* elementary_fprime(const char* name, Expr* g) {
    /* --- trigonometric --- */
    if (!strcmp(name, "Sin")) return mk_fn1("Cos", expr_copy(g));
    if (!strcmp(name, "Cos")) return mk_neg(mk_fn1("Sin", expr_copy(g)));
    if (!strcmp(name, "Tan")) return mk_fn2("Power", mk_fn1("Sec", expr_copy(g)), mk_int(2));
    if (!strcmp(name, "Cot")) return mk_neg(mk_fn2("Power", mk_fn1("Csc", expr_copy(g)), mk_int(2)));
    if (!strcmp(name, "Sec")) return mk_fn2("Times", mk_fn1("Sec", expr_copy(g)), mk_fn1("Tan", expr_copy(g)));
    if (!strcmp(name, "Csc")) return mk_neg(mk_fn2("Times", mk_fn1("Csc", expr_copy(g)), mk_fn1("Cot", expr_copy(g))));

    /* --- inverse trigonometric --- */
    /* D[ArcSin[g], g] = 1 / Sqrt[1 - g^2]; we encode 1/Sqrt[u] as Power[u, -1/2]. */
    if (!strcmp(name, "ArcSin")) {
        Expr* one_minus_g2 = mk_fn2("Plus", mk_int(1),
                                    mk_neg(mk_fn2("Power", expr_copy(g), mk_int(2))));
        return mk_fn2("Power", mk_fn1("Sqrt", one_minus_g2), mk_int(-1));
    }
    if (!strcmp(name, "ArcCos")) {
        Expr* one_minus_g2 = mk_fn2("Plus", mk_int(1),
                                    mk_neg(mk_fn2("Power", expr_copy(g), mk_int(2))));
        return mk_neg(mk_fn2("Power", mk_fn1("Sqrt", one_minus_g2), mk_int(-1)));
    }
    if (!strcmp(name, "ArcTan")) {
        Expr* one_plus_g2 = mk_fn2("Plus", mk_int(1), mk_fn2("Power", expr_copy(g), mk_int(2)));
        return mk_fn2("Power", one_plus_g2, mk_int(-1));
    }
    if (!strcmp(name, "ArcCot")) {
        Expr* one_plus_g2 = mk_fn2("Plus", mk_int(1), mk_fn2("Power", expr_copy(g), mk_int(2)));
        return mk_neg(mk_fn2("Power", one_plus_g2, mk_int(-1)));
    }
    if (!strcmp(name, "ArcSec")) {
        /* 1 / (g^2 * Sqrt[1 - 1/g^2]) */
        Expr* g2 = mk_fn2("Power", expr_copy(g), mk_int(2));
        Expr* inv_g2 = mk_fn2("Power", expr_copy(g), mk_int(-2));
        Expr* radicand = mk_fn2("Plus", mk_int(1), mk_neg(inv_g2));
        Expr* denom = mk_fn2("Times", g2, mk_fn1("Sqrt", radicand));
        return mk_fn2("Power", denom, mk_int(-1));
    }
    if (!strcmp(name, "ArcCsc")) {
        Expr* g2 = mk_fn2("Power", expr_copy(g), mk_int(2));
        Expr* inv_g2 = mk_fn2("Power", expr_copy(g), mk_int(-2));
        Expr* radicand = mk_fn2("Plus", mk_int(1), mk_neg(inv_g2));
        Expr* denom = mk_fn2("Times", g2, mk_fn1("Sqrt", radicand));
        return mk_neg(mk_fn2("Power", denom, mk_int(-1)));
    }

    /* --- hyperbolic --- */
    if (!strcmp(name, "Sinh")) return mk_fn1("Cosh", expr_copy(g));
    if (!strcmp(name, "Cosh")) return mk_fn1("Sinh", expr_copy(g));
    if (!strcmp(name, "Tanh")) return mk_fn2("Power", mk_fn1("Sech", expr_copy(g)), mk_int(2));
    if (!strcmp(name, "Coth")) return mk_neg(mk_fn2("Power", mk_fn1("Csch", expr_copy(g)), mk_int(2)));
    if (!strcmp(name, "Sech")) return mk_neg(mk_fn2("Times", mk_fn1("Sech", expr_copy(g)), mk_fn1("Tanh", expr_copy(g))));
    if (!strcmp(name, "Csch")) return mk_neg(mk_fn2("Times", mk_fn1("Csch", expr_copy(g)), mk_fn1("Coth", expr_copy(g))));

    /* --- inverse hyperbolic --- */
    if (!strcmp(name, "ArcSinh")) {
        Expr* radicand = mk_fn2("Plus", mk_int(1), mk_fn2("Power", expr_copy(g), mk_int(2)));
        return mk_fn2("Power", mk_fn1("Sqrt", radicand), mk_int(-1));
    }
    if (!strcmp(name, "ArcCosh")) {
        Expr* radicand = mk_fn2("Plus", mk_fn2("Power", expr_copy(g), mk_int(2)), mk_int(-1));
        return mk_fn2("Power", mk_fn1("Sqrt", radicand), mk_int(-1));
    }
    if (!strcmp(name, "ArcTanh")) {
        Expr* denom = mk_fn2("Plus", mk_int(1), mk_neg(mk_fn2("Power", expr_copy(g), mk_int(2))));
        return mk_fn2("Power", denom, mk_int(-1));
    }
    if (!strcmp(name, "ArcCoth")) {
        Expr* denom = mk_fn2("Plus", mk_int(1), mk_neg(mk_fn2("Power", expr_copy(g), mk_int(2))));
        return mk_fn2("Power", denom, mk_int(-1));
    }
    if (!strcmp(name, "ArcSech")) {
        Expr* g2 = mk_fn2("Power", expr_copy(g), mk_int(2));
        Expr* inv_g2 = mk_fn2("Power", expr_copy(g), mk_int(-2));
        Expr* radicand = mk_fn2("Plus", inv_g2, mk_int(-1));
        Expr* denom = mk_fn2("Times", g2, mk_fn1("Sqrt", radicand));
        return mk_neg(mk_fn2("Power", denom, mk_int(-1)));
    }
    if (!strcmp(name, "ArcCsch")) {
        Expr* g2 = mk_fn2("Power", expr_copy(g), mk_int(2));
        Expr* inv_g2 = mk_fn2("Power", expr_copy(g), mk_int(-2));
        Expr* radicand = mk_fn2("Plus", mk_int(1), inv_g2);
        Expr* denom = mk_fn2("Times", g2, mk_fn1("Sqrt", radicand));
        return mk_neg(mk_fn2("Power", denom, mk_int(-1)));
    }

    /* --- exp/log and sqrt --- */
    if (!strcmp(name, "Exp")) return mk_fn1("Exp", expr_copy(g));
    if (!strcmp(name, "Log")) return mk_fn2("Power", expr_copy(g), mk_int(-1));
    if (!strcmp(name, "Sqrt")) {
        /* d/dg Sqrt[g] = 1/(2 Sqrt[g]) */
        return mk_fn2("Power", mk_fn2("Times", mk_int(2), mk_fn1("Sqrt", expr_copy(g))), mk_int(-1));
    }

    return NULL; /* not a recognised elementary unary */
}

/* ---------------------------------------------------------------------- */
/* Core recursive derivative                                               */
/* ---------------------------------------------------------------------- */

/* Forward declarations. */
static Expr* compute_deriv(Expr* f, Expr* x);

/* Shortcut: return D[g, x] as a fresh tree. When x is NULL, return
 * Dt[g] instead. For symbols and numeric atoms the answer is folded
 * immediately; for compound expressions we recurse through
 * compute_deriv so constants short-circuit there too. */
static Expr* deriv_of(Expr* g, Expr* x) {
    return compute_deriv(g, x);
}

/* The chain rule applied to an ``f[g1, g2, ..., gn]`` expression whose
 * head is not one of the explicitly-handled elementary heads. Produces
 *
 *   Sum_{k : D[gk, x] != 0} Derivative[0..1_k..0][f][g1..gn] * D[gk, x]
 *
 * If all partials vanish the result is 0. If exactly one contributes,
 * the outer Plus is elided. */
static Expr* chain_rule_unknown(Expr* f, Expr* x) {
    Expr* head = f->data.function.head;
    size_t n = f->data.function.arg_count;
    Expr** args = f->data.function.args;

    /* Collect (nonzero-derivative, partial-index) terms. */
    Expr** terms = malloc(sizeof(Expr*) * n);
    size_t nterms = 0;
    for (size_t k = 0; k < n; k++) {
        Expr* dk = deriv_of(args[k], x);
        if (is_lit_zero(dk)) { expr_free(dk); continue; }

        /* Build the indicator Derivative[0,..,1,..,0]. */
        Expr** idx = malloc(sizeof(Expr*) * n);
        for (size_t j = 0; j < n; j++) idx[j] = mk_int(j == k ? 1 : 0);
        Expr* op = mk_fnN_adopt("Derivative", idx, n);          /* Derivative[...]   */
        Expr* op_f = mk_fn_head1(op, expr_copy(head));          /* Derivative[...][f] */

        /* Apply to the original argument list. */
        Expr** gcopy = malloc(sizeof(Expr*) * n);
        for (size_t j = 0; j < n; j++) gcopy[j] = expr_copy(args[j]);
        Expr* applied = mk_fn_head_adopt(op_f, gcopy, n);

        terms[nterms++] = mk_fn2("Times", applied, dk);
    }

    if (nterms == 0) { free(terms); return mk_int(0); }
    if (nterms == 1) { Expr* r = terms[0]; free(terms); return r; }
    return mk_fnN_adopt("Plus", terms, nterms);
}

/* Differentiate an expression of the form
 *
 *     Derivative[n1, ..., nm][f] [ g1, ..., gn ]
 *
 * Each nonzero partial derivative advances the corresponding index in
 * the Derivative head. Returns NULL if the expression doesn't match
 * this shape (so the caller can try other dispatch rules). */
static Expr* deriv_of_derivative_form(Expr* f, Expr* x) {
    Expr* head = f->data.function.head;
    if (head->type != EXPR_FUNCTION) return NULL;
    if (head->data.function.arg_count != 1) return NULL;   /* need exactly one f */

    Expr* outer = head->data.function.head;                /* Derivative[n1..nm] */
    if (outer->type != EXPR_FUNCTION) return NULL;
    if (!is_sym(outer->data.function.head, "Derivative")) return NULL;

    size_t m = outer->data.function.arg_count;
    size_t n = f->data.function.arg_count;
    if (m != n) return NULL;                               /* arity mismatch */

    Expr* f0 = head->data.function.args[0];                /* the underlying f */
    Expr** args = f->data.function.args;

    Expr** terms = malloc(sizeof(Expr*) * n);
    size_t nterms = 0;
    for (size_t k = 0; k < n; k++) {
        Expr* dk = deriv_of(args[k], x);
        if (is_lit_zero(dk)) { expr_free(dk); continue; }

        /* Advance the k-th order index. If it's a plain integer, bump it;
         * otherwise wrap it in a Plus so the evaluator can simplify. */
        Expr** new_idx = malloc(sizeof(Expr*) * m);
        for (size_t j = 0; j < m; j++) {
            Expr* orig = outer->data.function.args[j];
            if (j == k) {
                if (orig->type == EXPR_INTEGER) {
                    new_idx[j] = mk_int(orig->data.integer + 1);
                } else {
                    new_idx[j] = mk_fn2("Plus", expr_copy(orig), mk_int(1));
                }
            } else {
                new_idx[j] = expr_copy(orig);
            }
        }
        Expr* new_op = mk_fnN_adopt("Derivative", new_idx, m);
        Expr* new_op_f = mk_fn_head1(new_op, expr_copy(f0));

        Expr** gcopy = malloc(sizeof(Expr*) * n);
        for (size_t j = 0; j < n; j++) gcopy[j] = expr_copy(args[j]);
        Expr* applied = mk_fn_head_adopt(new_op_f, gcopy, n);

        terms[nterms++] = mk_fn2("Times", applied, dk);
    }

    if (nterms == 0) { free(terms); return mk_int(0); }
    if (nterms == 1) { Expr* r = terms[0]; free(terms); return r; }
    return mk_fnN_adopt("Plus", terms, nterms);
}

/* Core dispatch. If `x` is non-NULL we are computing D[f, x]; otherwise
 * we are computing Dt[f]. The two paths share 95% of their logic --
 * only constant-handling and the "unknown symbol" base-case differ. */
static Expr* compute_deriv(Expr* f, Expr* x) {
    /* --- Partial-derivative base cases. --- */
    if (x) {
        if (expr_free_of(f, x)) return mk_int(0);
        if (expr_eq(f, x)) return mk_int(1);
    } else {
        /* Dt mode: numeric atoms and distinguished constants vanish. */
        if (f->type == EXPR_INTEGER || f->type == EXPR_REAL || f->type == EXPR_BIGINT) {
            return mk_int(0);
        }
        if (is_dt_constant_symbol(f)) return mk_int(0);
        if (f->type == EXPR_SYMBOL) {
            /* Unknown symbols are treated as implicit functions of an
             * implicit parameter -- their total derivative is Dt[x]. */
            return mk_fn1("Dt", expr_copy(f));
        }
    }

    /* At this point f must be compound. */
    if (f->type != EXPR_FUNCTION) {
        /* Strings and other atoms we cannot differentiate -- bail. */
        return NULL;
    }

    Expr* head = f->data.function.head;
    size_t n = f->data.function.arg_count;
    Expr** args = f->data.function.args;

    /* ------------------------------------------------------------------ */
    /* Head is a plain symbol -- dispatch on the head name.               */
    /* ------------------------------------------------------------------ */
    if (head->type == EXPR_SYMBOL) {
        const char* h = head->data.symbol;

        /* --- Plus: sum of derivatives. --- */
        if (h == SYM_Plus) {
            Expr** ts = malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) ts[i] = deriv_of(args[i], x);
            return mk_fnN_adopt("Plus", ts, n);
        }

        /* --- Times: general product rule. For n factors this is
         *      sum_i (d/dx arg_i) * prod_{j!=i} arg_j.
         * We build each term directly; the evaluator handles
         * simplification of Times[0, ...] and trivial identities. --- */
        if (h == SYM_Times) {
            Expr** ts = malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) {
                Expr** factors = malloc(sizeof(Expr*) * n);
                for (size_t j = 0; j < n; j++) {
                    factors[j] = (i == j) ? deriv_of(args[j], x)
                                          : expr_copy(args[j]);
                }
                ts[i] = mk_fnN_adopt("Times", factors, n);
            }
            return mk_fnN_adopt("Plus", ts, n);
        }

        /* --- Power: d/dx f^g = f^(g-1) * (g*f' + f*Log[f]*g'). ---
         *
         * When g turns out to be independent of x we collapse to the
         * usual g * f^(g-1) * f'. This avoids introducing spurious
         * Log[f] terms that the outer evaluator would then have to
         * fight to remove (and which are problematic for non-positive
         * bases). We detect "g is constant" by computing d g / d x once
         * and checking for the literal 0. */
        if (h == SYM_Power && n == 2) {
            Expr* a = args[0];
            Expr* b = args[1];
            Expr* da = deriv_of(a, x);
            Expr* db = deriv_of(b, x);

            if (is_lit_zero(db)) {
                expr_free(db);
                Expr* bm1 = mk_fn2("Plus", expr_copy(b), mk_int(-1));
                Expr* f_pow = mk_fn2("Power", expr_copy(a), bm1);
                Expr** factors = malloc(sizeof(Expr*) * 3);
                factors[0] = expr_copy(b);
                factors[1] = f_pow;
                factors[2] = da;
                return mk_fnN_adopt("Times", factors, 3);
            }

            Expr* bm1 = mk_fn2("Plus", expr_copy(b), mk_int(-1));
            Expr* f_pow = mk_fn2("Power", expr_copy(a), bm1);
            Expr* t1 = mk_fn2("Times", expr_copy(b), da);
            Expr* t2_inner = mk_fn2("Times", expr_copy(a), mk_fn1("Log", expr_copy(a)));
            Expr* t2 = mk_fn2("Times", t2_inner, db);
            Expr* inner = mk_fn2("Plus", t1, t2);
            return mk_fn2("Times", f_pow, inner);
        }

        /* --- List: thread element-wise. --- */
        if (h == SYM_List) {
            Expr** ts = malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) ts[i] = deriv_of(args[i], x);
            return mk_fnN_adopt("List", ts, n);
        }

        /* --- RootSum[Function[t, p], Function[t, body]]: thread the
         *     derivative through the body Function.  The bound variable
         *     `t` ranges over roots of p[t] == 0, which are independent
         *     of the outer differentiation variable x; only the body
         *     can depend on x.  See src/root.c for the construct's
         *     definition. */
        if (h == SYM_RootSum && n == 2) {
            Expr* fn1 = args[0];
            Expr* fn2 = args[1];
            if (fn2->type == EXPR_FUNCTION
                && fn2->data.function.head->type == EXPR_SYMBOL
                && fn2->data.function.head->data.symbol == SYM_Function) {
                /* Two accepted body forms (the bound variable is
                 * always independent of x, so D threads through):
                 *   2-arg  Function[t, body]      — named bound var
                 *   1-arg  Function[body]         — Slot[_]-bound
                 *                                   ((... &) syntax)
                 * Both are produced by root_make_rootsum at different
                 * historical points; we accept either. */
                size_t fa = fn2->data.function.arg_count;
                if (fa == 2) {
                    Expr* bvar = fn2->data.function.args[0];
                    Expr* body = fn2->data.function.args[1];
                    Expr* dbody = deriv_of(body, x);
                    Expr* dbody_eval = eval_and_free(dbody);
                    Expr* new_fn2 = mk_fn2("Function",
                        expr_copy(bvar), dbody_eval);
                    return mk_fn2("RootSum", expr_copy(fn1), new_fn2);
                }
                if (fa == 1) {
                    Expr* body = fn2->data.function.args[0];
                    Expr* dbody = deriv_of(body, x);
                    Expr* dbody_eval = eval_and_free(dbody);
                    Expr* new_fn2 = mk_fn1("Function", dbody_eval);
                    return mk_fn2("RootSum", expr_copy(fn1), new_fn2);
                }
            }
            /* Unrecognised body form: leave RootSum unevaluated. */
            return NULL;
        }

        /* --- Log[b, f]: reduce to Log[f]/Log[b]. --- */
        if (h == SYM_Log && n == 2) {
            Expr* lb = mk_fn1("Log", expr_copy(args[0]));
            Expr* lf = mk_fn1("Log", expr_copy(args[1]));
            Expr* quot = mk_fn2("Times", lf, mk_fn2("Power", lb, mk_int(-1)));
            Expr* r = compute_deriv(quot, x);
            expr_free(quot);
            return r;
        }

        /* --- Known elementary unary function: F'(g) * D[g, x]. --- */
        if (n == 1) {
            Expr* fp = elementary_fprime(h, args[0]);
            if (fp) {
                Expr* dg = deriv_of(args[0], x);
                return mk_fn2("Times", fp, dg);
            }
            /* Unknown single-argument function -- standard chain rule:
             *     Derivative[1][f][g] * D[g, x]. */
            Expr* op = mk_fn1("Derivative", mk_int(1));         /* Derivative[1]      */
            Expr* op_f = mk_fn_head1(op, expr_copy(head));      /* Derivative[1][f]   */
            Expr* applied = mk_fn_head1(op_f, expr_copy(args[0])); /* Derivative[1][f][g] */
            Expr* dg = deriv_of(args[0], x);
            return mk_fn2("Times", applied, dg);
        }

        /* --- Unknown multi-argument function: full chain rule. --- */
        return chain_rule_unknown(f, x);
    }

    /* ------------------------------------------------------------------ */
    /* Head is itself a function. The only pattern we can reduce is       */
    /* Derivative[...][f][args...].                                       */
    /* ------------------------------------------------------------------ */
    if (head->type == EXPR_FUNCTION) {
        Expr* r = deriv_of_derivative_form(f, x);
        if (r) return r;
    }

    /* Give up -- caller keeps the expression unevaluated. */
    return NULL;
}

/* ---------------------------------------------------------------------- */
/* Higher-order: D[f, {x, n}]                                              */
/* ---------------------------------------------------------------------- */

/* Returns a new expression representing the n-th partial derivative of
 * f with respect to x. For n == 0 the input is deep-copied. */
static Expr* higher_order_partial(Expr* f, Expr* x, int64_t order) {
    if (order <= 0) return expr_copy(f);
    Expr* current = expr_copy(f);
    for (int64_t i = 0; i < order; i++) {
        Expr* nxt = compute_deriv(current, x);
        if (!nxt) {
            /* Cannot differentiate further -- leave as D[current, x]
             * and wrap any remaining orders as nested D calls. In
             * practice this only occurs for truly opaque heads. */
            Expr* wrap = mk_fn2("D", current, expr_copy(x));
            for (int64_t k = i + 1; k < order; k++) {
                wrap = mk_fn2("D", wrap, expr_copy(x));
            }
            return wrap;
        }
        /* The outer evaluator will also reduce the final expression,
         * but we need a simplified form between iterations so that
         * expr_free_of short-circuits work on the next pass.  Note that
         * evaluate() returns a new tree without freeing its argument;
         * we must explicitly free `nxt`. */
        Expr* reduced = evaluate(nxt);
        expr_free(nxt);
        expr_free(current);
        current = reduced;
    }
    return current;
}

/* ---------------------------------------------------------------------- */
/* Builtin: D                                                              */
/* ---------------------------------------------------------------------- */

/* Parse a second-argument spec which may be either
 *     x            -- a bare variable, order 1
 *     {x, n}       -- an integer-order spec
 *     {x, k}       -- a symbolic-order spec (k an Expr, not an Integer)
 * On success, *var is set to a NON-owned pointer into the spec.
 *   - For an integer order, *order is set to the integer and *order_expr
 *     is NULL.
 *   - For a symbolic order, *order is set to the sentinel -1 and
 *     *order_expr is a NON-owned pointer to the spec's symbolic order
 *     argument.
 * Returns false if the spec cannot be interpreted. */
static bool parse_var_spec(Expr* spec, Expr** var, int64_t* order,
                           Expr** order_expr) {
    if (order_expr) *order_expr = NULL;
    if (spec->type == EXPR_FUNCTION &&
        is_sym(spec->data.function.head, "List")) {
        /* A List spec must have the {var, n} form. Any other shape
         * (e.g., {x}, {x, y, z}) is a malformed multiple-derivative
         * specifier; reject so the caller can emit D::dvar and leave
         * the call unevaluated (matching Mathematica). */
        if (spec->data.function.arg_count != 2) {
            return false;
        }
        Expr* k = spec->data.function.args[1];
        if (k->type == EXPR_INTEGER) {
            *var = spec->data.function.args[0];
            *order = k->data.integer;
            return *order >= 0;
        }
        /* Reject non-integer numeric orders (Real, Rational, Complex);
         * MMA's D::dvar requires n to be symbolic or non-negative integer. */
        if (k->type == EXPR_REAL || k->type == EXPR_BIGINT ||
            (k->type == EXPR_FUNCTION &&
             k->data.function.head->type == EXPR_SYMBOL &&
             (k->data.function.head->data.symbol == SYM_Rational ||
              k->data.function.head->data.symbol == SYM_Complex))) {
            return false;
        }
        /* Symbolic order. */
        *var = spec->data.function.args[0];
        *order = -1;
        if (order_expr) *order_expr = k;
        return true;
    }
    *var = spec;
    *order = 1;
    return true;
}

/* Emit MMA-style D::dvar message for malformed multiple-derivative
 * specifiers. Only called when parse_var_spec rejects a List spec. */
static void emit_dvar_message(Expr* spec) {
    if (spec->type != EXPR_FUNCTION ||
        !is_sym(spec->data.function.head, "List")) {
        return;
    }
    char* s = expr_to_string(spec);
    printf("D::dvar: Multiple derivative specifier %s does not have "
           "the form {variable, n}, where n is symbolic or a "
           "non-negative integer.\n", s ? s : "?");
    free(s);
}

/* Compute D[f, {var, k}] for a symbolic order k. The algorithm pattern-
 * matches the small set of forms with closed-form symbolic-order
 * derivatives, recurses through additive structure, pulls var-free
 * factors out of products, and otherwise returns NULL (caller falls
 * back to leaving D[f, {var, k}] unevaluated).
 *
 *   D[c, {x, k}]              = 0                          (c free of x)
 *   D[c f, {x, k}]            = c D[f, {x, k}]              (c free of x)
 *   D[a + b + ..., {x, k}]    = D[a, {x, k}] + ...
 *   D[Power[x, n], {x, k}]    = FactorialPower[n, k] x^(n-k)  (n free of x)
 *   D[Power[b, x], {x, k}]    = b^x Log[b]^k                  (b free of x)
 *
 * Sound for non-negative integer k (MMA's Piecewise conventions).
 * Returns a fresh Expr* on success, NULL when no symbolic-k closed form
 * applies. */
static Expr* compute_deriv_symbolic_order(Expr* f, Expr* var, Expr* k);

static Expr* compute_deriv_symbolic_order(Expr* f, Expr* var, Expr* k) {
    if (!f || !var || !k) return NULL;

    /* Var-free constant: D[c, {x, k}] = 0 for k >= 1. MMA also gives 0;
     * we match. (Strictly Piecewise[{{c, k==0}}, 0], but our convention
     * matches the unsubscripted MMA result.) */
    if (expr_free_of(f, var)) return mk_int(0);

    /* The variable itself: D[x, {x, k}] = If[k == 1, 1, 0] but MMA
     * returns Piecewise; for the symbolic-k path we leave unevaluated
     * (no general single-Expr closed form). */
    if (expr_eq(f, var)) return NULL;

    /* Plus: distribute. */
    if (f->type == EXPR_FUNCTION && is_sym(f->data.function.head, "Plus")) {
        size_t n = f->data.function.arg_count;
        Expr** parts = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            Expr* d = compute_deriv_symbolic_order(f->data.function.args[i], var, k);
            if (!d) {
                for (size_t j = 0; j < i; j++) expr_free(parts[j]);
                free(parts);
                return NULL;
            }
            parts[i] = d;
        }
        Expr* sum = expr_new_function(expr_new_symbol("Plus"), parts, n);
        free(parts);
        return eval_and_free(sum);
    }

    /* Times: pull out var-free factors. If all but one factor is var-
     * free, the var-bearing factor's symbolic-k derivative scales by
     * the product of the var-free factors. */
    if (f->type == EXPR_FUNCTION && is_sym(f->data.function.head, "Times")) {
        size_t n = f->data.function.arg_count;
        Expr** free_factors = malloc(sizeof(Expr*) * n);
        size_t free_count = 0;
        Expr* var_factor = NULL;
        size_t var_factor_count = 0;
        for (size_t i = 0; i < n; i++) {
            if (expr_free_of(f->data.function.args[i], var)) {
                free_factors[free_count++] = expr_copy(f->data.function.args[i]);
            } else {
                var_factor = f->data.function.args[i];
                var_factor_count++;
            }
        }
        if (var_factor_count == 0 || var_factor_count > 1) {
            /* All var-free already handled above; >1 var-bearing factors
             * have no general symbolic-k form (Leibniz produces a
             * binomial sum that the user can request via Sum). */
            for (size_t i = 0; i < free_count; i++) expr_free(free_factors[i]);
            free(free_factors);
            return NULL;
        }
        Expr* d_var = compute_deriv_symbolic_order(var_factor, var, k);
        if (!d_var) {
            for (size_t i = 0; i < free_count; i++) expr_free(free_factors[i]);
            free(free_factors);
            return NULL;
        }
        if (free_count == 0) { free(free_factors); return d_var; }
        Expr** product_args = malloc(sizeof(Expr*) * (free_count + 1));
        for (size_t i = 0; i < free_count; i++) product_args[i] = free_factors[i];
        product_args[free_count] = d_var;
        Expr* prod = expr_new_function(expr_new_symbol("Times"),
                                       product_args, free_count + 1);
        free(product_args);
        free(free_factors);
        return eval_and_free(prod);
    }

    /* Power forms with symbolic-k closed-form derivatives. */
    if (f->type == EXPR_FUNCTION && is_sym(f->data.function.head, "Power") &&
        f->data.function.arg_count == 2) {
        Expr* base = f->data.function.args[0];
        Expr* exp  = f->data.function.args[1];

        /* Power[var, n] with n free of var:
         *   D[x^n, {x, k}] = FactorialPower[n, k] * Power[x, n - k]. */
        if (expr_eq(base, var) && expr_free_of(exp, var)) {
            Expr* fp = mk_fn2("FactorialPower", expr_copy(exp), expr_copy(k));
            Expr* new_exp = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_copy(exp), mk_neg(expr_copy(k)) }, 2);
            Expr* new_pow = mk_fn2("Power", expr_copy(var), new_exp);
            Expr* product = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ fp, new_pow }, 2);
            return eval_and_free(product);
        }

        /* Power[b, var] with b free of var:
         *   D[b^x, {x, k}] = b^x * Log[b]^k. */
        if (expr_eq(exp, var) && expr_free_of(base, var)) {
            Expr* same  = mk_fn2("Power", expr_copy(base), expr_copy(var));
            Expr* logb  = mk_fn1("Log", expr_copy(base));
            Expr* logbk = mk_fn2("Power", logb, expr_copy(k));
            Expr* product = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ same, logbk }, 2);
            return eval_and_free(product);
        }
    }

    /* Sin/Cos/Sinh/Cosh of a linear-in-var argument.
     *
     * For u(var) = a*var + b with a, b free of var, the chain rule gives
     *   D^k/dvar^k F[u] = a^k * F^(k)[u]
     * with closed-form k-th derivatives matching Mathematica:
     *   Sin^(k)[t]  = Sin[k Pi/2 + t]
     *   Cos^(k)[t]  = Cos[k Pi/2 + t]
     *   Cosh^(k)[t] = (-I)^k Cos[k Pi/2 - I t]
     *   Sinh^(k)[t] = I (-I)^k Sin[k Pi/2 - I t]
     *
     * Linearity is detected by computing du/dvar; if that derivative is
     * free of var then u is linear in var and the residual is exactly a.
     * Non-linear arguments fall through (caller leaves D unevaluated). */
    if (f->type == EXPR_FUNCTION && f->data.function.arg_count == 1 &&
        f->data.function.head->type == EXPR_SYMBOL) {
        const char* h = f->data.function.head->data.symbol;
        bool is_trig  = (h == SYM_Sin) || (h == SYM_Cos);
        bool is_hyper = (h == SYM_Sinh) || (h == SYM_Cosh);
        if (is_trig || is_hyper) {
            Expr* u = f->data.function.args[0];
            Expr* du_raw = compute_deriv(u, var);
            if (du_raw) {
                Expr* a = evaluate(du_raw);
                expr_free(du_raw);
                if (expr_free_of(a, var)) {
                    /* half_pi_k = k * Pi / 2 */
                    Expr* half_pi_k = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_copy(k),
                                   expr_new_symbol("Pi"),
                                   mk_fn2("Power", mk_int(2), mk_int(-1)) }, 3);
                    Expr* ak = mk_fn2("Power", a, expr_copy(k));
                    Expr* result;
                    if (is_trig) {
                        /* a^k * F[k*Pi/2 + u] */
                        Expr* shifted = mk_fn2("Plus", half_pi_k, expr_copy(u));
                        Expr* trig    = mk_fn1(h, shifted);
                        result = mk_fn2("Times", ak, trig);
                    } else {
                        /* a^k * (-I)^k * Cos[k*Pi/2 - I u]      (Cosh)
                         * a^k * I * (-I)^k * Sin[k*Pi/2 - I u]  (Sinh) */
                        Expr* minus_I_u = mk_fn2("Times",
                            mk_fn2("Times", mk_int(-1), expr_new_symbol("I")),
                            expr_copy(u));
                        Expr* arg     = mk_fn2("Plus", half_pi_k, minus_I_u);
                        const char* trig_name = (h == SYM_Cosh) ? "Cos" : "Sin";
                        Expr* trig    = mk_fn1(trig_name, arg);
                        Expr* neg_I   = mk_fn2("Times", mk_int(-1),
                                               expr_new_symbol("I"));
                        Expr* neg_I_k = mk_fn2("Power", neg_I, expr_copy(k));
                        if (h == SYM_Cosh) {
                            result = expr_new_function(expr_new_symbol("Times"),
                                (Expr*[]){ ak, neg_I_k, trig }, 3);
                        } else {
                            result = expr_new_function(expr_new_symbol("Times"),
                                (Expr*[]){ ak, expr_new_symbol("I"),
                                           neg_I_k, trig }, 4);
                        }
                    }
                    return eval_and_free(result);
                }
                expr_free(a);
            }
        }
    }

    /* No symbolic-k closed form recognised for this shape. */
    return NULL;
}

/* Build an unevaluated D[f, {var, k}] tree. Used as the fallback when
 * compute_deriv_symbolic_order can't reduce the expression. */
static Expr* build_unevaluated_d(Expr* f, Expr* var, Expr* k) {
    Expr* spec = expr_new_function(expr_new_symbol("List"),
        (Expr*[]){ expr_copy(var), expr_copy(k) }, 2);
    return expr_new_function(expr_new_symbol("D"),
        (Expr*[]){ expr_copy(f), spec }, 2);
}

Expr* builtin_d(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;                 /* D[f] stays unevaluated */

    Expr* f = res->data.function.args[0];
    Expr** specs = &res->data.function.args[1];
    size_t nspecs = argc - 1;

    /* Sequentially apply each spec (handling mixed partials D[f, x, y, ...]). */
    Expr* current = expr_copy(f);
    for (size_t i = 0; i < nspecs; i++) {
        Expr* var = NULL;
        int64_t order = 1;
        Expr* order_expr = NULL;
        if (!parse_var_spec(specs[i], &var, &order, &order_expr)) {
            emit_dvar_message(specs[i]);       /* MMA-style D::dvar */
            expr_free(current);
            return NULL;                       /* malformed spec */
        }
        Expr* stepped = NULL;
        if (order >= 0) {
            stepped = higher_order_partial(current, var, order);
        } else {
            /* Symbolic-order spec. */
            Expr* sym = compute_deriv_symbolic_order(current, var, order_expr);
            if (sym) stepped = sym;
            else stepped = build_unevaluated_d(current, var, order_expr);
        }
        expr_free(current);
        current = stepped;
    }

    return current;
}

/* ---------------------------------------------------------------------- */
/* Builtin: Dt                                                             */
/* ---------------------------------------------------------------------- */

Expr* builtin_dt(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return NULL;

    Expr* f = res->data.function.args[0];

    /* Dt[f]: total derivative. */
    if (argc == 1) {
        return compute_deriv(f, NULL);        /* may be NULL to stay unevaluated */
    }

    /* Dt[f, var_specs...] is identical to D[f, var_specs...]: it gives
     * the partial derivative. Forward to the D path. */
    Expr* current = expr_copy(f);
    for (size_t i = 1; i < argc; i++) {
        Expr* var = NULL;
        int64_t order = 1;
        Expr* order_expr = NULL;
        if (!parse_var_spec(res->data.function.args[i], &var, &order, &order_expr)) {
            emit_dvar_message(res->data.function.args[i]);
            expr_free(current);
            return NULL;
        }
        Expr* stepped = NULL;
        if (order >= 0) {
            stepped = higher_order_partial(current, var, order);
        } else {
            Expr* sym = compute_deriv_symbolic_order(current, var, order_expr);
            stepped = sym ? sym : build_unevaluated_d(current, var, order_expr);
        }
        expr_free(current);
        current = stepped;
    }

    return current;
}

/* ---------------------------------------------------------------------- */
/* Reduce Derivative[n1,...,nm][Function[...]] to a plain Function[body']  */
/* ---------------------------------------------------------------------- */

/* Build a fresh Slot[i] expression. */
static Expr* make_slot(int64_t i) {
    Expr* idx[1] = { mk_int(i) };
    return expr_new_function(mk_sym("Slot"), idx, 1);
}

/* See deriv.h for the contract. */
Expr* derivative_of_pure_function(Expr* deriv_head, Expr* pure_fn) {
    if (!deriv_head || !pure_fn) return NULL;
    if (deriv_head->type != EXPR_FUNCTION) return NULL;
    if (!is_sym(deriv_head->data.function.head, "Derivative")) return NULL;
    size_t m = deriv_head->data.function.arg_count;
    if (m == 0) return NULL;

    /* All derivative orders must be nonnegative integers. */
    int64_t* orders = malloc(sizeof(int64_t) * m);
    for (size_t i = 0; i < m; i++) {
        Expr* oi = deriv_head->data.function.args[i];
        if (oi->type != EXPR_INTEGER || oi->data.integer < 0) {
            free(orders);
            return NULL;
        }
        orders[i] = oi->data.integer;
    }

    if (pure_fn->type != EXPR_FUNCTION ||
        !is_sym(pure_fn->data.function.head, "Function")) {
        free(orders);
        return NULL;
    }

    size_t fargc = pure_fn->data.function.arg_count;
    if (fargc == 0) { free(orders); return NULL; }

    /* Determine body and the variables to differentiate with respect to.
     * For each supported Function signature we identify m variables
     * (matching the arity of Derivative) and remember whether the rebuilt
     * Function preserves the original parameter/attribute shape. */
    Expr* body = NULL;
    Expr** vars = malloc(sizeof(Expr*) * m);
    bool slot_form = false;

    if (fargc == 1) {
        /* Function[body]: slot form. Differentiate wrt Slot[1..m]. */
        body = pure_fn->data.function.args[0];
        for (size_t i = 0; i < m; i++) vars[i] = make_slot((int64_t)(i + 1));
        slot_form = true;
    } else {
        /* Function[params, body, ...]. */
        Expr* params = pure_fn->data.function.args[0];
        body = pure_fn->data.function.args[1];

        if (params->type == EXPR_SYMBOL &&
            params->data.symbol == SYM_Null) {
            /* Function[Null, body, attrs]: slot form with attributes. */
            for (size_t i = 0; i < m; i++) vars[i] = make_slot((int64_t)(i + 1));
        } else if (params->type == EXPR_SYMBOL) {
            /* Function[x, body]: single named parameter. Only Derivative[n]
             * (m == 1) is meaningful here. */
            if (m != 1) { free(vars); free(orders); return NULL; }
            vars[0] = expr_copy(params);
        } else if (params->type == EXPR_FUNCTION &&
                   is_sym(params->data.function.head, "List")) {
            /* Function[{x1,...,xk}, body]: arity of Derivative must match. */
            size_t k = params->data.function.arg_count;
            if (k != m) { free(vars); free(orders); return NULL; }
            for (size_t i = 0; i < m; i++) {
                vars[i] = expr_copy(params->data.function.args[i]);
            }
        } else {
            free(vars);
            free(orders);
            return NULL;
        }
    }

    /* Apply each partial derivative sequentially, simplifying between
     * iterations so subsequent passes can see a reduced expression. */
    Expr* current = expr_copy(body);
    for (size_t i = 0; i < m; i++) {
        for (int64_t k = 0; k < orders[i]; k++) {
            Expr* d = compute_deriv(current, vars[i]);
            expr_free(current);
            if (!d) {
                for (size_t j = 0; j < m; j++) expr_free(vars[j]);
                free(vars);
                free(orders);
                return NULL;
            }
            current = evaluate(d);
            expr_free(d);
        }
    }

    for (size_t j = 0; j < m; j++) expr_free(vars[j]);
    free(vars);
    free(orders);

    /* Rebuild a Function with the original parameter/attribute shape but
     * the differentiated body. */
    if (slot_form) {
        Expr* fn_args[1] = { current };
        return expr_new_function(mk_sym("Function"), fn_args, 1);
    }

    size_t new_argc = fargc;
    Expr** new_args = malloc(sizeof(Expr*) * new_argc);
    new_args[0] = expr_copy(pure_fn->data.function.args[0]);
    new_args[1] = current;
    for (size_t i = 2; i < new_argc; i++) {
        new_args[i] = expr_copy(pure_fn->data.function.args[i]);
    }
    Expr* out = expr_new_function(mk_sym("Function"), new_args, new_argc);
    free(new_args);
    return out;
}

/* ---------------------------------------------------------------------- */
/* Reduce Derivative[n1,...,nm][f] for a symbol f with DownValues          */
/* ---------------------------------------------------------------------- */

/* Counter used to mint unique temporary symbol names when synthesising
 * the Function[{t1,...,tm}, f[t1,...,tm]] wrapper. The temporaries never
 * have OwnValues installed and are not registered with the symbol table,
 * so they cannot collide with anything in the body produced by the
 * DownValue rewrite (which would have to literally contain the same
 * unique name to clash, and the names are generated fresh per call). */
static int64_t deriv_temp_counter = 0;

/* See deriv.h for the contract. */
Expr* derivative_of_symbol(Expr* deriv_head, Expr* fsym) {
    if (!deriv_head || !fsym) return NULL;
    if (deriv_head->type != EXPR_FUNCTION) return NULL;
    if (!is_sym(deriv_head->data.function.head, "Derivative")) return NULL;
    if (fsym->type != EXPR_SYMBOL) return NULL;

    size_t m = deriv_head->data.function.arg_count;
    if (m == 0) return NULL;

    /* All derivative orders must be nonnegative integers; otherwise the
     * downstream derivative_of_pure_function would reject them and we
     * would have done substantial work for nothing. */
    for (size_t i = 0; i < m; i++) {
        Expr* oi = deriv_head->data.function.args[i];
        if (oi->type != EXPR_INTEGER || oi->data.integer < 0) return NULL;
    }

    /* Only proceed when the symbol carries at least one DownValue. For a
     * bare symbol with no rules, attempting to substitute would just give
     * us back f[t1,...,tm] and we would correctly fall through; the
     * early exit here is purely a fast path that avoids the call/eval
     * roundtrip in the common no-definition case. */
    if (!symtab_get_down_values(fsym->data.symbol)) return NULL;

    /* Mint fresh temporary variable symbols. We use a static counter so
     * names are globally unique per process invocation; this avoids the
     * need to register temporaries (and clean them up) in the symbol
     * table just to guarantee freshness. */
    int64_t my_id = ++deriv_temp_counter;

    Expr** vars = malloc(sizeof(Expr*) * m);
    for (size_t i = 0; i < m; i++) {
        char namebuf[64];
        snprintf(namebuf, sizeof(namebuf),
                 "Derivative$%lld$%zu", (long long)my_id, i + 1);
        vars[i] = mk_sym(namebuf);
    }

    /* Build f[t1,...,tm] and evaluate it. The DownValue lookup happens
     * during evaluation, so we get back the substituted body when a
     * single-symbol-argument rule like f[x_] := body matches. */
    Expr** call_args = malloc(sizeof(Expr*) * m);
    for (size_t i = 0; i < m; i++) call_args[i] = expr_copy(vars[i]);
    Expr* call = expr_new_function(expr_copy(fsym), call_args, m);
    free(call_args);

    Expr* body = evaluate(call);
    expr_free(call);

    /* If the call did not rewrite (no matching DownValue), abort. We
     * detect this by checking that the result is still structurally
     * f[t1,...,tm]. */
    bool unchanged = (body->type == EXPR_FUNCTION
                      && body->data.function.head->type == EXPR_SYMBOL
                      && body->data.function.head->data.symbol == fsym->data.symbol
                      && body->data.function.arg_count == m);
    if (unchanged) {
        for (size_t i = 0; i < m; i++) {
            if (!expr_eq(body->data.function.args[i], vars[i])) {
                unchanged = false;
                break;
            }
        }
    }
    if (unchanged) {
        expr_free(body);
        for (size_t i = 0; i < m; i++) expr_free(vars[i]);
        free(vars);
        return NULL;
    }

    /* Wrap the substituted body in Function[{t1,...,tm}, body] so we can
     * reuse the existing pure-function differentiation pipeline. */
    Expr** list_args = malloc(sizeof(Expr*) * m);
    for (size_t i = 0; i < m; i++) list_args[i] = expr_copy(vars[i]);
    Expr* var_list = mk_fnN_adopt("List", list_args, m);

    Expr* fn_args[2] = { var_list, body };
    Expr* pure_fn = expr_new_function(mk_sym("Function"), fn_args, 2);

    Expr* result = derivative_of_pure_function(deriv_head, pure_fn);
    expr_free(pure_fn);

    for (size_t i = 0; i < m; i++) expr_free(vars[i]);
    free(vars);

    return result;
}

/* ---------------------------------------------------------------------- */
/* Builtin: Derivative                                                     */
/* ---------------------------------------------------------------------- */

/*
 * The Derivative symbol primarily serves as a tag; the actual
 * differentiation happens inside the D dispatch. However, there is one
 * useful simplification we can perform at the head level:
 *
 *     Derivative[0, 0, ..., 0]  -->  Identity-like behaviour via the
 *     calling evaluator (the expression Derivative[0,..,0][f][args]
 *     is recognised by compute_deriv and does not need a builtin
 *     rewrite here).
 *
 * We therefore simply return NULL, leaving ``Derivative[n]`` in its
 * canonical unevaluated form. The value of this builtin is that it
 * lets us register attributes on the Derivative symbol.
 */
Expr* builtin_derivative(Expr* res) {
    (void)res;
    return NULL;
}

/* ---------------------------------------------------------------------- */
/* Registration                                                            */
/* ---------------------------------------------------------------------- */

void deriv_init(void) {
    symtab_add_builtin("D", builtin_d);
    symtab_add_builtin("Dt", builtin_dt);
    symtab_add_builtin("Derivative", builtin_derivative);

    /* Match the original deriv.m attribute set:
     *     SetAttributes[D, {Protected, ReadProtected}]
     *     SetAttributes[Dt, {Protected, ReadProtected}]
     *     SetAttributes[Derivative, {Protected, ReadProtected}] */
    symtab_get_def("D")->attributes          |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_get_def("Dt")->attributes         |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_get_def("Derivative")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
}
