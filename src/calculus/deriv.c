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

/* True if ``sym`` (an EXPR_SYMBOL) is listed in ``nonconsts`` (a
 * List[...] of symbols, or NULL). Used to recognise variables that the
 * caller has declared implicitly dependent on the differentiation
 * variable via the NonConstants option. */
static bool nonconsts_contains(const Expr* nonconsts, const Expr* sym) {
    if (!nonconsts || !sym) return false;
    if (sym->type != EXPR_SYMBOL) return false;
    if (nonconsts->type != EXPR_FUNCTION) return false;
    if (!is_sym(nonconsts->data.function.head, "List")) return false;
    size_t n = nonconsts->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        const Expr* v = nonconsts->data.function.args[i];
        if (v && v->type == EXPR_SYMBOL && v->data.symbol == sym->data.symbol) {
            return true;
        }
    }
    return false;
}

/* True if ``f`` contains any symbol listed in ``nonconsts``. With
 * NonConstants set we must NOT short-circuit constant subtrees that
 * happen to mention an implicitly-dependent symbol -- otherwise
 * D[2 y, x, NonConstants -> {y}] would wrongly fold to 0. */
static bool expr_contains_nonconst(const Expr* f, const Expr* nonconsts) {
    if (!nonconsts || nonconsts->type != EXPR_FUNCTION) return false;
    if (!is_sym(nonconsts->data.function.head, "List")) return false;
    size_t n = nonconsts->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        if (!expr_free_of(f, nonconsts->data.function.args[i])) return true;
    }
    return false;
}

/* Build the canonical unevaluated form D[sym, x, NonConstants -> {...}].
 * ``nonconsts`` is the already-canonicalised List of implicit-dependence
 * symbols; we deep-copy it into the result so the caller retains
 * ownership of the original. */
static Expr* build_unevaluated_d_nonconsts(const Expr* sym, const Expr* x,
                                           const Expr* nonconsts) {
    Expr* opt = mk_fn2("Rule",
                       mk_sym("NonConstants"),
                       expr_copy((Expr*)nonconsts));
    Expr** args = malloc(sizeof(Expr*) * 3);
    args[0] = expr_copy((Expr*)sym);
    args[1] = expr_copy((Expr*)x);
    args[2] = opt;
    Expr* r = expr_new_function(mk_sym("D"), args, 3);
    free(args);
    return r;
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
           s == SYM_GoldenAngle || s == SYM_Glaisher ||
           s == SYM_Khinchin ||
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

    /* --- Airy Ai: d/dg AiryAi[g] = AiryAiPrime[g]. --- */
    if (!strcmp(name, "AiryAi")) {
        return mk_fn1("AiryAiPrime", expr_copy(g));
    }

    /* --- Airy Ai': d/dg AiryAiPrime[g] = g AiryAi[g]  (from Ai'' = z Ai). --- */
    if (!strcmp(name, "AiryAiPrime")) {
        return mk_fn2("Times", expr_copy(g), mk_fn1("AiryAi", expr_copy(g)));
    }

    /* --- Airy Bi: d/dg AiryBi[g] = AiryBiPrime[g]. --- */
    if (!strcmp(name, "AiryBi")) {
        return mk_fn1("AiryBiPrime", expr_copy(g));
    }

    /* --- Airy Bi': d/dg AiryBiPrime[g] = g AiryBi[g]  (from Bi'' = z Bi). --- */
    if (!strcmp(name, "AiryBiPrime")) {
        return mk_fn2("Times", expr_copy(g), mk_fn1("AiryBi", expr_copy(g)));
    }

    /* --- error function: d/dg Erf[g] = (2/Sqrt[Pi]) E^(-g^2). --- */
    if (!strcmp(name, "Erf")) {
        Expr* coeff = mk_fn2("Times", mk_int(2),
                             mk_fn2("Power", mk_fn1("Sqrt", mk_sym("Pi")), mk_int(-1)));
        Expr* gaussian = mk_fn1("Exp",
                             mk_neg(mk_fn2("Power", expr_copy(g), mk_int(2))));
        return mk_fn2("Times", coeff, gaussian);
    }

    /* --- imaginary error function: d/dg Erfi[g] = (2/Sqrt[Pi]) E^(g^2). --- */
    if (!strcmp(name, "Erfi")) {
        Expr* coeff = mk_fn2("Times", mk_int(2),
                             mk_fn2("Power", mk_fn1("Sqrt", mk_sym("Pi")), mk_int(-1)));
        Expr* gaussian = mk_fn1("Exp",
                             mk_fn2("Power", expr_copy(g), mk_int(2)));
        return mk_fn2("Times", coeff, gaussian);
    }

    /* --- complementary error function: d/dg Erfc[g] = -(2/Sqrt[Pi]) E^(-g^2). --- */
    if (!strcmp(name, "Erfc")) {
        Expr* coeff = mk_fn2("Times", mk_int(-2),
                             mk_fn2("Power", mk_fn1("Sqrt", mk_sym("Pi")), mk_int(-1)));
        Expr* gaussian = mk_fn1("Exp",
                             mk_neg(mk_fn2("Power", expr_copy(g), mk_int(2))));
        return mk_fn2("Times", coeff, gaussian);
    }

    /* --- inverse error function: d/dg InverseErf[g] = (Sqrt[Pi]/2) E^(InverseErf[g]^2). --- */
    if (!strcmp(name, "InverseErf")) {
        Expr* coeff = mk_fn2("Times", mk_fn1("Sqrt", mk_sym("Pi")),
                             mk_fn2("Power", mk_int(2), mk_int(-1)));
        Expr* sq = mk_fn2("Power", mk_fn1("InverseErf", expr_copy(g)), mk_int(2));
        return mk_fn2("Times", coeff, mk_fn1("Exp", sq));
    }

    /* --- inverse complementary error function:
     *     d/dg InverseErfc[g] = -(Sqrt[Pi]/2) E^(InverseErfc[g]^2). --- */
    if (!strcmp(name, "InverseErfc")) {
        Expr* coeff = mk_fn2("Times", mk_int(-1),
                             mk_fn2("Times", mk_fn1("Sqrt", mk_sym("Pi")),
                                    mk_fn2("Power", mk_int(2), mk_int(-1))));
        Expr* sq = mk_fn2("Power", mk_fn1("InverseErfc", expr_copy(g)), mk_int(2));
        return mk_fn2("Times", coeff, mk_fn1("Exp", sq));
    }

    /* --- exponential integral: d/dg ExpIntegralEi[g] = E^g / g. --- */
    if (!strcmp(name, "ExpIntegralEi")) {
        return mk_fn2("Times", mk_fn1("Exp", expr_copy(g)),
                      mk_fn2("Power", expr_copy(g), mk_int(-1)));
    }

    /* --- logarithmic integral: d/dg LogIntegral[g] = 1 / Log[g]. --- */
    if (!strcmp(name, "LogIntegral")) {
        return mk_fn2("Power", mk_fn1("Log", expr_copy(g)), mk_int(-1));
    }

    /* --- exp/log and sqrt --- */
    if (!strcmp(name, "Exp")) return mk_fn1("Exp", expr_copy(g));
    if (!strcmp(name, "Log")) return mk_fn2("Power", expr_copy(g), mk_int(-1));
    if (!strcmp(name, "Sqrt")) {
        /* d/dg Sqrt[g] = 1/(2 Sqrt[g]) */
        return mk_fn2("Power", mk_fn2("Times", mk_int(2), mk_fn1("Sqrt", expr_copy(g))), mk_int(-1));
    }

    /* --- Fibonacci number: derivative w.r.t. its single argument. ---
     *   d/dn Fibonacci[n] = (1/Sqrt[5]) GoldenRatio^-n
     *       (GoldenRatio^(2 n) Log[GoldenRatio] + Cos[n Pi] Log[GoldenRatio]
     *        + Pi Sin[n Pi]). */
    if (!strcmp(name, "Fibonacci")) {
        Expr* logphi = mk_fn1("Log", mk_sym("GoldenRatio"));
        Expr* t1 = mk_fn2("Times",
                          mk_fn2("Power", mk_sym("GoldenRatio"),
                                 mk_fn2("Times", mk_int(2), expr_copy(g))),
                          expr_copy(logphi));
        Expr* t2 = mk_fn2("Times",
                          mk_fn1("Cos", mk_fn2("Times", expr_copy(g), mk_sym("Pi"))),
                          logphi);                          /* consumes logphi */
        Expr* t3 = mk_fn2("Times", mk_sym("Pi"),
                          mk_fn1("Sin", mk_fn2("Times", expr_copy(g), mk_sym("Pi"))));
        Expr** terms = malloc(sizeof(Expr*) * 3);
        terms[0] = t1; terms[1] = t2; terms[2] = t3;
        Expr* inner = mk_fnN_adopt("Plus", terms, 3);
        Expr* pref = mk_fn2("Times",
                            mk_fn2("Power", mk_fn1("Sqrt", mk_int(5)), mk_int(-1)),
                            mk_fn2("Power", mk_sym("GoldenRatio"),
                                   mk_fn2("Times", mk_int(-1), expr_copy(g))));
        return mk_fn2("Times", pref, inner);
    }

    /* --- LucasL number: derivative w.r.t. its single argument. ---
     *   d/dn LucasL[n] = GoldenRatio^-n (GoldenRatio^(2 n) Log[GoldenRatio]
     *       - Cos[n Pi] Log[GoldenRatio] - Pi Sin[n Pi]). */
    if (!strcmp(name, "LucasL")) {
        Expr* logphi = mk_fn1("Log", mk_sym("GoldenRatio"));
        Expr* t1 = mk_fn2("Times",
                          mk_fn2("Power", mk_sym("GoldenRatio"),
                                 mk_fn2("Times", mk_int(2), expr_copy(g))),
                          expr_copy(logphi));
        Expr* t2 = mk_neg(mk_fn2("Times",
                          mk_fn1("Cos", mk_fn2("Times", expr_copy(g), mk_sym("Pi"))),
                          logphi));                         /* consumes logphi */
        Expr* t3 = mk_neg(mk_fn2("Times", mk_sym("Pi"),
                          mk_fn1("Sin", mk_fn2("Times", expr_copy(g), mk_sym("Pi")))));
        Expr** terms = malloc(sizeof(Expr*) * 3);
        terms[0] = t1; terms[1] = t2; terms[2] = t3;
        Expr* inner = mk_fnN_adopt("Plus", terms, 3);
        Expr* pref = mk_fn2("Power", mk_sym("GoldenRatio"),
                            mk_fn2("Times", mk_int(-1), expr_copy(g)));
        return mk_fn2("Times", pref, inner);
    }

    return NULL; /* not a recognised elementary unary */
}

/* ---------------------------------------------------------------------- */
/* Core recursive derivative                                               */
/* ---------------------------------------------------------------------- */

/* Forward declarations. ``nonconsts`` may be NULL (no NonConstants
 * option in effect) or a fully canonicalised List[sym1, ...] of
 * symbols treated as implicit functions of x. */
static Expr* compute_deriv(Expr* f, Expr* x, Expr* nonconsts);

/* Shortcut: return D[g, x] as a fresh tree. When x is NULL, return
 * Dt[g] instead. For symbols and numeric atoms the answer is folded
 * immediately; for compound expressions we recurse through
 * compute_deriv so constants short-circuit there too. */
static Expr* deriv_of(Expr* g, Expr* x, Expr* nonconsts) {
    return compute_deriv(g, x, nonconsts);
}

/* The chain rule applied to an ``f[g1, g2, ..., gn]`` expression whose
 * head is not one of the explicitly-handled elementary heads. Produces
 *
 *   Sum_{k : D[gk, x] != 0} Derivative[0..1_k..0][f][g1..gn] * D[gk, x]
 *
 * If all partials vanish the result is 0. If exactly one contributes,
 * the outer Plus is elided. */
static Expr* chain_rule_unknown(Expr* f, Expr* x, Expr* nonconsts) {
    Expr* head = f->data.function.head;
    size_t n = f->data.function.arg_count;
    Expr** args = f->data.function.args;

    /* Collect (nonzero-derivative, partial-index) terms. */
    Expr** terms = malloc(sizeof(Expr*) * n);
    size_t nterms = 0;
    for (size_t k = 0; k < n; k++) {
        Expr* dk = deriv_of(args[k], x, nonconsts);
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
static Expr* deriv_of_derivative_form(Expr* f, Expr* x, Expr* nonconsts) {
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
        Expr* dk = deriv_of(args[k], x, nonconsts);
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
 * only constant-handling and the "unknown symbol" base-case differ.
 * ``nonconsts`` (may be NULL) names symbols that should be treated as
 * implicit functions of x; these short-circuit neither expr_free_of
 * nor the symbol-base-case. */
static Expr* compute_deriv(Expr* f, Expr* x, Expr* nonconsts) {
    /* --- Partial-derivative base cases. --- */
    if (x) {
        if (expr_eq(f, x)) return mk_int(1);
        /* D[Indeterminate, x] = Indeterminate: an indeterminate quantity has
         * no well-defined derivative. Matches Mathematica and keeps higher
         * derivatives of UnitStep (whose first derivative carries an
         * Indeterminate Piecewise value) stable rather than collapsing the
         * value to 0. */
        if (f->type == EXPR_SYMBOL && f->data.symbol == SYM_Indeterminate)
            return expr_copy(f);
        /* A symbol declared as NonConstant produces the canonical
         * unevaluated implicit-derivative form. */
        if (f->type == EXPR_SYMBOL && nonconsts_contains(nonconsts, f)) {
            return build_unevaluated_d_nonconsts(f, x, nonconsts);
        }
        /* Don't short-circuit `Equal[...]` or `Inequality[...]` even when
         * they have no x: the dispatch below distributes D over each
         * value slot and lets the outer evaluator simplify the residue,
         * matching Mathematica's semantics for D[a == b, x] and
         * D[a < b < c, x]. */
        bool f_is_equal_or_ineq = (f->type == EXPR_FUNCTION &&
                                   f->data.function.head->type == EXPR_SYMBOL &&
                                   (f->data.function.head->data.symbol == SYM_Equal
                                    || f->data.function.head->data.symbol == SYM_Inequality));
        if (!f_is_equal_or_ineq &&
            expr_free_of(f, x) && !expr_contains_nonconst(f, nonconsts)) {
            return mk_int(0);
        }
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
            for (size_t i = 0; i < n; i++) ts[i] = deriv_of(args[i], x, nonconsts);
            return mk_fnN_adopt("Plus", ts, n);
        }

        /* --- Equal (a == b [== c ...]): distribute the derivative
         *     through each side. Mathematica returns Equal[D[a,x],
         *     D[b,x], ...], collapsing arms whose derivative is 0
         *     via the outer evaluator. The same convention applies
         *     for partial, total, and NonConstants-tagged derivatives. */
        if (h == SYM_Equal) {
            Expr** ts = malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) ts[i] = deriv_of(args[i], x, nonconsts);
            return mk_fnN_adopt("Equal", ts, n);
        }

        /* --- Inequality[v0, op0, v1, op1, ...]: differentiate each value
         *     slot, keep the operator symbols verbatim. The outer
         *     evaluator may then collapse the residual (e.g. if all
         *     derivatives are 0, it becomes Inequality[0, op, 0, op, 0]
         *     and reduces in builtin_inequality). */
        if (h == SYM_Inequality && (n & 1u) == 1) {
            Expr** ts = malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) {
                if ((i & 1u) == 1) ts[i] = expr_copy(args[i]);    /* op symbol */
                else               ts[i] = deriv_of(args[i], x, nonconsts);
            }
            return mk_fnN_adopt("Inequality", ts, n);
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
                    factors[j] = (i == j) ? deriv_of(args[j], x, nonconsts)
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
            Expr* da = deriv_of(a, x, nonconsts);
            Expr* db = deriv_of(b, x, nonconsts);

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
            for (size_t i = 0; i < n; i++) ts[i] = deriv_of(args[i], x, nonconsts);
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
                    Expr* dbody = deriv_of(body, x, nonconsts);
                    Expr* dbody_eval = eval_and_free(dbody);
                    Expr* new_fn2 = mk_fn2("Function",
                        expr_copy(bvar), dbody_eval);
                    return mk_fn2("RootSum", expr_copy(fn1), new_fn2);
                }
                if (fa == 1) {
                    Expr* body = fn2->data.function.args[0];
                    Expr* dbody = deriv_of(body, x, nonconsts);
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
            Expr* r = compute_deriv(quot, x, nonconsts);
            expr_free(quot);
            return r;
        }

        /* --- HypergeometricPFQ[{a..}, {b..}, z]: derivative w.r.t. z is
         *   (prod a_i / prod b_j) pFq[{a_i+1}, {b_j+1}, z], times D[z, x].
         * Only the z-derivative is handled here; if the parameter lists
         * depend on x, defer to the generic chain rule (return NULL). */
        if (h == SYM_HypergeometricPFQ && n == 3
            && args[0]->type == EXPR_FUNCTION && args[1]->type == EXPR_FUNCTION
            && expr_free_of(args[0], x) && expr_free_of(args[1], x)) {
            Expr* al = args[0];
            Expr* bl = args[1];
            Expr* z  = args[2];
            Expr* dz = deriv_of(z, x, nonconsts);
            if (is_lit_zero(dz)) { expr_free(dz); return mk_int(0); }
            size_t p = al->data.function.arg_count;
            size_t q = bl->data.function.arg_count;

            /* prefactor = prod a_i / prod b_j */
            Expr* pref;
            {
                size_t w = 0;
                Expr** pf = malloc(sizeof(Expr*) * (p + 1));
                for (size_t i = 0; i < p; i++) pf[w++] = expr_copy(al->data.function.args[i]);
                if (q > 0) {
                    Expr** db = malloc(sizeof(Expr*) * q);
                    for (size_t j = 0; j < q; j++) db[j] = expr_copy(bl->data.function.args[j]);
                    Expr* bprod = mk_fnN_adopt("Times", db, q);
                    pf[w++] = mk_fn2("Power", bprod, mk_int(-1));
                }
                if (w == 0) { free(pf); pref = mk_int(1); }
                else pref = mk_fnN_adopt("Times", pf, w);
            }

            /* shifted parameter lists {a_i + 1}, {b_j + 1} */
            Expr** na = p ? malloc(sizeof(Expr*) * p) : NULL;
            for (size_t i = 0; i < p; i++)
                na[i] = mk_fn2("Plus", expr_copy(al->data.function.args[i]), mk_int(1));
            Expr* nal = expr_new_function(mk_sym("List"), na, p);
            free(na);
            Expr** nb = q ? malloc(sizeof(Expr*) * q) : NULL;
            for (size_t j = 0; j < q; j++)
                nb[j] = mk_fn2("Plus", expr_copy(bl->data.function.args[j]), mk_int(1));
            Expr* nbl = expr_new_function(mk_sym("List"), nb, q);
            free(nb);

            Expr* newp = expr_new_function(mk_sym("HypergeometricPFQ"),
                             (Expr*[]){ nal, nbl, expr_copy(z) }, 3);
            Expr** ff = malloc(sizeof(Expr*) * 3);
            ff[0] = pref; ff[1] = newp; ff[2] = dz;
            return mk_fnN_adopt("Times", ff, 3);
        }

        /* --- Fibonacci[A, B]: chain rule through both arguments.
         *   dF/dB = (2 A F[A-1,B] + (A-1) B F[A,B]) / (4 + B^2)
         *   dF/dA = (Pi + 2 ArcSinh[B/2] Cot[A Pi]) Csc[A Pi] F[-A,B]
         *           + (Pi Cot[A Pi]
         *              + ArcSinh[B/2] (Cot[A Pi]^2 + Csc[A Pi]^2)) F[A,B]
         * The total derivative is dF/dA D[A,x] + dF/dB D[B,x]; zero-derivative
         * arms are dropped (so D[F[n,x],x] keeps only the dF/dB term, and
         * D[F[n,x],n] keeps only dF/dA). */
        if (h == SYM_Fibonacci && n == 2) {
            Expr* A = args[0];
            Expr* B = args[1];
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* dB = deriv_of(B, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dB)) {
                Expr* fAm1 = mk_fn2("Fibonacci",
                                    mk_fn2("Plus", expr_copy(A), mk_int(-1)),
                                    expr_copy(B));
                Expr* fA = mk_fn2("Fibonacci", expr_copy(A), expr_copy(B));
                Expr* t1 = mk_fn2("Times",
                                  mk_fn2("Times", mk_int(2), expr_copy(A)), fAm1);
                Expr* t2 = mk_fn2("Times",
                                  mk_fn2("Times",
                                         mk_fn2("Plus", expr_copy(A), mk_int(-1)),
                                         expr_copy(B)),
                                  fA);
                Expr* numer = mk_fn2("Plus", t1, t2);
                Expr* denom = mk_fn2("Plus", mk_int(4),
                                     mk_fn2("Power", expr_copy(B), mk_int(2)));
                Expr* dFdB = mk_fn2("Times", numer,
                                    mk_fn2("Power", denom, mk_int(-1)));
                terms[nt++] = mk_fn2("Times", dFdB, dB);
            } else {
                expr_free(dB);
            }

            if (!is_lit_zero(dA)) {
                Expr* asinh = mk_fn1("ArcSinh",
                                     mk_fn2("Times", expr_copy(B),
                                            mk_fn2("Power", mk_int(2), mk_int(-1))));
                Expr* cot = mk_fn1("Cot",
                                   mk_fn2("Times", expr_copy(A), mk_sym("Pi")));
                Expr* csc = mk_fn1("Csc",
                                   mk_fn2("Times", expr_copy(A), mk_sym("Pi")));
                Expr* fnegA = mk_fn2("Fibonacci",
                                     mk_fn2("Times", mk_int(-1), expr_copy(A)),
                                     expr_copy(B));
                Expr* fA = mk_fn2("Fibonacci", expr_copy(A), expr_copy(B));

                /* coef1 = (Pi + 2 ArcSinh[B/2] Cot[A Pi]) Csc[A Pi] */
                Expr* c1inner = mk_fn2("Plus", mk_sym("Pi"),
                                       mk_fn2("Times",
                                              mk_fn2("Times", mk_int(2), expr_copy(asinh)),
                                              expr_copy(cot)));
                Expr* coef1 = mk_fn2("Times", c1inner, expr_copy(csc));

                /* coef2 = Pi Cot[A Pi] + ArcSinh[B/2] (Cot[A Pi]^2 + Csc[A Pi]^2) */
                Expr* cot2 = mk_fn2("Power", expr_copy(cot), mk_int(2));
                Expr* csc2 = mk_fn2("Power", csc, mk_int(2));
                Expr* c2b = mk_fn2("Times", asinh, mk_fn2("Plus", cot2, csc2));
                Expr* c2a = mk_fn2("Times", mk_sym("Pi"), cot);
                Expr* coef2 = mk_fn2("Plus", c2a, c2b);

                Expr* dFdA = mk_fn2("Plus",
                                    mk_fn2("Times", coef1, fnegA),
                                    mk_fn2("Times", coef2, fA));
                terms[nt++] = mk_fn2("Times", dFdA, dA);
            } else {
                expr_free(dA);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- LucasL[A, B]: chain rule through both arguments.
         *   dL/dB = (A (2 L[A-1,B] + B L[A,B])) / (4 + B^2)
         *   dL/dA = L[A+1,B] (2 ArcSinh[B/2] + Pi Tan[A Pi]) / Sqrt[4+B^2]
         *           + L[A,B] (-2 B ArcSinh[B/2]
         *                     - Pi (B + Sqrt[4+B^2]) Tan[A Pi]) / (2 Sqrt[4+B^2])
         * The total derivative is dL/dA D[A,x] + dL/dB D[B,x]; zero-derivative
         * arms are dropped (so D[L[n,x],x] keeps only the dL/dB term, and
         * D[L[n,x],n] keeps only dL/dA). */
        if (h == SYM_LucasL && n == 2) {
            Expr* A = args[0];
            Expr* B = args[1];
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* dB = deriv_of(B, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dB)) {
                Expr* lAm1 = mk_fn2("LucasL",
                                    mk_fn2("Plus", expr_copy(A), mk_int(-1)),
                                    expr_copy(B));
                Expr* lA = mk_fn2("LucasL", expr_copy(A), expr_copy(B));
                Expr* inner = mk_fn2("Plus",
                                     mk_fn2("Times", mk_int(2), lAm1),
                                     mk_fn2("Times", expr_copy(B), lA));
                Expr* numer = mk_fn2("Times", expr_copy(A), inner);
                Expr* denom = mk_fn2("Plus", mk_int(4),
                                     mk_fn2("Power", expr_copy(B), mk_int(2)));
                Expr* dLdB = mk_fn2("Times", numer,
                                    mk_fn2("Power", denom, mk_int(-1)));
                terms[nt++] = mk_fn2("Times", dLdB, dB);
            } else {
                expr_free(dB);
            }

            if (!is_lit_zero(dA)) {
                Expr* S = mk_fn1("Sqrt",
                                 mk_fn2("Plus", mk_int(4),
                                        mk_fn2("Power", expr_copy(B), mk_int(2))));
                Expr* asinh = mk_fn1("ArcSinh",
                                     mk_fn2("Times", expr_copy(B),
                                            mk_fn2("Power", mk_int(2), mk_int(-1))));
                Expr* tan = mk_fn1("Tan",
                                   mk_fn2("Times", expr_copy(A), mk_sym("Pi")));
                Expr* lAp1 = mk_fn2("LucasL",
                                    mk_fn2("Plus", expr_copy(A), mk_int(1)),
                                    expr_copy(B));
                Expr* lA = mk_fn2("LucasL", expr_copy(A), expr_copy(B));

                /* term1 = L[A+1,B] (2 ArcSinh[B/2] + Pi Tan[A Pi]) / Sqrt[4+B^2] */
                Expr* c1 = mk_fn2("Plus",
                                  mk_fn2("Times", mk_int(2), expr_copy(asinh)),
                                  mk_fn2("Times", mk_sym("Pi"), expr_copy(tan)));
                Expr* term1 = mk_fn2("Times",
                                     mk_fn2("Times", lAp1, c1),
                                     mk_fn2("Power", expr_copy(S), mk_int(-1)));

                /* term2 = L[A,B] (-2 B ArcSinh[B/2] - Pi (B + S) Tan[A Pi])
                 *         / (2 Sqrt[4+B^2]) */
                Expr* c2a = mk_neg(mk_fn2("Times",
                                   mk_fn2("Times", mk_int(2), expr_copy(B)),
                                   asinh));                  /* consumes asinh */
                Expr* c2b = mk_neg(mk_fn2("Times",
                                   mk_fn2("Times", mk_sym("Pi"),
                                          mk_fn2("Plus", expr_copy(B), S)),  /* consumes S */
                                   tan));                    /* consumes tan */
                Expr* c2 = mk_fn2("Plus", c2a, c2b);
                Expr* denom2 = mk_fn2("Times", mk_int(2),
                                      mk_fn1("Sqrt",
                                             mk_fn2("Plus", mk_int(4),
                                                    mk_fn2("Power", expr_copy(B),
                                                           mk_int(2)))));
                Expr* term2 = mk_fn2("Times",
                                     mk_fn2("Times", lA, c2),
                                     mk_fn2("Power", denom2, mk_int(-1)));

                Expr* dLdA = mk_fn2("Plus", term1, term2);
                terms[nt++] = mk_fn2("Times", dLdA, dA);
            } else {
                expr_free(dA);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- Zeta[S, A] (Hurwitz zeta): chain rule on both args.
         *   d/dA Zeta[S, A] = -S Zeta[1+S, A]   (elementary)
         *   d/dS Zeta[S, A] = Derivative[1,0][Zeta][S,A]  (generic; the
         *               derivative with respect to the order s has no
         *               elementary closed form).
         * Zero-derivative arms are dropped, so D[Zeta[s,a],a] keeps only the
         * dA term, matching Mathematica's -s Zeta[1+s, a]. */
        if (h == SYM_Zeta && n == 2) {
            Expr* S = args[0];
            Expr* A = args[1];
            Expr* dS = deriv_of(S, x, nonconsts);
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dA)) {
                /* -S Zeta[1+S, A] * dA */
                Expr* s1 = mk_fn2("Plus", mk_int(1), expr_copy(S));
                Expr* z  = mk_fn2("Zeta", s1, expr_copy(A));
                Expr* dZdA = mk_neg(mk_fn2("Times", expr_copy(S), z));
                terms[nt++] = mk_fn2("Times", dZdA, dA);
            } else {
                expr_free(dA);
            }

            if (!is_lit_zero(dS)) {
                /* Derivative[1, 0][Zeta][S, A] * dS */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("Zeta"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(S), expr_copy(A) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dS);
            } else {
                expr_free(dS);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- Zeta[S] (one-argument Riemann zeta): d/dS Zeta[S] has no
         * elementary closed form, so emit Derivative[1][Zeta][S] * D[S, x]. */
        if (h == SYM_Zeta && n == 1) {
            Expr* S = args[0];
            Expr* dS = deriv_of(S, x, nonconsts);
            if (is_lit_zero(dS)) { expr_free(dS); return mk_int(0); }
            Expr* op = expr_new_function(mk_sym("Derivative"),
                          (Expr*[]){ mk_int(1) }, 1);
            Expr* op_g = mk_fn_head1(op, mk_sym("Zeta"));
            Expr* applied = expr_new_function(op_g, (Expr*[]){ expr_copy(S) }, 1);
            return mk_fn2("Times", applied, dS);
        }

        /* --- PolyLog[N, Z] (polylogarithm): chain rule on both args.
         *   d/dZ PolyLog[N, Z] = PolyLog[N-1, Z]/Z   (elementary)
         *   d/dN PolyLog[N, Z] = Derivative[1,0][PolyLog][N, Z]  (no elementary
         *               closed form in the order N).
         * Zero-derivative arms are dropped, so D[PolyLog[n,z],z] keeps only the
         * dZ term, matching Mathematica's PolyLog[n-1, z]/z. */
        if (h == SYM_PolyLog && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* PolyLog[N-1, Z] Z^-1 * dZ */
                Expr* nm1 = mk_fn2("Plus", expr_copy(N), mk_int(-1));
                Expr* pl  = mk_fn2("PolyLog", nm1, expr_copy(Z));
                Expr* dPdZ = mk_fn2("Times", pl,
                                    mk_fn2("Power", expr_copy(Z), mk_int(-1)));
                terms[nt++] = mk_fn2("Times", dPdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1, 0][PolyLog][N, Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("PolyLog"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- BesselJ[N, Z]: chain rule on both args.
         *   d/dZ BesselJ[N, Z] = (BesselJ[N-1, Z] - BesselJ[N+1, Z]) / 2  (DLMF 10.6.1)
         *   d/dN BesselJ[N, Z] = Derivative[1,0][BesselJ][N,Z]  (the derivative
         *               with respect to the order has no elementary form).
         * Zero-derivative arms are dropped, so D[BesselJ[n,x],x] keeps only the
         * dZ term, matching Mathematica's 1/2 (BesselJ[-1+n,x]-BesselJ[1+n,x]). */
        if (h == SYM_BesselJ && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* (BesselJ[N-1, Z] - BesselJ[N+1, Z]) / 2 * dZ */
                Expr* bm1 = mk_fn2("BesselJ",
                              mk_fn2("Plus", expr_copy(N), mk_int(-1)), expr_copy(Z));
                Expr* bp1 = mk_fn2("BesselJ",
                              mk_fn2("Plus", expr_copy(N), mk_int(1)), expr_copy(Z));
                Expr* diff = mk_fn2("Plus", bm1, mk_neg(bp1));
                Expr* dBdZ = mk_fn2("Times",
                              mk_fn2("Power", mk_int(2), mk_int(-1)), diff);
                terms[nt++] = mk_fn2("Times", dBdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1,0][BesselJ][N,Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("BesselJ"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- BesselK[N, Z]: chain rule on both args.
         *   d/dZ BesselK[N, Z] = -(BesselK[N-1, Z] + BesselK[N+1, Z]) / 2  (DLMF 10.29.5)
         *   d/dN BesselK[N, Z] = Derivative[1,0][BesselK][N,Z]  (no elementary form).
         * Note the sign differs from BesselJ. */
        if (h == SYM_BesselK && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* -(BesselK[N-1, Z] + BesselK[N+1, Z]) / 2 * dZ */
                Expr* bm1 = mk_fn2("BesselK",
                              mk_fn2("Plus", expr_copy(N), mk_int(-1)), expr_copy(Z));
                Expr* bp1 = mk_fn2("BesselK",
                              mk_fn2("Plus", expr_copy(N), mk_int(1)), expr_copy(Z));
                Expr* sum = mk_fn2("Plus", bm1, bp1);
                Expr* dBdZ = mk_fn2("Times",
                              mk_fn2("Power", mk_int(-2), mk_int(-1)), sum);
                terms[nt++] = mk_fn2("Times", dBdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1,0][BesselK][N,Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("BesselK"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- BesselI[N, Z]: chain rule on both args.
         *   d/dZ BesselI[N, Z] = (BesselI[N-1, Z] + BesselI[N+1, Z]) / 2  (DLMF 10.29.5)
         *   d/dN BesselI[N, Z] = Derivative[1,0][BesselI][N,Z]  (no elementary form).
         * Like BesselK the two-term sum carries a '+', but (like BesselJ) the
         * overall coefficient is +1/2 rather than -1/2. */
        if (h == SYM_BesselI && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* (BesselI[N-1, Z] + BesselI[N+1, Z]) / 2 * dZ */
                Expr* bm1 = mk_fn2("BesselI",
                              mk_fn2("Plus", expr_copy(N), mk_int(-1)), expr_copy(Z));
                Expr* bp1 = mk_fn2("BesselI",
                              mk_fn2("Plus", expr_copy(N), mk_int(1)), expr_copy(Z));
                Expr* sum = mk_fn2("Plus", bm1, bp1);
                Expr* dBdZ = mk_fn2("Times",
                              mk_fn2("Power", mk_int(2), mk_int(-1)), sum);
                terms[nt++] = mk_fn2("Times", dBdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1,0][BesselI][N,Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("BesselI"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- Gamma[A, Z] (upper incomplete gamma): chain rule on both args.
         *   dGamma/dZ = -Z^(A-1) E^-Z    (elementary)
         *   dGamma/dA = Derivative[1,0][Gamma][A,Z]  (no closed form without
         *               PolyGamma; left as the generic partial).
         * Zero-derivative arms are dropped, so D[Gamma[a,z],z] keeps only the
         * dZ term, matching Mathematica's -E^-z z^(a-1). */
        if (h == SYM_Gamma && n == 2) {
            Expr* A = args[0];
            Expr* Z = args[1];
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* -Z^(A-1) E^-Z * dZ */
                Expr* zpow = mk_fn2("Power", expr_copy(Z),
                                    mk_fn2("Plus", expr_copy(A), mk_int(-1)));
                Expr* enz  = mk_fn1("Exp", mk_neg(expr_copy(Z)));
                Expr* dGdZ = mk_neg(mk_fn2("Times", zpow, enz));
                terms[nt++] = mk_fn2("Times", dGdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dA)) {
                /* Derivative[1,0][Gamma][A,Z] * dA */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("Gamma"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(A), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dA);
            } else {
                expr_free(dA);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- PolyGamma[N, Z] (the polygamma family): chain rule on both args.
         *   d/dZ PolyGamma[N, Z] = PolyGamma[N+1, Z]   (raises the order)
         *   d/dN PolyGamma[N, Z] = Derivative[1,0][PolyGamma][N,Z]  (generic; the
         *               derivative with respect to the order has no elementary form).
         * Zero-derivative arms are dropped, so D[PolyGamma[n,x],x] keeps only the
         * dZ term, matching Mathematica's PolyGamma[1+n, x]. */
        if (h == SYM_PolyGamma && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* PolyGamma[N+1, Z] * dZ */
                Expr* np1 = mk_fn2("Plus", expr_copy(N), mk_int(1));
                Expr* pg  = mk_fn2("PolyGamma", np1, expr_copy(Z));
                terms[nt++] = mk_fn2("Times", pg, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1, 0][PolyGamma][N, Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("PolyGamma"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- Beta[A, B] (Euler beta): chain rule on both arguments.
         *   d/dA Beta[A,B] = Beta[A,B] (PolyGamma[0,A] - PolyGamma[0,A+B])
         *   d/dB Beta[A,B] = Beta[A,B] (PolyGamma[0,B] - PolyGamma[0,A+B])
         * Higher derivatives compose automatically by re-differentiating the
         * Beta and PolyGamma factors. Zero-derivative arms are dropped. */
        if (h == SYM_Beta && n == 2) {
            Expr* A = args[0];
            Expr* B = args[1];
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* dB = deriv_of(B, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dA)) {
                Expr* psiA = mk_fn2("PolyGamma", mk_int(0), expr_copy(A));
                Expr* psiS = mk_fn2("PolyGamma", mk_int(0),
                                    mk_fn2("Plus", expr_copy(A), expr_copy(B)));
                Expr* fac  = mk_fn2("Times",
                                    mk_fn2("Beta", expr_copy(A), expr_copy(B)),
                                    mk_fn2("Subtract", psiA, psiS));
                terms[nt++] = mk_fn2("Times", fac, dA);
            } else { expr_free(dA); }

            if (!is_lit_zero(dB)) {
                Expr* psiB = mk_fn2("PolyGamma", mk_int(0), expr_copy(B));
                Expr* psiS = mk_fn2("PolyGamma", mk_int(0),
                                    mk_fn2("Plus", expr_copy(A), expr_copy(B)));
                Expr* fac  = mk_fn2("Times",
                                    mk_fn2("Beta", expr_copy(A), expr_copy(B)),
                                    mk_fn2("Subtract", psiB, psiS));
                terms[nt++] = mk_fn2("Times", fac, dB);
            } else { expr_free(dB); }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- Beta[Z, A, B] (incomplete beta): chain rule on all three args.
         *   d/dZ Beta[Z,A,B] = Z^(A-1) (1-Z)^(B-1)   (the integrand; elementary)
         *   d/dA, d/dB        = generic Derivative[0,1,0]/[0,0,1][Beta][Z,A,B]
         *               (no elementary form). Zero-derivative arms are dropped,
         * so D[Beta[z,a,b],z] keeps only Z^(a-1)(1-z)^(b-1). */
        if (h == SYM_Beta && n == 3) {
            Expr* Z = args[0];
            Expr* A = args[1];
            Expr* B = args[2];
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* dB = deriv_of(B, x, nonconsts);
            Expr* terms[3];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                Expr* za  = mk_fn2("Power", expr_copy(Z),
                                   mk_fn2("Plus", expr_copy(A), mk_int(-1)));
                Expr* omz = mk_fn2("Power",
                                   mk_fn2("Subtract", mk_int(1), expr_copy(Z)),
                                   mk_fn2("Plus", expr_copy(B), mk_int(-1)));
                terms[nt++] = mk_fn2("Times", mk_fn2("Times", za, omz), dZ);
            } else { expr_free(dZ); }

            if (!is_lit_zero(dA)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(1), mk_int(0) }, 3);
                Expr* opb = mk_fn_head1(op, mk_sym("Beta"));
                Expr* applied = expr_new_function(opb,
                              (Expr*[]){ expr_copy(Z), expr_copy(A), expr_copy(B) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dA);
            } else { expr_free(dA); }

            if (!is_lit_zero(dB)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(0), mk_int(1) }, 3);
                Expr* opb = mk_fn_head1(op, mk_sym("Beta"));
                Expr* applied = expr_new_function(opb,
                              (Expr*[]){ expr_copy(Z), expr_copy(A), expr_copy(B) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dB);
            } else { expr_free(dB); }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return expr_new_function(mk_sym("Plus"), terms, nt);
        }

        /* --- Beta[Z0, Z1, A, B] (generalized incomplete beta).
         *   d/dZ1 = Z1^(A-1) (1-Z1)^(B-1),  d/dZ0 = -Z0^(A-1) (1-Z0)^(B-1)
         *   d/dA, d/dB = generic Derivative[0,0,1,0]/[0,0,0,1][Beta][...]. */
        if (h == SYM_Beta && n == 4) {
            Expr* Z0 = args[0];
            Expr* Z1 = args[1];
            Expr* A  = args[2];
            Expr* B  = args[3];
            Expr* dZ0 = deriv_of(Z0, x, nonconsts);
            Expr* dZ1 = deriv_of(Z1, x, nonconsts);
            Expr* dA  = deriv_of(A,  x, nonconsts);
            Expr* dB  = deriv_of(B,  x, nonconsts);
            Expr* terms[4];
            size_t nt = 0;

            if (!is_lit_zero(dZ1)) {
                Expr* za  = mk_fn2("Power", expr_copy(Z1),
                                   mk_fn2("Plus", expr_copy(A), mk_int(-1)));
                Expr* omz = mk_fn2("Power",
                                   mk_fn2("Subtract", mk_int(1), expr_copy(Z1)),
                                   mk_fn2("Plus", expr_copy(B), mk_int(-1)));
                terms[nt++] = mk_fn2("Times", mk_fn2("Times", za, omz), dZ1);
            } else { expr_free(dZ1); }

            if (!is_lit_zero(dZ0)) {
                Expr* za  = mk_fn2("Power", expr_copy(Z0),
                                   mk_fn2("Plus", expr_copy(A), mk_int(-1)));
                Expr* omz = mk_fn2("Power",
                                   mk_fn2("Subtract", mk_int(1), expr_copy(Z0)),
                                   mk_fn2("Plus", expr_copy(B), mk_int(-1)));
                terms[nt++] = mk_fn2("Times", mk_neg(mk_fn2("Times", za, omz)), dZ0);
            } else { expr_free(dZ0); }

            if (!is_lit_zero(dA)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(0), mk_int(1), mk_int(0) }, 4);
                Expr* opb = mk_fn_head1(op, mk_sym("Beta"));
                Expr* applied = expr_new_function(opb,
                              (Expr*[]){ expr_copy(Z0), expr_copy(Z1),
                                         expr_copy(A), expr_copy(B) }, 4);
                terms[nt++] = mk_fn2("Times", applied, dA);
            } else { expr_free(dA); }

            if (!is_lit_zero(dB)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(0), mk_int(0), mk_int(1) }, 4);
                Expr* opb = mk_fn_head1(op, mk_sym("Beta"));
                Expr* applied = expr_new_function(opb,
                              (Expr*[]){ expr_copy(Z0), expr_copy(Z1),
                                         expr_copy(A), expr_copy(B) }, 4);
                terms[nt++] = mk_fn2("Times", applied, dB);
            } else { expr_free(dB); }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return expr_new_function(mk_sym("Plus"), terms, nt);
        }

        /* --- Gamma[A] (one-argument): d/dA Gamma[A] = Gamma[A] PolyGamma[0, A].
         * Repeated differentiation then composes via the product rule, e.g.
         * D[Gamma[z], {z, 2}] = Gamma[z] PolyGamma[0,z]^2 + Gamma[z] PolyGamma[1,z]. */
        if (h == SYM_Gamma && n == 1) {
            Expr* A = args[0];
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* g  = mk_fn1("Gamma", expr_copy(A));
            Expr* pg = mk_fn2("PolyGamma", mk_int(0), expr_copy(A));
            Expr* dGdA = mk_fn2("Times", g, pg);
            return mk_fn2("Times", dGdA, dA);
        }

        /* --- LogGamma[A]: d/dA LogGamma[A] = PolyGamma[0, A] (the digamma).
         * Higher derivatives then raise the PolyGamma order automatically. */
        if (h == SYM_LogGamma && n == 1) {
            Expr* A = args[0];
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* pg = mk_fn2("PolyGamma", mk_int(0), expr_copy(A));
            return mk_fn2("Times", pg, dA);
        }

        /* --- Piecewise[{{v1, c1}, ...}, default]: differentiate each value
         *     expression and the default; the conditions ride through
         *     unchanged. D[Piecewise[{{vi, ci}}, d], x]
         *       = Piecewise[{{D[vi, x], ci}}, D[d, x]].
         *     Without this rule a second derivative of UnitStep (whose first
         *     derivative is a Piecewise) falls through to the generic chain
         *     rule and produces garbage. The default is 0 when omitted. --- */
        if (h == SYM_Piecewise && n >= 1 &&
            args[0]->type == EXPR_FUNCTION &&
            args[0]->data.function.head->type == EXPR_SYMBOL &&
            args[0]->data.function.head->data.symbol == SYM_List) {
            Expr* pairs = args[0];
            size_t np = pairs->data.function.arg_count;
            Expr** new_pairs = malloc(sizeof(Expr*) * (np ? np : 1));
            bool ok = true;
            size_t built = 0;
            for (size_t i = 0; i < np; i++) {
                Expr* pr = pairs->data.function.args[i];
                if (!(pr->type == EXPR_FUNCTION &&
                      pr->data.function.head->type == EXPR_SYMBOL &&
                      pr->data.function.head->data.symbol == SYM_List &&
                      pr->data.function.arg_count == 2)) { ok = false; break; }
                /* Piecewise has HoldAll, so its branch values are not
                 * re-evaluated by the outer evaluator; reduce each derivative
                 * here (e.g. 2 x^(2-1) -> 2 x). Conditions ride through
                 * verbatim. */
                Expr* dv   = eval_and_free(deriv_of(pr->data.function.args[0], x, nonconsts));
                Expr* cond = expr_copy(pr->data.function.args[1]);
                new_pairs[built++] = mk_fn2("List", dv, cond);
            }
            if (ok) {
                Expr* pairs_list = mk_fnN_adopt("List", new_pairs, np);
                Expr* ddefault = (n >= 2) ? eval_and_free(deriv_of(args[1], x, nonconsts))
                                          : mk_int(0);
                return mk_fn2("Piecewise", pairs_list, ddefault);
            }
            for (size_t i = 0; i < built; i++) expr_free(new_pairs[i]);
            free(new_pairs);
        }

        /* --- UnitStep[a1, ..., an]: acts as the product of the per-argument
         *     UnitStep[ai], so the derivative follows the product rule. Each
         *     factor contributes UnitStep'[ai] = Piecewise[{{Indeterminate,
         *     ai == 0}}, 0], multiplied by D[ai, x] and the UnitStep of the
         *     remaining arguments. For n == 1 the trailing UnitStep factor is
         *     absent and the result collapses to the single Piecewise. --- */
        if (h == SYM_UnitStep && n >= 1) {
            Expr** terms = malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) {
                Expr* cond  = mk_fn2("Equal", expr_copy(args[i]), mk_int(0));
                Expr* pair  = mk_fn2("List", mk_sym("Indeterminate"), cond);
                Expr* pairs = mk_fn1("List", pair);
                Expr* pw    = mk_fn2("Piecewise", pairs, mk_int(0));
                Expr* dai   = deriv_of(args[i], x, nonconsts);

                size_t fc = 0;
                Expr** factors = malloc(sizeof(Expr*) * 3);
                factors[fc++] = pw;
                factors[fc++] = dai;
                if (n > 1) {
                    Expr** rest = malloc(sizeof(Expr*) * (n - 1));
                    size_t rc = 0;
                    for (size_t j = 0; j < n; j++)
                        if (j != i) rest[rc++] = expr_copy(args[j]);
                    factors[fc++] = mk_fnN_adopt("UnitStep", rest, n - 1);
                }
                terms[i] = mk_fnN_adopt("Times", factors, fc);
            }
            return mk_fnN_adopt("Plus", terms, n);
        }

        /* --- Known elementary unary function: F'(g) * D[g, x]. --- */
        if (n == 1) {
            Expr* fp = elementary_fprime(h, args[0]);
            if (fp) {
                Expr* dg = deriv_of(args[0], x, nonconsts);
                return mk_fn2("Times", fp, dg);
            }
            /* Unknown single-argument function -- standard chain rule:
             *     Derivative[1][f][g] * D[g, x]. */
            Expr* op = mk_fn1("Derivative", mk_int(1));         /* Derivative[1]      */
            Expr* op_f = mk_fn_head1(op, expr_copy(head));      /* Derivative[1][f]   */
            Expr* applied = mk_fn_head1(op_f, expr_copy(args[0])); /* Derivative[1][f][g] */
            Expr* dg = deriv_of(args[0], x, nonconsts);
            return mk_fn2("Times", applied, dg);
        }

        /* --- Unknown multi-argument function: full chain rule. --- */
        return chain_rule_unknown(f, x, nonconsts);
    }

    /* ------------------------------------------------------------------ */
    /* Head is itself a function. The only pattern we can reduce is       */
    /* Derivative[...][f][args...].                                       */
    /* ------------------------------------------------------------------ */
    if (head->type == EXPR_FUNCTION) {
        Expr* r = deriv_of_derivative_form(f, x, nonconsts);
        if (r) return r;

        /* Applied InterpolatingFunction object: differentiate by the generic
         * chain rule, producing Derivative[0..1_k..0][InterpolatingFunction[...]]
         * [g1..gn] * D[gk, x]. The evaluator then reduces the Derivative head
         * to a derivative-annotated InterpolatingFunction (see interp.c), so
         * D[ifun[x], x] and ifun'[x] agree. */
        if (is_sym(head->data.function.head, "InterpolatingFunction")) {
            return chain_rule_unknown(f, x, nonconsts);
        }
    }

    /* Give up -- caller keeps the expression unevaluated. */
    return NULL;
}

/* ---------------------------------------------------------------------- */
/* Higher-order: D[f, {x, n}]                                              */
/* ---------------------------------------------------------------------- */

/* Returns a new expression representing the n-th partial derivative of
 * f with respect to x. For n == 0 the input is deep-copied. */
static Expr* higher_order_partial(Expr* f, Expr* x, int64_t order,
                                  Expr* nonconsts) {
    if (order <= 0) return expr_copy(f);
    Expr* current = expr_copy(f);
    for (int64_t i = 0; i < order; i++) {
        Expr* nxt = compute_deriv(current, x, nonconsts);
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
 *     {{x1,...,xN}}      -- an array (vector) spec, order 1 (gradient)
 *     {{x1,...,xN}, n}   -- an array spec, integer order n
 * On success, *var is set to a NON-owned pointer into the spec.
 *   - For an integer order, *order is set to the integer and *order_expr
 *     is NULL.
 *   - For a symbolic order, *order is set to the sentinel -1 and
 *     *order_expr is a NON-owned pointer to the spec's symbolic order
 *     argument.
 *   - For an array spec, *is_array is set to true and *var points at the
 *     inner array List expression (so its args are the variables).
 * Returns false if the spec cannot be interpreted. */
static bool parse_var_spec(Expr* spec, Expr** var, int64_t* order,
                           Expr** order_expr, bool* is_array) {
    if (order_expr) *order_expr = NULL;
    if (is_array) *is_array = false;
    if (spec->type == EXPR_FUNCTION &&
        is_sym(spec->data.function.head, "List")) {
        size_t ac = spec->data.function.arg_count;
        Expr* first = (ac >= 1) ? spec->data.function.args[0] : NULL;
        bool first_is_list = (first && first->type == EXPR_FUNCTION &&
                              is_sym(first->data.function.head, "List"));

        /* Array specs: {{x1,...,xN}} or {{x1,...,xN}, n}. */
        if (first_is_list) {
            if (ac == 1) {
                if (is_array) *is_array = true;
                *var = first;
                *order = 1;
                return true;
            }
            if (ac == 2) {
                Expr* k = spec->data.function.args[1];
                if (k->type == EXPR_INTEGER) {
                    if (k->data.integer < 0) return false;
                    if (is_array) *is_array = true;
                    *var = first;
                    *order = k->data.integer;
                    return true;
                }
                /* Symbolic / non-integer order with array spec: not handled. */
                return false;
            }
            return false;
        }

        /* Scalar specs: {var, n} or {var, k}. Any other shape
         * (e.g., {x}, {x, y, z}) is a malformed multiple-derivative
         * specifier; reject so the caller can emit D::dvar and leave
         * the call unevaluated (matching Mathematica). */
        if (ac != 2) {
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

/* ---------------------------------------------------------------------- */
/* Array (Outer-style) derivatives                                         */
/* ---------------------------------------------------------------------- */

/* Forward declaration. */
static Expr* outer_d_step(Expr* f, Expr* array, Expr* nonconsts);

/* For a scalar (non-List) ``leaf_of_f`` and an ``array`` that may itself
 * be a nested List, build the tensor whose leaves are D[leaf_of_f, v]
 * for each leaf v of ``array``. Preserves ``array``'s nested structure.
 * When the inner D cannot be reduced, we leave a literal D[leaf, v]
 * wrapper so the outer evaluator can re-attempt later. */
static Expr* d_at_array_leaves(Expr* leaf_of_f, Expr* array, Expr* nonconsts) {
    if (array->type == EXPR_FUNCTION &&
        is_sym(array->data.function.head, "List")) {
        size_t n = array->data.function.arg_count;
        Expr** items = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            items[i] = d_at_array_leaves(leaf_of_f,
                                         array->data.function.args[i],
                                         nonconsts);
        }
        return mk_fnN_adopt("List", items, n);
    }
    Expr* d = compute_deriv(leaf_of_f, array, nonconsts);
    if (!d) {
        /* Fallback: emit an unevaluated D[leaf, var] so the outer
         * evaluator can take another pass. Matches MMA's behaviour for
         * heads with no known derivative rule. */
        return mk_fn2("D", expr_copy(leaf_of_f), expr_copy(array));
    }
    return d;
}

/* Compute one array-derivative step: First[Outer[D, {f}, array]].
 *
 * This is equivalent to Outer[D, f, array] (the {f} wrapper plus First
 * cancel cleanly), so we recurse into f at each List level, and at every
 * non-List leaf we lay down a copy of ``array``'s skeleton populated
 * with D[leaf, v_i] entries. The result therefore has shape
 *     shape(f) ++ shape(array)
 * matching Mathematica's definition. */
static Expr* outer_d_step(Expr* f, Expr* array, Expr* nonconsts) {
    if (f->type == EXPR_FUNCTION &&
        is_sym(f->data.function.head, "List")) {
        size_t n = f->data.function.arg_count;
        Expr** items = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            items[i] = outer_d_step(f->data.function.args[i], array, nonconsts);
        }
        return mk_fnN_adopt("List", items, n);
    }
    return d_at_array_leaves(f, array, nonconsts);
}

/* Repeated array derivative: D[f, {array, n}] == apply outer_d_step n
 * times. Between iterations we evaluate so that the next pass sees a
 * simplified expression (mirrors higher_order_partial). */
static Expr* array_higher_order(Expr* f, Expr* array, int64_t order,
                                Expr* nonconsts) {
    if (order <= 0) return expr_copy(f);
    Expr* current = expr_copy(f);
    for (int64_t i = 0; i < order; i++) {
        Expr* nxt = outer_d_step(current, array, nonconsts);
        Expr* reduced = evaluate(nxt);
        expr_free(nxt);
        expr_free(current);
        current = reduced;
    }
    return current;
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
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), parts, n);
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
        Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
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
            Expr* new_exp = expr_new_function(expr_new_symbol(SYM_Plus),
                (Expr*[]){ expr_copy(exp), mk_neg(expr_copy(k)) }, 2);
            Expr* new_pow = mk_fn2("Power", expr_copy(var), new_exp);
            Expr* product = expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ fp, new_pow }, 2);
            return eval_and_free(product);
        }

        /* Power[b, var] with b free of var:
         *   D[b^x, {x, k}] = b^x * Log[b]^k. */
        if (expr_eq(exp, var) && expr_free_of(base, var)) {
            Expr* same  = mk_fn2("Power", expr_copy(base), expr_copy(var));
            Expr* logb  = mk_fn1("Log", expr_copy(base));
            Expr* logbk = mk_fn2("Power", logb, expr_copy(k));
            Expr* product = expr_new_function(expr_new_symbol(SYM_Times),
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
            Expr* du_raw = compute_deriv(u, var, NULL);
            if (du_raw) {
                Expr* a = evaluate(du_raw);
                expr_free(du_raw);
                if (expr_free_of(a, var)) {
                    /* half_pi_k = k * Pi / 2 */
                    Expr* half_pi_k = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_copy(k),
                                   expr_new_symbol(SYM_Pi),
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
                            mk_fn2("Times", mk_int(-1), expr_new_symbol(SYM_I)),
                            expr_copy(u));
                        Expr* arg     = mk_fn2("Plus", half_pi_k, minus_I_u);
                        const char* trig_name = (h == SYM_Cosh) ? "Cos" : "Sin";
                        Expr* trig    = mk_fn1(trig_name, arg);
                        Expr* neg_I   = mk_fn2("Times", mk_int(-1),
                                               expr_new_symbol(SYM_I));
                        Expr* neg_I_k = mk_fn2("Power", neg_I, expr_copy(k));
                        if (h == SYM_Cosh) {
                            result = expr_new_function(expr_new_symbol(SYM_Times),
                                (Expr*[]){ ak, neg_I_k, trig }, 3);
                        } else {
                            result = expr_new_function(expr_new_symbol(SYM_Times),
                                (Expr*[]){ ak, expr_new_symbol(SYM_I),
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
    Expr* spec = expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){ expr_copy(var), expr_copy(k) }, 2);
    return expr_new_function(expr_new_symbol(SYM_D),
        (Expr*[]){ expr_copy(f), spec }, 2);
}

/* Detect Rule[NonConstants, val] / RuleDelayed[NonConstants, val].
 * Returns true when ``a`` is such an option, and writes a fresh
 * canonical List[sym1, ...] to ``*out`` (caller owns; we wrap a bare
 * single symbol into a singleton List to match Mathematica's output
 * canonicalisation). */
static bool parse_nonconsts_option(Expr* a, Expr** out) {
    if (!a || a->type != EXPR_FUNCTION) return false;
    if (a->data.function.arg_count != 2) return false;
    Expr* head = a->data.function.head;
    if (head->type != EXPR_SYMBOL) return false;
    if (head->data.symbol != SYM_Rule && head->data.symbol != SYM_RuleDelayed) {
        return false;
    }
    Expr* lhs = a->data.function.args[0];
    Expr* rhs = a->data.function.args[1];
    if (lhs->type != EXPR_SYMBOL || lhs->data.symbol != SYM_NonConstants) {
        return false;
    }
    if (rhs->type == EXPR_FUNCTION && is_sym(rhs->data.function.head, "List")) {
        *out = expr_copy(rhs);
    } else {
        *out = mk_fn1("List", expr_copy(rhs));
    }
    return true;
}

Expr* builtin_d(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;                 /* D[f] stays unevaluated */

    Expr* f = res->data.function.args[0];

    /* Split trailing arguments into option Rules and var-specs. Options
     * may appear interleaved (Mathematica conventionally places them at
     * the end, but we accept either order to match how users write
     * D[expr, x, NonConstants -> y]). */
    Expr** specs = malloc(sizeof(Expr*) * (argc - 1));
    size_t nspecs = 0;
    Expr* nonconsts = NULL;
    for (size_t i = 1; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        Expr* parsed = NULL;
        if (parse_nonconsts_option(a, &parsed)) {
            if (nonconsts) expr_free(nonconsts);
            nonconsts = parsed;             /* later wins, mirroring MMA */
            continue;
        }
        specs[nspecs++] = a;
    }
    if (nspecs == 0) {
        if (nonconsts) expr_free(nonconsts);
        free(specs);
        return NULL;                        /* D[f, opts...] stays unevaluated */
    }

    /* Fixed-point form D[sym, var, NonConstants -> {... sym ...}].
     * The internal compute_deriv emits exactly this shape when it hits
     * an implicit-dependence symbol; we return the canonical
     * List-wrapped form so that NonConstants -> y normalises to
     * NonConstants -> {y}. If the input is already canonical,
     * returning a structurally-equal result lets the outer fixed-point
     * detector stop without further work. */
    if (nonconsts && nspecs == 1 && f->type == EXPR_SYMBOL &&
        nonconsts_contains(nonconsts, f) &&
        specs[0]->type == EXPR_SYMBOL && !expr_eq(f, specs[0])) {
        Expr* r = build_unevaluated_d_nonconsts(f, specs[0], nonconsts);
        expr_free(nonconsts);
        free(specs);
        return r;
    }

    /* Sequentially apply each spec (handling mixed partials D[f, x, y, ...]). */
    Expr* current = expr_copy(f);
    for (size_t i = 0; i < nspecs; i++) {
        Expr* var = NULL;
        int64_t order = 1;
        Expr* order_expr = NULL;
        bool is_array = false;
        if (!parse_var_spec(specs[i], &var, &order, &order_expr, &is_array)) {
            emit_dvar_message(specs[i]);       /* MMA-style D::dvar */
            if (nonconsts) expr_free(nonconsts);
            free(specs);
            expr_free(current);
            return NULL;                       /* malformed spec */
        }
        Expr* stepped = NULL;
        if (is_array) {
            /* {{x1,...,xN}} or {{x1,...,xN}, n}: array (Outer) derivative. */
            stepped = array_higher_order(current, var, order, nonconsts);
        } else if (order >= 0) {
            stepped = higher_order_partial(current, var, order, nonconsts);
        } else {
            /* Symbolic-order spec. NonConstants is not threaded through
             * the closed-form symbolic-k path; it falls through to the
             * unevaluated form for any mixed case. */
            Expr* sym = compute_deriv_symbolic_order(current, var, order_expr);
            if (sym) stepped = sym;
            else stepped = build_unevaluated_d(current, var, order_expr);
        }
        expr_free(current);
        current = stepped;
    }

    if (nonconsts) expr_free(nonconsts);
    free(specs);
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
        return compute_deriv(f, NULL, NULL);  /* may be NULL to stay unevaluated */
    }

    /* Dt[f, var_specs...] is identical to D[f, var_specs...]: it gives
     * the partial derivative. Forward to the D path. */
    Expr* current = expr_copy(f);
    for (size_t i = 1; i < argc; i++) {
        Expr* var = NULL;
        int64_t order = 1;
        Expr* order_expr = NULL;
        bool is_array = false;
        if (!parse_var_spec(res->data.function.args[i], &var, &order, &order_expr,
                            &is_array)) {
            emit_dvar_message(res->data.function.args[i]);
            expr_free(current);
            return NULL;
        }
        Expr* stepped = NULL;
        if (is_array) {
            stepped = array_higher_order(current, var, order, NULL);
        } else if (order >= 0) {
            stepped = higher_order_partial(current, var, order, NULL);
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
            Expr* d = compute_deriv(current, vars[i], NULL);
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
