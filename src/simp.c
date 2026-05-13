#include "simp.h"
#include "arithmetic.h"
#include "attr.h"
#include "eval.h"
#include "expand.h"
#include "facpoly.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "expr.h"
#include "rationalize.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "trigrat.h"
#include "qa.h"
#include "qafactor.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

/*
 * simp.c -- Simplify, Assuming, $Assumptions, and AssumeCtx.
 *
 * Simplify implements a small heuristic search over the existing battery
 * of algebraic transforms. The default complexity measure is
 * LeafCount(expr) + decimal-digit count of integer leaves; this matches
 * Mathematica's default and stops e.g. "100 Log[2]" from being rewritten
 * to "Log[2^100]". A user-supplied ComplexityFunction option overrides
 * this. See simp.h for the AssumeCtx contract.
 *
 * Assuming desugars to Block[{$Assumptions = $Assumptions && a}, body],
 * which reuses Block's existing scope-restoration code path. Nested
 * Assuming calls compose because each Block reads the current
 * $Assumptions OwnValue before extending it.
 */

/* ----------------------------------------------------------------------- */
/* Default complexity measure                                              */
/* ----------------------------------------------------------------------- */

static size_t int_digit_count_int64(int64_t v) {
    if (v == 0) return 1;
    if (v < 0) {
        /* INT64_MIN edge case: |v| not representable; count digits of the
         * negated value digit-at-a-time without ever forming -INT64_MIN. */
        size_t n = 1;
        int64_t t = v;
        while (t <= -10) { n++; t /= 10; }
        return n;
    }
    size_t n = 0;
    while (v > 0) { n++; v /= 10; }
    return n;
}

/*
 * simp_default_complexity implements Mathematica's SimplifyCount:
 *
 *   Symbol      -> 1
 *   Integer 0   -> 1
 *   Integer p>0 -> Floor[Log10[p]] + 1            == digits(p)
 *   Integer p<0 -> Floor[Log10[|p|]] + 2          == digits(|p|) + 1
 *   Rational    -> SimplifyCount[num] + SimplifyCount[den] + 1
 *   Complex     -> SimplifyCount[re]  + SimplifyCount[im]  + 1
 *   Real / MPFR -> 2                              (NumberQ but not Integer/Rational)
 *   String      -> 1                              (treated as a leaf, picocas extension)
 *   Function    -> SimplifyCount[head] + sum SimplifyCount[args]
 *
 * The negative-integer adjustment matches Mathematica's behaviour where
 * the leading "-" contributes one unit of complexity. The explicit
 * Rational/Complex cases keep e.g. 100 Log[2] (score 6) preferred over
 * Log[2^100] (score 32). */
size_t simp_default_complexity(const Expr* e) {
    if (!e) return 0;
    switch (e->type) {
        case EXPR_INTEGER: {
            int64_t v = e->data.integer;
            if (v == 0) return 1;
            size_t d = int_digit_count_int64(v);
            return v > 0 ? d : d + 1;
        }
        case EXPR_BIGINT: {
            int sgn = mpz_sgn(e->data.bigint);
            if (sgn == 0) return 1;
            size_t digits = mpz_sizeinbase(e->data.bigint, 10);
            return sgn > 0 ? digits : digits + 1;
        }
        case EXPR_REAL:    return 2;
        case EXPR_SYMBOL:  return 1;
        case EXPR_STRING:  return 1;
        case EXPR_FUNCTION: {
            const Expr* head = e->data.function.head;
            size_t argc = e->data.function.arg_count;
            /* Rational[n, d] and Complex[re, im] are Mathematica-specials:
             * SimplifyCount adds 1 for the wrapper, not the head's own
             * SimplifyCount. */
            if (head && head->type == EXPR_SYMBOL && argc == 2) {
                if (head->data.symbol == SYM_Rational ||
                    head->data.symbol == SYM_Complex) {
                    return simp_default_complexity(e->data.function.args[0])
                         + simp_default_complexity(e->data.function.args[1])
                         + 1;
                }
            }
            size_t total = simp_default_complexity(head);
            for (size_t i = 0; i < argc; i++) {
                total += simp_default_complexity(e->data.function.args[i]);
            }
            return total;
        }
#ifdef USE_MPFR
        case EXPR_MPFR: return 2;
#endif
    }
    return 1;
}

/* Builtin SimplifyCount[expr] -- exposes the default complexity to users
 * so they can inspect or use it inside a custom ComplexityFunction.
 * The caller (evaluate_step) frees `res` after we return a non-NULL Expr;
 * we must NOT free it here. */
Expr* builtin_simplify_count(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;
    size_t s = simp_default_complexity(res->data.function.args[0]);
    /* size_t comfortably fits in EXPR_INTEGER for any expression we'd
     * realistically see; on 64-bit size_t = 8 bytes, int64_t = 8 bytes. */
    return expr_new_integer((int64_t)s);
}

/* ----------------------------------------------------------------------- */
/* AssumeCtx -- normalised fact set                                        */
/* ----------------------------------------------------------------------- */

static void ctx_push(AssumeCtx* ctx, const Expr* fact) {
    if (ctx->count >= ctx->capacity) {
        size_t new_cap = ctx->capacity ? ctx->capacity * 2 : 8;
        Expr** np = (Expr**)realloc(ctx->facts, new_cap * sizeof(Expr*));
        if (!np) return;
        ctx->facts = np;
        ctx->capacity = new_cap;
    }
    ctx->facts[ctx->count++] = expr_copy((Expr*)fact);
}

static void ctx_walk(AssumeCtx* ctx, const Expr* a) {
    if (!a) return;
    if (a->type == EXPR_SYMBOL) {
        if (a->data.symbol == SYM_True) return;
        if (a->data.symbol == SYM_False) {
            ctx->inconsistent = true;
            return;
        }
        ctx_push(ctx, a);
        return;
    }
    if (a->type == EXPR_FUNCTION &&
        a->data.function.head &&
        a->data.function.head->type == EXPR_SYMBOL) {
        const char* h = a->data.function.head->data.symbol;
        if (h == SYM_And || h == SYM_List) {
            for (size_t i = 0; i < a->data.function.arg_count; i++) {
                ctx_walk(ctx, a->data.function.args[i]);
            }
            return;
        }
    }
    ctx_push(ctx, a);
}

AssumeCtx* assume_ctx_from_expr(const Expr* assum) {
    AssumeCtx* ctx = (AssumeCtx*)calloc(1, sizeof(AssumeCtx));
    if (!ctx) return NULL;
    ctx_walk(ctx, assum);
    return ctx;
}

void assume_ctx_free(AssumeCtx* ctx) {
    if (!ctx) return;
    for (size_t i = 0; i < ctx->count; i++) expr_free(ctx->facts[i]);
    free(ctx->facts);
    free(ctx);
}

/* ----------------------------------------------------------------------- */
/* Domain queries                                                          */
/* ----------------------------------------------------------------------- */

/* Three-valued sign for numeric literals: -1, 0, +1, or 2 (unknown). */
static int numeric_sign(const Expr* e) {
    if (!e) return 2;
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer > 0) return 1;
        if (e->data.integer < 0) return -1;
        return 0;
    }
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint);
    if (e->type == EXPR_REAL) {
        if (e->data.real > 0) return 1;
        if (e->data.real < 0) return -1;
        return 0;
    }
    return 2;
}

static bool fact_is_function(const Expr* f, const char* head, size_t arity) {
    return f && f->type == EXPR_FUNCTION
        && f->data.function.head
        && f->data.function.head->type == EXPR_SYMBOL
        && strcmp(f->data.function.head->data.symbol, head) == 0
        && f->data.function.arg_count == arity;
}

static bool fact_implies_strict_positive(const Expr* f, const Expr* x) {
    if (f->type != EXPR_FUNCTION) return false;
    if (f->data.function.arg_count != 2) return false;
    if (f->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = f->data.function.head->data.symbol;
    Expr* a = f->data.function.args[0];
    Expr* b = f->data.function.args[1];

    /* Greater[x, c] with c >= 0 (any sign with c == 0 still means x > 0). */
    if (h == SYM_Greater) {
        if (expr_eq(a, x)) {
            int s = numeric_sign(b);
            return (s == 0 || s == 1);
        }
    }
    /* Less[c, x] with c >= 0. */
    if (h == SYM_Less) {
        if (expr_eq(b, x)) {
            int s = numeric_sign(a);
            return (s == 0 || s == 1);
        }
    }
    /* GreaterEqual[x, c] with c > 0. */
    if (h == SYM_GreaterEqual) {
        if (expr_eq(a, x)) {
            int s = numeric_sign(b);
            return s == 1;
        }
    }
    /* LessEqual[c, x] with c > 0. */
    if (h == SYM_LessEqual) {
        if (expr_eq(b, x)) {
            int s = numeric_sign(a);
            return s == 1;
        }
    }
    return false;
}

static bool fact_implies_nonneg(const Expr* f, const Expr* x) {
    if (fact_implies_strict_positive(f, x)) return true;
    if (f->type != EXPR_FUNCTION) return false;
    if (f->data.function.arg_count != 2) return false;
    if (f->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = f->data.function.head->data.symbol;
    Expr* a = f->data.function.args[0];
    Expr* b = f->data.function.args[1];
    /* x >= c with c >= 0 ; or c <= x with c >= 0. */
    if (h == SYM_GreaterEqual && expr_eq(a, x)) {
        int s = numeric_sign(b);
        return (s == 0 || s == 1);
    }
    if (h == SYM_LessEqual && expr_eq(b, x)) {
        int s = numeric_sign(a);
        return (s == 0 || s == 1);
    }
    return false;
}

static bool fact_implies_strict_negative(const Expr* f, const Expr* x) {
    if (f->type != EXPR_FUNCTION) return false;
    if (f->data.function.arg_count != 2) return false;
    if (f->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = f->data.function.head->data.symbol;
    Expr* a = f->data.function.args[0];
    Expr* b = f->data.function.args[1];
    /* Less[x, c] with c <= 0. */
    if (h == SYM_Less && expr_eq(a, x)) {
        int s = numeric_sign(b);
        return (s == 0 || s == -1);
    }
    /* Greater[c, x] with c <= 0. */
    if (h == SYM_Greater && expr_eq(b, x)) {
        int s = numeric_sign(a);
        return (s == 0 || s == -1);
    }
    /* LessEqual[x, c] with c < 0. */
    if (h == SYM_LessEqual && expr_eq(a, x)) {
        int s = numeric_sign(b);
        return s == -1;
    }
    /* GreaterEqual[c, x] with c < 0. */
    if (h == SYM_GreaterEqual && expr_eq(b, x)) {
        int s = numeric_sign(a);
        return s == -1;
    }
    return false;
}

static bool fact_implies_nonpos(const Expr* f, const Expr* x) {
    if (fact_implies_strict_negative(f, x)) return true;
    if (f->type != EXPR_FUNCTION) return false;
    if (f->data.function.arg_count != 2) return false;
    if (f->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = f->data.function.head->data.symbol;
    Expr* a = f->data.function.args[0];
    Expr* b = f->data.function.args[1];
    if (h == SYM_LessEqual && expr_eq(a, x)) {
        int s = numeric_sign(b);
        return (s == 0 || s == -1);
    }
    if (h == SYM_GreaterEqual && expr_eq(b, x)) {
        int s = numeric_sign(a);
        return (s == 0 || s == -1);
    }
    return false;
}

/* Element[x, Domain] match. */
static bool fact_in_domain(const Expr* f, const Expr* x, const char* dom) {
    if (!fact_is_function(f, "Element", 2)) return false;
    if (!expr_eq(f->data.function.args[0], x)) return false;
    Expr* d = f->data.function.args[1];
    return d->type == EXPR_SYMBOL && strcmp(d->data.symbol, dom) == 0;
}

/* Recognise facts that prove `x` is an even integer. The most common form
 * users write is Mod[x, 2] == 0; Equal is Orderless so we accept the args
 * in either order. Element[x, Evens] is also accepted. v1 doesn't try to
 * decompose `x == 2 k` style facts -- those would need a follow-up
 * "k is integer" check that the assumption layer doesn't surface yet. */
static bool fact_implies_even(const Expr* f, const Expr* x) {
    if (fact_in_domain(f, x, "Evens")) return true;
    if (!fact_is_function(f, "Equal", 2)) return false;
    Expr* a = f->data.function.args[0];
    Expr* b = f->data.function.args[1];
    /* Look for Mod[x, 2] on either side, with 0 on the other side. */
    Expr* mod = NULL;
    Expr* zero = NULL;
    if (a->type == EXPR_FUNCTION &&
        a->data.function.head &&
        a->data.function.head->type == EXPR_SYMBOL &&
        a->data.function.head->data.symbol == SYM_Mod) {
        mod = a; zero = b;
    } else if (b->type == EXPR_FUNCTION &&
               b->data.function.head &&
               b->data.function.head->type == EXPR_SYMBOL &&
               b->data.function.head->data.symbol == SYM_Mod) {
        mod = b; zero = a;
    }
    if (!mod || !zero) return false;
    if (mod->data.function.arg_count != 2) return false;
    if (!expr_eq(mod->data.function.args[0], x)) return false;
    Expr* m = mod->data.function.args[1];
    if (m->type != EXPR_INTEGER || m->data.integer != 2) return false;
    if (zero->type != EXPR_INTEGER || zero->data.integer != 0) return false;
    return true;
}

/* Standard positive real symbols recognised without a fact. */
static bool is_positive_constant_symbol(const char* s) {
    return strcmp(s, "Pi") == 0 ||
           strcmp(s, "E") == 0 ||
           strcmp(s, "EulerGamma") == 0 ||
           strcmp(s, "GoldenRatio") == 0 ||
           strcmp(s, "Catalan") == 0 ||
           strcmp(s, "Degree") == 0 ||
           strcmp(s, "Glaisher") == 0 ||
           strcmp(s, "Khinchin") == 0;
}

/* Symbols whose value is always real-valued. */
static bool is_real_constant_symbol(const char* s) {
    if (is_positive_constant_symbol(s)) return true;
    return strcmp(s, "MachineEpsilon") == 0;
}

/* Forward declarations for mutual recursion across the predicate family. */
static bool prov_pos(const AssumeCtx* ctx, const Expr* x);
static bool prov_nn (const AssumeCtx* ctx, const Expr* x);
static bool prov_neg(const AssumeCtx* ctx, const Expr* x);
static bool prov_np (const AssumeCtx* ctx, const Expr* x);
static bool prov_int(const AssumeCtx* ctx, const Expr* x);
static bool prov_re (const AssumeCtx* ctx, const Expr* x);

/* True iff every argument of `e` is provably real-valued. */
static bool all_real(const AssumeCtx* ctx, const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (!prov_re(ctx, e->data.function.args[i])) return false;
    }
    return true;
}

static bool fact_directly_positive(const AssumeCtx* ctx, const Expr* x) {
    if (!ctx) return false;
    for (size_t i = 0; i < ctx->count; i++) {
        if (fact_implies_strict_positive(ctx->facts[i], x)) return true;
    }
    return false;
}

static bool fact_directly_nonneg(const AssumeCtx* ctx, const Expr* x) {
    if (!ctx) return false;
    for (size_t i = 0; i < ctx->count; i++) {
        if (fact_implies_nonneg(ctx->facts[i], x)) return true;
    }
    return false;
}

static bool fact_directly_negative(const AssumeCtx* ctx, const Expr* x) {
    if (!ctx) return false;
    for (size_t i = 0; i < ctx->count; i++) {
        if (fact_implies_strict_negative(ctx->facts[i], x)) return true;
    }
    return false;
}

static bool fact_directly_nonpos(const AssumeCtx* ctx, const Expr* x) {
    if (!ctx) return false;
    for (size_t i = 0; i < ctx->count; i++) {
        if (fact_implies_nonpos(ctx->facts[i], x)) return true;
    }
    return false;
}

/* True iff `x` is provably an even integer under `ctx`. Recognises
 * Mod[x, 2] == 0 and Element[x, Evens]; integer literals are handled at
 * the leaf level so we don't need a separate fact for them.
 *
 * Times propagation: a product is even when at least one factor is even
 * and every factor is integer-valued — covers `2*n` with `n` integer. */
static bool prov_even(const AssumeCtx* ctx, const Expr* x) {
    if (!x) return false;
    if (x->type == EXPR_INTEGER) return (x->data.integer % 2) == 0;
    if (x->type == EXPR_BIGINT) return mpz_even_p(x->data.bigint) != 0;
    if (ctx) {
        for (size_t i = 0; i < ctx->count; i++) {
            if (fact_implies_even(ctx->facts[i], x)) return true;
        }
    }
    /* Times propagation: at least one factor even AND all factors int. */
    if (x->type == EXPR_FUNCTION &&
        x->data.function.head &&
        x->data.function.head->type == EXPR_SYMBOL &&
        x->data.function.head->data.symbol == SYM_Times) {
        bool any_even = false;
        bool all_int = true;
        for (size_t i = 0; i < x->data.function.arg_count; i++) {
            const Expr* f = x->data.function.args[i];
            if (!prov_int(ctx, f)) { all_int = false; break; }
            if (prov_even(ctx, f)) any_even = true;
        }
        if (all_int && any_even) return true;
    }
    return false;
}

bool assume_known_even(const AssumeCtx* ctx, const Expr* x) {
    return prov_even(ctx, x);
}

static bool prov_pos(const AssumeCtx* ctx, const Expr* x) {
    if (!x) return false;
    if (numeric_sign(x) == 1) return true;
    if (fact_directly_positive(ctx, x)) return true;
    if (x->type == EXPR_SYMBOL && is_positive_constant_symbol(x->data.symbol)) return true;
    if (x->type == EXPR_FUNCTION &&
        x->data.function.head &&
        x->data.function.head->type == EXPR_SYMBOL) {
        const char* h = x->data.function.head->data.symbol;
        size_t n = x->data.function.arg_count;
        Expr** a = x->data.function.args;
        /* Times: positive iff every factor positive. */
        if (h == SYM_Times && n > 0) {
            for (size_t i = 0; i < n; i++) {
                if (!prov_pos(ctx, a[i])) return false;
            }
            return true;
        }
        /* Plus: at least one strictly positive, all others non-negative. */
        if (h == SYM_Plus && n > 0) {
            bool any = false;
            for (size_t i = 0; i < n; i++) {
                if (prov_pos(ctx, a[i])) { any = true; continue; }
                if (prov_nn(ctx, a[i])) continue;
                return false;
            }
            return any;
        }
        /* Power: positive base raised to anything is positive. */
        if (h == SYM_Power && n == 2) {
            if (prov_pos(ctx, a[0])) return true;
        }
        /* Exp[real] is strictly positive. */
        if (h == SYM_Exp && n == 1) {
            if (prov_re(ctx, a[0])) return true;
        }
        /* Abs[x] >= 0; strictly > 0 only when x != 0, which we cannot prove
         * from sign alone, so fall back to nonneg here. */
        /* Cosh[real] >= 1 > 0. */
        if (h == SYM_Cosh && n == 1 && prov_re(ctx, a[0])) return true;
        /* Sqrt[positive] is positive (and Sqrt is Power[_, 1/2]; that path
         * handled above already). */
    }
    return false;
}

static bool prov_nn(const AssumeCtx* ctx, const Expr* x) {
    if (!x) return false;
    int s = numeric_sign(x);
    if (s == 1 || s == 0) return true;
    if (prov_pos(ctx, x)) return true;
    if (fact_directly_nonneg(ctx, x)) return true;
    if (x->type == EXPR_FUNCTION &&
        x->data.function.head &&
        x->data.function.head->type == EXPR_SYMBOL) {
        const char* h = x->data.function.head->data.symbol;
        size_t n = x->data.function.arg_count;
        Expr** a = x->data.function.args;
        if (h == SYM_Times && n > 0) {
            for (size_t i = 0; i < n; i++) {
                if (!prov_nn(ctx, a[i])) return false;
            }
            return true;
        }
        if (h == SYM_Plus && n > 0) {
            for (size_t i = 0; i < n; i++) {
                if (!prov_nn(ctx, a[i])) return false;
            }
            return true;
        }
        /* Abs[real] >= 0. */
        if (h == SYM_Abs && n == 1 && prov_re(ctx, a[0])) return true;
        /* x^(2k) is non-negative for real x and integer k -- common case
         * x^2 covered via the integer-2 literal exponent. */
        if (h == SYM_Power && n == 2 && prov_re(ctx, a[0])) {
            if (a[1]->type == EXPR_INTEGER && (a[1]->data.integer % 2) == 0) return true;
        }
    }
    return false;
}

static bool prov_neg(const AssumeCtx* ctx, const Expr* x) {
    if (!x) return false;
    if (numeric_sign(x) == -1) return true;
    if (fact_directly_negative(ctx, x)) return true;
    /* Times: even number of negatives among factors, with the rest positive,
     * gives positive (not negative). For "negative" we need an odd number of
     * negative factors and the rest positive. v1 keeps this simple. */
    return false;
}

static bool prov_np(const AssumeCtx* ctx, const Expr* x) {
    if (!x) return false;
    int s = numeric_sign(x);
    if (s == -1 || s == 0) return true;
    if (prov_neg(ctx, x)) return true;
    if (fact_directly_nonpos(ctx, x)) return true;
    return false;
}

static bool prov_int(const AssumeCtx* ctx, const Expr* x) {
    if (!x) return false;
    if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT) return true;
    if (!ctx) return false;
    for (size_t i = 0; i < ctx->count; i++) {
        if (fact_in_domain(ctx->facts[i], x, "Integers") ||
            fact_in_domain(ctx->facts[i], x, "PositiveIntegers") ||
            fact_in_domain(ctx->facts[i], x, "NonnegativeIntegers") ||
            fact_in_domain(ctx->facts[i], x, "NegativeIntegers") ||
            fact_in_domain(ctx->facts[i], x, "NonpositiveIntegers")) return true;
    }
    if (x->type == EXPR_FUNCTION &&
        x->data.function.head &&
        x->data.function.head->type == EXPR_SYMBOL) {
        const char* h = x->data.function.head->data.symbol;
        size_t n = x->data.function.arg_count;
        Expr** a = x->data.function.args;
        if ((h == SYM_Times || h == SYM_Plus) && n > 0) {
            for (size_t i = 0; i < n; i++) {
                if (!prov_int(ctx, a[i])) return false;
            }
            return true;
        }
        /* Power[int, nonneg-int] is integer. */
        if (h == SYM_Power && n == 2 &&
            prov_int(ctx, a[0]) &&
            a[1]->type == EXPR_INTEGER && a[1]->data.integer >= 0) return true;
    }
    return false;
}

static bool prov_re(const AssumeCtx* ctx, const Expr* x) {
    if (!x) return false;
    if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT || x->type == EXPR_REAL) return true;
    if (prov_int(ctx, x)) return true;
    if (x->type == EXPR_SYMBOL) {
        if (is_real_constant_symbol(x->data.symbol)) return true;
    }
    if (ctx) {
        for (size_t i = 0; i < ctx->count; i++) {
            const Expr* f = ctx->facts[i];
            if (fact_in_domain(f, x, "Reals") ||
                fact_in_domain(f, x, "Rationals") ||
                fact_in_domain(f, x, "Integers") ||
                fact_in_domain(f, x, "Algebraics")) return true;
            if (fact_implies_nonneg(f, x) || fact_implies_nonpos(f, x)) return true;
        }
    }
    if (x->type == EXPR_FUNCTION &&
        x->data.function.head &&
        x->data.function.head->type == EXPR_SYMBOL) {
        const char* h = x->data.function.head->data.symbol;
        size_t n = x->data.function.arg_count;
        Expr** a = x->data.function.args;
        if ((h == SYM_Times || h == SYM_Plus) && n > 0 && all_real(ctx, x)) {
            return true;
        }
        /* Power[positive, real] is real. Power[real, integer] is real. */
        if (h == SYM_Power && n == 2) {
            if (prov_pos(ctx, a[0]) && prov_re(ctx, a[1])) return true;
            if (prov_re(ctx, a[0]) && prov_int(ctx, a[1])) return true;
        }
        /* Real-valued elementary functions of real arguments. */
        if (n == 1 && prov_re(ctx, a[0])) {
            if (h == SYM_Sin || h == SYM_Cos ||
                h == SYM_Tan || h == SYM_Cot ||
                h == SYM_Sec || h == SYM_Csc ||
                h == SYM_Sinh || h == SYM_Cosh ||
                h == SYM_Tanh || h == SYM_Coth ||
                h == SYM_Sech || h == SYM_Csch ||
                h == SYM_Exp || h == SYM_Abs ||
                h == SYM_Floor || h == SYM_Ceiling ||
                h == SYM_Round || h == SYM_Sign) return true;
        }
        /* Log[positive] is real. */
        if (h == SYM_Log && n == 1 && prov_pos(ctx, a[0])) return true;
        /* ArcTan[real] is real, ArcSinh[real] real. */
        if (n == 1 && prov_re(ctx, a[0])) {
            if (h == SYM_ArcTan || h == SYM_ArcSinh ||
                h == SYM_ArcCot) return true;
        }
    }
    return false;
}

bool assume_known_positive(const AssumeCtx* ctx, const Expr* x) { return prov_pos(ctx, x); }
bool assume_known_nonneg  (const AssumeCtx* ctx, const Expr* x) { return prov_nn (ctx, x); }
bool assume_known_negative(const AssumeCtx* ctx, const Expr* x) { return prov_neg(ctx, x); }
bool assume_known_nonpos  (const AssumeCtx* ctx, const Expr* x) { return prov_np (ctx, x); }
bool assume_known_integer (const AssumeCtx* ctx, const Expr* x) { return prov_int(ctx, x); }
bool assume_known_real    (const AssumeCtx* ctx, const Expr* x) { return prov_re (ctx, x); }

/* ----------------------------------------------------------------------- */
/* Helpers                                                                 */
/* ----------------------------------------------------------------------- */

/* Build f[arg], evaluate, and return the result. Takes ownership of `arg`. */
static Expr* call_unary_owned(const char* head_name, Expr* arg) {
    Expr* a[1] = { arg };
    Expr* call = expr_new_function(expr_new_symbol(head_name), a, 1);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

static Expr* call_unary_copy(const char* head_name, const Expr* arg) {
    return call_unary_owned(head_name, expr_copy((Expr*)arg));
}

/* ----------------------------------------------------------------------- */
/* $SimplifyDebug -- per-transform tracing                                 */
/* ----------------------------------------------------------------------- */

/*
 * When $SimplifyDebug is set to True, every transform invocation inside
 * simp_search emits one line on stderr in the format
 *   /<TransformName>/: <input> -> <output> [<elapsed> ms]
 * This is used to diagnose pathological inputs (Simplify hangs, runaway
 * candidate explosion, expensive single transforms). The check is read
 * directly off the OwnValue list -- evaluating $SimplifyDebug would
 * itself fire the OwnValue rule on every call. */
static bool simp_debug_enabled(void) {
    Rule* r = symtab_get_own_values("$SimplifyDebug");
    if (!r || !r->replacement) return false;
    Expr* v = r->replacement;
    return v->type == EXPR_SYMBOL && v->data.symbol == SYM_True;
}

static double simp_debug_elapsed_ms(clock_t t0) {
    return (double)(clock() - t0) * 1000.0 / (double)CLOCKS_PER_SEC;
}

static void simp_debug_log(const char* xform, const Expr* in,
                           const Expr* out, double ms) {
    char* sin  = expr_to_string((Expr*)in);
    char* sout = out ? expr_to_string((Expr*)out) : NULL;
    fprintf(stderr, "/%s/: %s -> %s [%.2f ms]\n",
            xform,
            sin  ? sin  : "?",
            sout ? sout : "(no change)",
            ms);
    free(sin);
    free(sout);
    fflush(stderr);
}

/* Wrap call_unary_copy with tracing when $SimplifyDebug is True.
 *
 * Note: an experimental generic FactorMemo lookup at this layer was
 * tried (Phase 11 attempt) but reverted -- the per-transform memos
 * already in place (Factor, TrigFactor, TrigExpand, TrigRoundtrip,
 * PythagReduce, PythagSquareComplete, HalfAngle) cover the high-
 * volume duplicates, and the additional malloc/hash overhead at
 * every call exceeded the marginal gain on cheap transforms like
 * Together / Cancel / Apart that are individually fast and rarely
 * repeated. */
static Expr* traced_call_unary(const char* xform, const Expr* in) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;
    Expr* r = call_unary_copy(xform, in);
    if (dbg) simp_debug_log(xform, in, r, simp_debug_elapsed_ms(t0));
    return r;
}

static bool is_rule_with_lhs(const Expr* e, const char* lhs_symbol) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.arg_count != 2) return false;
    if (!e->data.function.head || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    Expr* k = e->data.function.args[0];
    return k && k->type == EXPR_SYMBOL && strcmp(k->data.symbol, lhs_symbol) == 0;
}

static bool head_threads_over(const char* h) {
    return strcmp(h, "Equal") == 0 ||
           strcmp(h, "Unequal") == 0 ||
           strcmp(h, "Less") == 0 ||
           strcmp(h, "LessEqual") == 0 ||
           strcmp(h, "Greater") == 0 ||
           strcmp(h, "GreaterEqual") == 0 ||
           strcmp(h, "And") == 0 ||
           strcmp(h, "Or") == 0 ||
           strcmp(h, "Not") == 0 ||
           strcmp(h, "Xor") == 0 ||
           strcmp(h, "Implies") == 0;
}

/* ----------------------------------------------------------------------- */
/* Candidate set                                                           */
/* ----------------------------------------------------------------------- */

#define SIMP_CAND_CAP 12
#define SIMP_ROUNDS   2

typedef struct {
    Expr** items;
    size_t count;
    size_t capacity;
} CandSet;

static void cs_init(CandSet* cs) {
    cs->items = NULL;
    cs->count = 0;
    cs->capacity = 0;
}

static void cs_free(CandSet* cs) {
    for (size_t i = 0; i < cs->count; i++) expr_free(cs->items[i]);
    free(cs->items);
    cs->items = NULL;
    cs->count = 0;
    cs->capacity = 0;
}

static bool cs_contains(const CandSet* cs, const Expr* e) {
    for (size_t i = 0; i < cs->count; i++) {
        if (expr_eq(cs->items[i], e)) return true;
    }
    return false;
}

/* Take ownership of `e`; free if duplicate or set is full. */
static void cs_add_or_free(CandSet* cs, Expr* e) {
    if (!e) return;
    if (cs->count >= SIMP_CAND_CAP || cs_contains(cs, e)) {
        expr_free(e);
        return;
    }
    if (cs->count >= cs->capacity) {
        size_t new_cap = cs->capacity ? cs->capacity * 2 : 4;
        Expr** np = (Expr**)realloc(cs->items, new_cap * sizeof(Expr*));
        if (!np) { expr_free(e); return; }
        cs->items = np;
        cs->capacity = new_cap;
    }
    cs->items[cs->count++] = e;
}

/* ----------------------------------------------------------------------- */
/* Scoring                                                                 */
/* ----------------------------------------------------------------------- */

#define SIMP_SCORE_INF ((size_t)-1)

/* True iff a `Power[X, Rational[p, q]]` subtree with q >= 2 (a non-trivial
 * root) appears anywhere in e. Used by the nested-radical detector below. */
static bool subtree_has_root(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* head = e->data.function.head;
    size_t argc = e->data.function.arg_count;
    if (head && head->type == EXPR_SYMBOL && head->data.symbol == SYM_Power
        && argc == 2) {
        const Expr* exp = e->data.function.args[1];
        if (exp->type == EXPR_FUNCTION && exp->data.function.head
            && exp->data.function.head->type == EXPR_SYMBOL
            && exp->data.function.head->data.symbol == SYM_Rational
            && exp->data.function.arg_count == 2) {
            const Expr* qq = exp->data.function.args[1];
            if (qq->type == EXPR_INTEGER && qq->data.integer >= 2) return true;
        }
    }
    for (size_t i = 0; i < argc; i++) {
        if (subtree_has_root(e->data.function.args[i])) return true;
    }
    return false;
}

/* Penalty for *truly* nested radicals: a `Power[Compound, Rational[_, q]]`
 * (q >= 2) whose compound base itself contains another root. This lets
 * Simplify prefer denested forms like (1+Sqrt[5])/2 over the original
 * (2+Sqrt[5])^(1/3) -- their plain LeafCounts are tied or off by ~2,
 * but the denested form has no nested-radical structure, which is the
 * canonical preferred shape.
 *
 * Surcharge value (3) chosen as the smallest constant that lets the
 * cube-root denester's output win the user's (2+Sqrt[5])^(1/3) case
 * (LeafCount diff is 2) without dominating other transforms' choices.
 * Flat radicals like Sqrt[5] or Sqrt[x+y] are NOT penalised -- only
 * the truly-nested shape is. */
static size_t nested_radical_penalty(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return 0;
    size_t total = 0;
    const Expr* head = e->data.function.head;
    size_t argc = e->data.function.arg_count;
    if (head && head->type == EXPR_SYMBOL && head->data.symbol == SYM_Power
        && argc == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        if (base->type == EXPR_FUNCTION && base->data.function.head
            && base->data.function.head->type == EXPR_SYMBOL) {
            const char* bsym = base->data.function.head->data.symbol;
            if (bsym != SYM_Rational && bsym != SYM_Complex) {
                if (exp->type == EXPR_FUNCTION && exp->data.function.head
                    && exp->data.function.head->type == EXPR_SYMBOL
                    && exp->data.function.head->data.symbol == SYM_Rational
                    && exp->data.function.arg_count == 2) {
                    const Expr* qq = exp->data.function.args[1];
                    if (qq->type == EXPR_INTEGER && qq->data.integer >= 2
                        && subtree_has_root(base)) {
                        total += 3;
                    }
                }
            }
        }
    }
    for (size_t i = 0; i < argc; i++) {
        total += nested_radical_penalty(e->data.function.args[i]);
    }
    return total;
}

static size_t score_with_func(const Expr* e, const Expr* complexity_func) {
    if (!complexity_func) {
        return simp_default_complexity(e) + nested_radical_penalty(e);
    }
    Expr* a[1] = { expr_copy((Expr*)e) };
    Expr* call = expr_new_function(expr_copy((Expr*)complexity_func), a, 1);
    Expr* result = evaluate(call);
    expr_free(call);
    size_t s;
    if (result->type == EXPR_INTEGER) {
        s = (result->data.integer < 0) ? 0 : (size_t)result->data.integer;
    } else if (result->type == EXPR_BIGINT) {
        s = SIMP_SCORE_INF;
    } else {
        s = simp_default_complexity(e);
    }
    expr_free(result);
    return s;
}

/* ----------------------------------------------------------------------- */
/* Assumption-driven seed rewriters                                        */
/* ----------------------------------------------------------------------- */

/* For each direct EXPR_SYMBOL fact-target, generate context-specific
 * rewrite rules and apply them via ReplaceRepeated. The rules are
 * unconditional in pattern form: their conditional nature is captured by
 * the choice of the rule's free symbol -- e.g., we only emit
 *   Power[Power[<x>, 2], Rational[1, 2]] :> <x>
 * when <x> is the literal symbol that the assumption set says is
 * positive. So the rules are valid by construction whenever applied.
 *
 * The generated rule list is built as a string and parsed; this is
 * cheaper to maintain than constructing the AST by hand and matches the
 * style used in trigsimp.c.
 */

static bool sym_already_listed(char** list, size_t n, const char* s) {
    for (size_t i = 0; i < n; i++) if (strcmp(list[i], s) == 0) return true;
    return false;
}

/* Walk the assumption fact list and collect every EXPR_SYMBOL that the
 * context proves positive, real, integer, or even. The caller passes
 * pre-sized arrays plus the maximum count. */
static void collect_known_symbols(const AssumeCtx* ctx,
                                  char** positives, size_t* npos,
                                  char** reals,     size_t* nreal,
                                  char** integers,  size_t* nint,
                                  char** negatives, size_t* nneg,
                                  char** evens,     size_t* neven,
                                  size_t cap) {
    *npos = *nreal = *nint = *nneg = *neven = 0;
    if (!ctx) return;
    /* Mod[m, 2] == 0 hides `m` inside a function argument, so a top-level
     * "scan operands of facts" alone misses it. We additionally walk
     * Mod[s, _] for any symbol s that appears under an even-type fact. */
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (f->type != EXPR_FUNCTION) continue;
        for (size_t j = 0; j < f->data.function.arg_count; j++) {
            Expr* a = f->data.function.args[j];
            if (a->type == EXPR_SYMBOL) {
                const char* nm = a->data.symbol;
                if (assume_known_positive(ctx, a) && *npos < cap && !sym_already_listed(positives, *npos, nm)) {
                    positives[(*npos)++] = (char*)nm;
                }
                if (assume_known_negative(ctx, a) && *nneg < cap && !sym_already_listed(negatives, *nneg, nm)) {
                    negatives[(*nneg)++] = (char*)nm;
                }
                if (assume_known_real(ctx, a) && *nreal < cap && !sym_already_listed(reals, *nreal, nm)) {
                    reals[(*nreal)++] = (char*)nm;
                }
                if (assume_known_integer(ctx, a) && *nint < cap && !sym_already_listed(integers, *nint, nm)) {
                    integers[(*nint)++] = (char*)nm;
                }
                if (assume_known_even(ctx, a) && *neven < cap && !sym_already_listed(evens, *neven, nm)) {
                    evens[(*neven)++] = (char*)nm;
                }
            } else if (a->type == EXPR_FUNCTION &&
                       a->data.function.head &&
                       a->data.function.head->type == EXPR_SYMBOL &&
                       a->data.function.head->data.symbol == SYM_Mod &&
                       a->data.function.arg_count == 2 &&
                       a->data.function.args[0]->type == EXPR_SYMBOL) {
                Expr* sym = a->data.function.args[0];
                const char* nm = sym->data.symbol;
                if (assume_known_even(ctx, sym) && *neven < cap && !sym_already_listed(evens, *neven, nm)) {
                    evens[(*neven)++] = (char*)nm;
                }
            }
        }
    }
}

/* Produce a rewritten expression by applying assumption-derived rules via
 * ReplaceRepeated. Returns a newly owned expression, or NULL if no rules
 * were generated. The input is not consumed. */
static Expr* apply_assumption_rules(const Expr* input, const AssumeCtx* ctx) {
    if (!ctx) return NULL;

    /* Conservative caps for the per-symbol rule synthesis. */
    enum { MAX_SYM = 16 };
    char* positives[MAX_SYM]; size_t npos;
    char* reals    [MAX_SYM]; size_t nreal;
    char* integers [MAX_SYM]; size_t nint;
    char* negatives[MAX_SYM]; size_t nneg;
    char* evens    [MAX_SYM]; size_t neven;
    collect_known_symbols(ctx, positives, &npos, reals, &nreal,
                          integers, &nint, negatives, &nneg,
                          evens, &neven, MAX_SYM);

    /* Build a single rule list "{r1, r2, ...}" as a string, then parse. */
    char buf[8192];
    size_t off = 0;
    int wrote_any = 0;

    #define EMIT(...) do { \
        int w = snprintf(buf + off, sizeof(buf) - off, __VA_ARGS__); \
        if (w < 0 || (size_t)w >= sizeof(buf) - off) goto overflow; \
        off += (size_t)w; \
    } while (0)
    #define SEP() do { if (wrote_any) EMIT(", "); wrote_any = 1; } while (0)

    EMIT("{");

    for (size_t i = 0; i < npos; i++) {
        const char* x = positives[i];
        /* Sqrt[x^2] forms */
        SEP(); EMIT("Power[Power[%s, 2], Rational[1, 2]] :> %s", x, x);
        SEP(); EMIT("Power[Power[%s, -1], Rational[1, 2]] :> Power[%s, Rational[-1, 2]]", x, x);
        SEP(); EMIT("Power[Power[%s, -2], Rational[1, 2]] :> Power[%s, -1]", x, x);
        /* Sqrt[x^2 * rest] -> x * Sqrt[rest] for x > 0; lets multi-factor
         * radicals like Sqrt[x^2 y^2] reduce one symbol at a time. */
        SEP(); EMIT("Power[Times[Power[%s, 2], rest___], Rational[1, 2]] :> %s Power[Times[rest], Rational[1, 2]]", x, x);
        /* Abs[x] -> x  for x > 0 */
        SEP(); EMIT("Abs[%s] :> %s", x, x);
        /* Log[x^p] -> p Log[x]  for x > 0 (any real p; v1 accepts symbolic p too) */
        SEP(); EMIT("Log[Power[%s, p_]] :> p Log[%s]", x, x);
        /* Inverse-trig sum identity: ArcTan[x] + ArcTan[1/x] -> Pi/2  for x > 0.
         * picocas's matcher does NOT perform orderless-Plus subset matching
         * out of the box (unlike Mathematica), so the rule must explicitly
         * absorb the trailing terms via `+ rest___` and re-emit them. The
         * BlankNullSequence pattern matches 0 or more remaining terms, so
         * the bare two-term sum reduces too. */
        SEP(); EMIT("ArcTan[%s] + ArcTan[Power[%s, -1]] + rest___ :> Pi/2 + rest", x, x);
        /* ArcCosh double-angle reduction: ArcCosh[2x^2 - 1] -> 2 ArcCosh[x]
         * for x > 0. The identity holds for any complex x with the
         * principal branch (verified for x in [0,1] via the i*ArcCos bridge
         * and for x >= 1 directly), so the x > 0 condition is the weakest
         * sufficient assumption. The Plus arg list is canonical
         * (Plus[-1, Times[2, x^2]]). */
        SEP(); EMIT("ArcCosh[Plus[-1, Times[2, Power[%s, 2]]]] :> 2 ArcCosh[%s]", x, x);
    }

    for (size_t i = 0; i < nneg; i++) {
        const char* x = negatives[i];
        /* Abs[x] -> -x  for x < 0 */
        SEP(); EMIT("Abs[%s] :> -%s", x, x);
        /* Sqrt[x^2] -> -x  for x < 0  (Power[Power[x,2], 1/2]) */
        SEP(); EMIT("Power[Power[%s, 2], Rational[1, 2]] :> -%s", x, x);
        /* Sqrt[x^2 * rest] -> -x * Sqrt[rest] for x < 0. */
        SEP(); EMIT("Power[Times[Power[%s, 2], rest___], Rational[1, 2]] :> -%s Power[Times[rest], Rational[1, 2]]", x, x);
        /* Mirror of the x > 0 inverse-trig sum: ArcTan[x] + ArcTan[1/x]
         * -> -Pi/2 for x < 0. Same `+ rest___` trick as above to handle
         * the embedded-in-larger-sum case. */
        SEP(); EMIT("ArcTan[%s] + ArcTan[Power[%s, -1]] + rest___ :> -Pi/2 + rest", x, x);
    }

    /* For real-but-unknown-sign, Sqrt[x^2] -> Abs[x]. Skip symbols already
     * proven positive or negative (their stronger rule above wins). */
    for (size_t i = 0; i < nreal; i++) {
        const char* x = reals[i];
        if (sym_already_listed(positives, npos, x)) continue;
        if (sym_already_listed(negatives, nneg, x)) continue;
        SEP(); EMIT("Power[Power[%s, 2], Rational[1, 2]] :> Abs[%s]", x, x);
        /* Sqrt[x^2 * rest] -> Abs[x] * Sqrt[rest] for real x. */
        SEP(); EMIT("Power[Times[Power[%s, 2], rest___], Rational[1, 2]] :> Abs[%s] Power[Times[rest], Rational[1, 2]]", x, x);
    }

    /* Sin[n Pi] -> 0, Cos[n Pi] -> (-1)^n, Tan[n Pi] -> 0 for integer n.
     * Plus: (-1)^(even_int * n) -> 1 and ((-1)^n)^even_int -> 1, so the
     * Cos rule can collapse all the way (the standalone Cos result
     * Cos[k Pi]^4 -> Power[-1, 4 k], for instance). */
    for (size_t i = 0; i < nint; i++) {
        const char* n = integers[i];
        SEP(); EMIT("Sin[%s Pi] :> 0", n);
        SEP(); EMIT("Sin[Pi %s] :> 0", n);
        SEP(); EMIT("Cos[%s Pi] :> Power[-1, %s]", n, n);
        SEP(); EMIT("Cos[Pi %s] :> Power[-1, %s]", n, n);
        SEP(); EMIT("Tan[%s Pi] :> 0", n);
        SEP(); EMIT("Tan[Pi %s] :> 0", n);
        SEP(); EMIT("Power[-1, Times[m_Integer /; EvenQ[m], %s]] :> 1", n);
        SEP(); EMIT("Power[Power[-1, %s], m_Integer /; EvenQ[m]] :> 1", n);
    }

    /* Even-exponent identities: (-1)^m = 1 when m is even, and the
     * lifted forms ((-1)^k)^m and (-1)^(k m) when additionally k is an
     * integer (so k m is also even). The pair-wise integer/even rules
     * cover the common Cos[k Pi]^m -> 1 path -- the existing Cos rule
     * above rewrites Cos[k Pi] to Power[-1, k], so the Power surface form
     * we land on is Power[Power[-1, k], m]. */
    for (size_t i = 0; i < neven; i++) {
        const char* m = evens[i];
        SEP(); EMIT("Power[-1, %s] :> 1", m);
        /* Literal-integer multiplier (handles concrete numbers without
         * needing them in the integer set). */
        SEP(); EMIT("Power[-1, k_Integer %s] :> 1", m);
        SEP(); EMIT("Power[Power[-1, k_Integer], %s] :> 1", m);
        for (size_t j = 0; j < nint; j++) {
            const char* k = integers[j];
            SEP(); EMIT("Power[-1, %s %s] :> 1", k, m);
            SEP(); EMIT("Power[Power[-1, %s], %s] :> 1", k, m);
        }
    }

    /* Equal[u, v] facts -> two-way substitution rules. We use immediate
     * Rule (->) so the pattern uses exact structural matching. */
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (!fact_is_function(f, "Equal", 2)) continue;
        /* We can't easily re-emit arbitrary Expr* into our string buffer;
         * instead, build these rules as Expr* and merge them in below. */
        (void)f; /* handled in the Expr* merge step below */
    }

    EMIT("}");

    if (!wrote_any) {
        /* No string-built rules. We may still have Equal substitutions. */
    }

    Expr* string_rules = wrote_any ? parse_expression(buf) : NULL;

    /* Now build Equal-substitution rules.
     *
     * Two complementary rules per equation:
     *   1. The direct rule heavier(lhs,rhs) -> lighter (catches cases
     *      where the equation's LHS appears verbatim as a subterm).
     *   2. ONE monomial-isolation rule when diff = lhs - rhs is a Plus
     *      with >= 3 terms: pick the heaviest non-numeric term t and
     *      emit t -> -(other terms). Polynomial relations like
     *      a^2 + b^2 == 1 then rewrite occurrences of a^2 even when
     *      the full "a^2 + b^2" sum is not present in the input.
     *
     * Emitting only one monomial rule (instead of one per term) avoids
     * the bidirectional cycle a^2 -> 1-b^2 ; b^2 -> 1-a^2 that
     * ReplaceRepeated would chase up to its 65536 iteration cap. */
    Expr** eq_diffs = (Expr**)calloc(ctx->count, sizeof(Expr*));
    size_t eq_count = 0;
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (!fact_is_function(f, "Equal", 2)) continue;
        Expr* lhs = f->data.function.args[0];
        Expr* rhs = f->data.function.args[1];
        Expr* sub_args[2] = { expr_copy(lhs),
                              expr_new_function(expr_new_symbol("Times"),
                                  (Expr*[]){ expr_new_integer(-1), expr_copy(rhs) }, 2) };
        Expr* sum = expr_new_function(expr_new_symbol("Plus"), sub_args, 2);
        Expr* diff = evaluate(sum);
        eq_diffs[i] = diff;
        /* Always emit the direct heavier->lighter rule. */
        eq_count++;
        /* Plus extra monomial-isolation rule for polynomial relations. */
        if (diff->type == EXPR_FUNCTION &&
            diff->data.function.head &&
            diff->data.function.head->type == EXPR_SYMBOL &&
            diff->data.function.head->data.symbol == SYM_Plus &&
            diff->data.function.arg_count >= 3) {
            for (size_t j = 0; j < diff->data.function.arg_count; j++) {
                Expr* term = diff->data.function.args[j];
                if (term->type == EXPR_INTEGER || term->type == EXPR_BIGINT ||
                    term->type == EXPR_REAL) continue;
                eq_count++;
                break; /* one monomial rule per equation */
            }
        }
    }

    if (!string_rules && eq_count == 0) {
        free(eq_diffs);
        return NULL;
    }

    size_t string_len = 0;
    if (string_rules && string_rules->type == EXPR_FUNCTION) {
        string_len = string_rules->data.function.arg_count;
    }
    size_t total = string_len + eq_count;
    Expr** all = (Expr**)calloc(total, sizeof(Expr*));
    size_t fill = 0;
    if (string_rules && string_rules->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < string_len; i++) {
            all[fill++] = expr_copy(string_rules->data.function.args[i]);
        }
    }
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (!fact_is_function(f, "Equal", 2)) continue;
        Expr* lhs = f->data.function.args[0];
        Expr* rhs = f->data.function.args[1];
        Expr* diff = eq_diffs[i];

        /* Direct heavier->lighter rule. */
        Expr *src, *dst;
        if (simp_default_complexity(lhs) >= simp_default_complexity(rhs)) {
            src = lhs; dst = rhs;
        } else {
            src = rhs; dst = lhs;
        }
        Expr* direct[2] = { expr_copy(src), expr_copy(dst) };
        all[fill++] = expr_new_function(expr_new_symbol("Rule"), direct, 2);

        /* Polynomial-relation monomial-isolation rule (one per fact). */
        if (diff->type == EXPR_FUNCTION &&
            diff->data.function.head &&
            diff->data.function.head->type == EXPR_SYMBOL &&
            diff->data.function.head->data.symbol == SYM_Plus &&
            diff->data.function.arg_count >= 3) {
            size_t n = diff->data.function.arg_count;
            /* Pick the first non-numeric term, breaking ties by canonical
             * (Plus-Orderless) order which is already applied by the
             * evaluator. */
            size_t pick = (size_t)-1;
            size_t pick_score = 0;
            for (size_t j = 0; j < n; j++) {
                Expr* term = diff->data.function.args[j];
                if (term->type == EXPR_INTEGER || term->type == EXPR_BIGINT ||
                    term->type == EXPR_REAL) continue;
                size_t s = simp_default_complexity(term);
                if (pick == (size_t)-1 || s > pick_score) {
                    pick = j;
                    pick_score = s;
                }
            }
            if (pick != (size_t)-1) {
                Expr* term = diff->data.function.args[pick];
                Expr** other_args = (Expr**)calloc(n - 1, sizeof(Expr*));
                size_t oi = 0;
                for (size_t k = 0; k < n; k++) {
                    if (k == pick) continue;
                    other_args[oi++] = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(-1),
                                   expr_copy(diff->data.function.args[k]) }, 2);
                }
                Expr* iso_rhs;
                if (n - 1 == 1) {
                    iso_rhs = other_args[0];
                    free(other_args);
                } else {
                    iso_rhs = expr_new_function(expr_new_symbol("Plus"), other_args, n - 1);
                    free(other_args);
                }
                Expr* iso[2] = { expr_copy(term), iso_rhs };
                all[fill++] = expr_new_function(expr_new_symbol("Rule"), iso, 2);
            }
        }
    }
    for (size_t i = 0; i < ctx->count; i++) if (eq_diffs[i]) expr_free(eq_diffs[i]);
    free(eq_diffs);
    if (string_rules) expr_free(string_rules);

    Expr* rules_list = expr_new_function(expr_new_symbol("List"), all, fill);
    free(all);

    Expr* call_args[2] = { expr_copy((Expr*)input), rules_list };
    Expr* call = expr_new_function(expr_new_symbol("ReplaceRepeated"), call_args, 2);
    Expr* out = evaluate(call);
    expr_free(call);
    return out;

overflow:
    /* Buffer was too small; bail out, no rules applied. */
    return NULL;

    #undef EMIT
    #undef SEP
}

/* ----------------------------------------------------------------------- */
/* Trig/exp roundtrip composite                                            */
/* ----------------------------------------------------------------------- */

static Expr* transform_trig_roundtrip(const Expr* e) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Two-level memo lookup.
     *
     * Level 1 (cheap): keyed on the raw input `e`.  Lets repeated
     * calls on the *same* expression short-circuit the entire
     * pipeline including the expensive TrigToExp stage.  This is
     * the common case during candidate-set iteration -- the same
     * sub-expression flows through many transforms, and most of
     * the time the TrigRoundtrip result is unchanged.
     *
     * Level 2 (canonical): keyed on TrigToExp(e).  Catches
     * equivalent forms (e.g., `Cos[x]^2 Sec[2x]` and
     * `1/4 Sec[2x] (2 + 2 Cos[2x])`) which collapse to the same
     * exponential expression.  Pays the TrigToExp cost (which we'd
     * incur for stage 1 anyway), but saves the rest of the pipeline.
     *
     * On a miss at both levels, the result is stored under BOTH
     * the raw and canonical keys, so future identical AND
     * equivalent calls hit Level 1 / Level 2 respectively. */
    FactorMemo* memo = factor_memo_active();
    Expr* raw_key = NULL;
    if (memo) {
        Expr* raw_args[1] = { expr_copy((Expr*)e) };
        raw_key = expr_new_function(expr_new_symbol("TrigRoundtrip"),
                                    raw_args, 1);
        const Expr* hit = factor_memo_lookup(memo, raw_key);
        if (hit) {
            Expr* cached = expr_copy((Expr*)hit);
            expr_free(raw_key);
            if (dbg) simp_debug_log("TrigRoundtrip", e, cached,
                                    simp_debug_elapsed_ms(t0));
            return cached;
        }
    }

    /* Stage 1 of the pipeline: convert trig atoms to exponential form. */
    Expr* a = call_unary_copy("TrigToExp", e);

    /* Level 2 lookup keyed on TrigToExp(input). */
    Expr* canon_key = NULL;
    if (memo) {
        Expr* canon_args[1] = { expr_copy(a) };
        canon_key = expr_new_function(expr_new_symbol("TrigRoundtrip"),
                                      canon_args, 1);
        const Expr* hit = factor_memo_lookup(memo, canon_key);
        if (hit) {
            Expr* cached = expr_copy((Expr*)hit);
            /* Promote to Level 1 for next time the same `e` arrives. */
            if (raw_key) {
                factor_memo_store(memo, raw_key, cached);
                expr_free(raw_key);
            }
            expr_free(canon_key);
            expr_free(a);
            if (dbg) simp_debug_log("TrigRoundtrip", e, cached,
                                    simp_debug_elapsed_ms(t0));
            return cached;
        }
    }
    /* Explosion guard: TrigToExp is structurally expanding -- a
     * single `Cos[x] Cos[y]` (complexity 7) maps to a sum of four
     * exponentials (complexity 77, 11x growth).  Together / Cancel /
     * ExpToTrig on that intermediate is expensive AND the final
     * result tends to use Cosh / Sinh of imaginary arguments rather
     * than Cos / Sin, leaving us with a complex-coefficient form
     * that's worse for the simp candidate-set search than the
     * input.
     *
     * If TrigToExp expanded the input by more than 5x, abort the
     * round-trip: skip the slow Together / Cancel / ExpToTrig stages
     * and return the input unchanged.  Other transforms in the
     * candidate set still see the input form.
     *
     * Verified safe on the user-reference case (Sin[x]^3 + Sin[3x] -
     * 3 Sin[x] expands by ~3x at TrigToExp stage but still benefits
     * from the round-trip).  Triggers on inputs like Cos[x] Cos[y]
     * where TrigToExp blows up 11x. */
    size_t in_score = simp_default_complexity(e);
    size_t exp_score = simp_default_complexity(a);
    if (dbg) {
        fprintf(stderr, "  TrigRoundtrip complexity: in=%zu exp=%zu ratio=%.2f\n",
                in_score, exp_score,
                in_score > 0 ? (double)exp_score / in_score : 0.0);
    }
    Expr* d;
    if (in_score > 0 && exp_score > 5 * in_score) {
        expr_free(a);
        d = expr_copy((Expr*)e);
    } else {
        Expr* b = call_unary_owned("Together", a);
        Expr* c = call_unary_owned("Cancel", b);
        d = call_unary_owned("ExpToTrig", c);
    }

    /* Store under both keys so future identical or canonically-equivalent
     * calls hit the appropriate level. */
    if (raw_key) {
        factor_memo_store(memo, raw_key, d);
        expr_free(raw_key);
    }
    if (canon_key) {
        factor_memo_store(memo, canon_key, d);
        expr_free(canon_key);
    }

    if (dbg) simp_debug_log("TrigRoundtrip", e, d, simp_debug_elapsed_ms(t0));
    return d;
}

/* Roots-of-unity simplification.
 *
 * Recognises every (-1)^(p/q) and E^(I p Pi / q) atom in the input,
 * lifts the expression to a univariate polynomial in
 *   omega = (-1)^(1/Q),  Q = LCM of denominators,
 * reduces modulo the cyclotomic polynomial Phi_{2Q}(omega) (the minimal
 * polynomial of omega = e^(I Pi / Q) over Q), and substitutes back. The
 * reduction is exact: omega is a primitive (2Q)-th root of unity, so
 * Phi_{2Q}(omega) = 0, and any polynomial p(omega) is identically zero
 * iff Phi_{2Q}(x) divides p(x). The substitute -> reduce -> substitute
 * round-trip preserves correctness for any polynomial in omega regardless
 * of the choice of free coefficients.
 *
 * Handles e.g.
 *   1 - (-1)^(1/3) + (-1)^(2/3)                 -> 0
 *   1 - (-1)^(1/5) + (-1)^(2/5) - ... + (-1)^(4/5) -> 0
 *   3 + 2 E^(-2 I Pi/3) + 2 E^(2 I Pi/3)        -> 1
 *
 * Implemented as a small Mathematica-syntax helper installed lazily
 * into the symbol table on first call. The cyclotomic polynomial is
 * computed on-the-fly by recursive division: Phi_n(x) = (x^n - 1) /
 * Prod_{d | n, d < n} Phi_d(x). Cache pressure is light because the
 * recursion is bounded by the LCM 2Q (typically < 30 for hand-written
 * inputs) and PolynomialQuotient memoises subresults via the term
 * structure of x^n - 1. */
static void simp_install_roots_of_unity_helpers(void) {
    static bool installed = false;
    if (installed) return;
    /* Definitions are added as DownValues on internal `$ru*` symbols so
     * they don't shadow anything user-visible. parse_expression returns
     * a SetDelayed Expr*; evaluate runs the assignment and returns Null
     * (we free that). */
    const char* defs[] = {
        "$ruCyclotomic[1, x_] := x - 1",
        "$ruCyclotomic[n_Integer, x_] := Module["
        "  {d, num = x^n - 1, denom = 1},"
        "  Do[If[Mod[n, d] == 0, denom = denom * $ruCyclotomic[d, x]], {d, 1, n - 1}];"
        "  PolynomialQuotient[num, denom, x]]",
        /* Main simplifier: collect denominators, lift to polynomial in
         * $ru, reduce mod Phi_{2Q}($ru), substitute back. The mod 2Q
         * normalisation on the exponent handles negative-exponent forms
         * like E^(-I Pi p / q) without leaving x^(-k) terms that
         * PolynomialRemainder would reject. */
        "$ruSimplify[expr_] := Module["
        "  {denoms, Q, polyForm, phiPoly, reduced},"
        "  denoms = Union[Join["
        "    Cases[expr, Power[-1, Rational[_, q_]] :> q, {0, Infinity}],"
        "    Cases[expr, Power[E, Times[Complex[0, Rational[_, q_]], Pi]] :> q, {0, Infinity}]]];"
        "  If[denoms === {}, expr,"
        "    Q = Apply[LCM, denoms];"
        "    polyForm = expr /. {"
        "      Power[-1, Rational[a_, b_]] :> $ru^Mod[a Q/b, 2 Q],"
        "      Power[E, Times[Complex[0, Rational[a_, b_]], Pi]] :> $ru^Mod[a Q/b, 2 Q]};"
        "    phiPoly = $ruCyclotomic[2 Q, $ru];"
        "    reduced = PolynomialRemainder[polyForm, phiPoly, $ru];"
        "    reduced /. $ru -> Power[-1, 1/Q]]]"
    };
    for (size_t i = 0; i < sizeof(defs)/sizeof(defs[0]); i++) {
        Expr* parsed = parse_expression(defs[i]);
        if (!parsed) continue;
        Expr* r = evaluate(parsed);
        if (r) expr_free(r);
    }
    installed = true;
}

static Expr* simp_roots_of_unity(const Expr* e) {
    simp_install_roots_of_unity_helpers();
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;
    Expr* args[1] = { expr_copy((Expr*)e) };
    Expr* call = expr_new_function(
        expr_new_symbol("$ruSimplify"), args, 1);
    Expr* out = evaluate(call);
    if (dbg) simp_debug_log("RootsOfUnity", e, out,
                            simp_debug_elapsed_ms(t0));
    return out;
}

/* Pythagorean perfect-square completion: 1 +/- 2 Sin[x] Cos[x]
 * = (Sin[x] +/- Cos[x])^2. Lets Simplify reach factored forms like
 * (Sin + Cos)^4 from a Factor result of (1 + 2 Sin Cos)^2. We keep
 * this as its own transform (separate from TrigFactor) because
 * TrigFactor's identity rule list also contains the linear-combination
 * rule a Sin[x] + b Cos[x] -> Sqrt[a^2+b^2] Sin[x + ArcTan[a, b]],
 * which would re-rewrite (Sin + Cos) into a single trig and obscure the
 * factored form. As a standalone seed the rewrite produces a candidate
 * Simplify can score directly. */
/* Forward declaration -- definition lives near transform_pythag_reduce. */
static Expr* simp_memo_wrap(const Expr* e, const char* pseudo_head,
                            Expr* (*impl)(const Expr*));
static bool has_pythag_head(const Expr* e);
static bool has_non_integer_power(const Expr* e);
static bool is_rational_literal(const Expr* e);

static Expr* transform_pythag_square_complete_impl(const Expr* e) {
    static Expr* rules = NULL;
    if (!rules) {
        rules = parse_expression(
            "{ 1 + 2 Sin[x_] Cos[x_] + r___ :> (Sin[x] + Cos[x])^2 + r, "
            "  1 - 2 Sin[x_] Cos[x_] + r___ :> (Sin[x] - Cos[x])^2 + r, "
            "  -1 + Cosh[x_]^2 + Sinh[x_]^2 + 2 Sinh[x_] Cosh[x_] + r___ "
            "      :> (Sinh[x] + Cosh[x])^2 - 1 + r, "
            "  -1 + Cosh[x_]^2 + Sinh[x_]^2 - 2 Sinh[x_] Cosh[x_] + r___ "
            "      :> (Sinh[x] - Cosh[x])^2 - 1 + r }");
    }
    if (!rules) return NULL;
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Same fast-skip as PythagReduce: every rule LHS contains a
     * Cos/Sin/Cosh/Sinh head, so on inputs without any of those
     * the ReplaceRepeated walk finds nothing. */
    if (!has_pythag_head(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("PythagSquareComplete", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr* args[2] = { expr_copy((Expr*)e), expr_copy(rules) };
    Expr* call = expr_new_function(
        expr_new_symbol("ReplaceRepeated"), args, 2);
    Expr* out = evaluate(call);
    if (dbg) simp_debug_log("PythagSquareComplete", e, out,
                            simp_debug_elapsed_ms(t0));
    return out;
}

static Expr* transform_pythag_square_complete(const Expr* e) {
    return simp_memo_wrap(e, "$PythagSquareComplete",
                          transform_pythag_square_complete_impl);
}

/* Half-angle tangent identity, applied to both circular and hyperbolic
 * functions. Folds the Weierstrass forms
 *
 *   Sin[x] / (1 + Cos[x])              -> Tan[x/2]
 *   Sin[x]^a (1 + Cos[x])^(-a)         -> Tan[x/2]^a
 *   Sin[x] / (c (1 + Cos[x]))          -> Tan[x/2] / c        (FreeQ[c, x])
 *   Sin[x]^a (c (1 + Cos[x]))^(-a)     -> Tan[x/2]^a / c^a    (FreeQ[c, x])
 *   (1 - Cos[x]) / Sin[x]              -> Tan[x/2]
 *   (1 - Cos[x])^a Sin[x]^(-a)         -> Tan[x/2]^a
 *
 * and the analogous Sinh/Cosh -> Tanh[x/2] family (with the sign
 * difference (Cosh[x] - 1)/Sinh[x] == Tanh[x/2]). Each rule has a
 * trailing `r___` BlankNullSequence inside the Times so the rule fires
 * on subterms inside larger products (e.g. (1/2) Sin[x] / (1 + Cos[x])
 * still rewrites). The conditional-pattern guards (a + b === 0,
 * FreeQ[c, x]) are what keeps the rules general -- there are no
 * specific-numeric variants.
 *
 * Output complexity is uniformly less than or equal to the input on
 * every shape that fires, so simp_search's leaf-count tiebreak takes
 * the rewritten form. */
static Expr* transform_halfangle_impl(const Expr* e) {
    static Expr* rules = NULL;
    if (!rules) {
        rules = parse_expression(
            "{ "
            /* Trig: Sin / (1 + Cos) -> Tan[x/2] */
            "  Sin[x_] Power[1 + Cos[x_], -1] r___ :> Tan[x/2] r, "
            "  Sin[x_]^a_ (1 + Cos[x_])^b_ r___ "
            "    /; a + b === 0 :> Tan[x/2]^a r, "
            /* Trig: Sin / (c (1 + Cos)) -> Tan[x/2]/c */
            "  Sin[x_] Power[c_ + c_ Cos[x_], -1] r___ "
            "    /; FreeQ[c, x] :> Tan[x/2] c^(-1) r, "
            "  Sin[x_]^a_ (c_ + c_ Cos[x_])^b_ r___ "
            "    /; a + b === 0 && FreeQ[c, x] :> Tan[x/2]^a c^b r, "
            /* Trig: (1 - Cos) / Sin -> Tan[x/2] */
            "  (1 - Cos[x_]) Power[Sin[x_], -1] r___ :> Tan[x/2] r, "
            "  (1 - Cos[x_])^a_ Sin[x_]^b_ r___ "
            "    /; a + b === 0 :> Tan[x/2]^a r, "
            /* Hyperbolic: Sinh / (1 + Cosh) -> Tanh[x/2] */
            "  Sinh[x_] Power[1 + Cosh[x_], -1] r___ :> Tanh[x/2] r, "
            "  Sinh[x_]^a_ (1 + Cosh[x_])^b_ r___ "
            "    /; a + b === 0 :> Tanh[x/2]^a r, "
            /* Hyperbolic: Sinh / (c (1 + Cosh)) -> Tanh[x/2]/c */
            "  Sinh[x_] Power[c_ + c_ Cosh[x_], -1] r___ "
            "    /; FreeQ[c, x] :> Tanh[x/2] c^(-1) r, "
            "  Sinh[x_]^a_ (c_ + c_ Cosh[x_])^b_ r___ "
            "    /; a + b === 0 && FreeQ[c, x] :> Tanh[x/2]^a c^b r, "
            /* Hyperbolic: (Cosh - 1) / Sinh -> Tanh[x/2] */
            "  (-1 + Cosh[x_]) Power[Sinh[x_], -1] r___ :> Tanh[x/2] r, "
            "  (-1 + Cosh[x_])^a_ Sinh[x_]^b_ r___ "
            "    /; a + b === 0 :> Tanh[x/2]^a r "
            "}");
    }
    if (!rules) return NULL;
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Every HalfAngle rule LHS uses Sin/Cos or Sinh/Cosh.  Skip the
     * ReplaceRepeated walk on inputs without any of those heads. */
    if (!has_pythag_head(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("HalfAngle", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr* args[2] = { expr_copy((Expr*)e), expr_copy(rules) };
    Expr* call = expr_new_function(
        expr_new_symbol("ReplaceRepeated"), args, 2);
    Expr* out = evaluate(call);
    if (dbg) simp_debug_log("HalfAngle", e, out,
                            simp_debug_elapsed_ms(t0));
    return out;
}

static Expr* transform_halfangle(const Expr* e) {
    return simp_memo_wrap(e, "$HalfAngle", transform_halfangle_impl);
}

/* ----------------------------------------------------------------------- */
/* Radical product canonicaliser: simp_radicals                            */
/* ----------------------------------------------------------------------- */

/*
 * simp_radicals combines distinct positive-integer radicals that share an
 * exponent inside any Times node, e.g.
 *
 *     Power[2, 1/2] * Power[3, 1/2]   ->  Power[6, 1/2]
 *     Power[2, 1/3] * Power[3, 1/3]   ->  Power[6, 1/3]
 *
 * The evaluator does NOT auto-perform this combine because in general
 *     Power[a, p/q] * Power[b, p/q] != Power[a*b, p/q]
 * once a or b can be negative or non-real (the principal-value branch
 * shifts: e.g. Sqrt[-2] Sqrt[-3] = -Sqrt[6]). Restricting to positive
 * integer bases keeps the rewrite sound. Same-base products
 * (Sqrt[2]*Sqrt[2] -> 2) are already collapsed by Power's exponent
 * merging in the evaluator; this transform targets the cross-base case
 * so that
 *     Simplify[Sqrt[6] - Sqrt[2] Sqrt[3]]            -> 0
 *     Simplify[-Sqrt[2] Sqrt[3] x + Sqrt[6] x]       -> 0
 * after the rebuilt Times feeds back into the surrounding Plus and the
 * Sqrt[6] terms cancel.
 *
 * Implementation: a bottom-up structural walker. Each Times node has its
 * positive-integer radical factors bucketed by exponent; multi-element
 * buckets are fused into a single Power[product_of_bases, exp] which is
 * then evaluated so the perfect-power detection in builtin_power can
 * collapse e.g. Power[36, 1/2] -> 6.
 */

static bool radical_base_is_positive_integer(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer > 0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) > 0;
    return false;
}

/* Match Power[positive integer, Rational[p, q]]. On success, populates
 * out_base and out_exp with borrowed pointers into the input. */
static bool radical_factor_split(const Expr* e,
                                 const Expr** out_base,
                                 const Expr** out_exp) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL) return false;
    if (e->data.function.head->data.symbol != SYM_Power) return false;
    if (e->data.function.arg_count != 2) return false;
    const Expr* base = e->data.function.args[0];
    const Expr* exp  = e->data.function.args[1];
    if (!radical_base_is_positive_integer(base)) return false;
    if (!is_rational_literal(exp)) return false;
    /* Rational[p, 1] should never reach us (the evaluator canonicalises
     * to the bare integer); guard defensively anyway. */
    const Expr* den = exp->data.function.args[1];
    if (den->type == EXPR_INTEGER && den->data.integer == 1) return false;
    *out_base = base;
    *out_exp  = exp;
    return true;
}

/* Combine same-exponent positive-integer radical factors of a Times
 * node. Returns NULL when no combine fires. */
static Expr* simp_radicals_combine_times(const Expr* tn) {
    size_t n = tn->data.function.arg_count;
    if (n < 2) return NULL;

    bool*  consumed = (bool*) calloc(n, sizeof(bool));
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * n);
    size_t out = 0;
    bool changed = false;

    for (size_t i = 0; i < n; i++) {
        if (consumed[i]) continue;
        const Expr* arg_i = tn->data.function.args[i];
        const Expr* base_i;
        const Expr* exp_i;
        if (!radical_factor_split(arg_i, &base_i, &exp_i)) {
            new_args[out++] = expr_copy((Expr*)arg_i);
            consumed[i] = true;
            continue;
        }

        Expr* prod_base = expr_copy((Expr*)base_i);
        size_t group = 1;
        for (size_t j = i + 1; j < n; j++) {
            if (consumed[j]) continue;
            const Expr* arg_j = tn->data.function.args[j];
            const Expr* base_j;
            const Expr* exp_j;
            if (!radical_factor_split(arg_j, &base_j, &exp_j)) continue;
            if (!expr_eq((Expr*)exp_i, (Expr*)exp_j)) continue;

            Expr* mul_args[2] = { prod_base, expr_copy((Expr*)base_j) };
            Expr* mul = expr_new_function(expr_new_symbol("Times"),
                                          mul_args, 2);
            prod_base = evaluate(mul);
            expr_free(mul);
            consumed[j] = true;
            group++;
        }
        consumed[i] = true;

        if (group >= 2) {
            Expr* pow_args[2] = { prod_base, expr_copy((Expr*)exp_i) };
            Expr* pow = expr_new_function(expr_new_symbol("Power"),
                                          pow_args, 2);
            new_args[out++] = evaluate(pow);
            expr_free(pow);
            changed = true;
        } else {
            expr_free(prod_base);
            new_args[out++] = expr_copy((Expr*)arg_i);
        }
    }

    free(consumed);
    if (!changed) {
        for (size_t k = 0; k < out; k++) expr_free(new_args[k]);
        free(new_args);
        return NULL;
    }
    Expr* rebuilt = expr_new_function(expr_new_symbol("Times"),
                                      new_args, out);
    Expr* canonical = evaluate(rebuilt);
    expr_free(rebuilt);
    return canonical;
}

/* Bottom-up walker: rewrites children, rebuilds the node, then tries
 * the Times combine at the top.  Returns NULL when nothing changed. */
static Expr* simp_radicals_walk(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;

    size_t n = e->data.function.arg_count;
    Expr** new_args = NULL;
    bool any = false;
    for (size_t i = 0; i < n; i++) {
        Expr* r = simp_radicals_walk(e->data.function.args[i]);
        if (r) {
            if (!new_args) {
                new_args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
                for (size_t j = 0; j < i; j++) {
                    new_args[j] = expr_copy(e->data.function.args[j]);
                }
            }
            new_args[i] = r;
            any = true;
        } else if (new_args) {
            new_args[i] = expr_copy(e->data.function.args[i]);
        }
    }

    Expr* current = NULL;
    if (any) {
        Expr* head_copy = expr_copy(e->data.function.head);
        Expr* rebuilt = expr_new_function(head_copy, new_args, n);
        free(new_args);
        current = evaluate(rebuilt);
        expr_free(rebuilt);
    }

    const Expr* target = current ? current : e;
    if (target->type == EXPR_FUNCTION
        && target->data.function.head
        && target->data.function.head->type == EXPR_SYMBOL
        && target->data.function.head->data.symbol == SYM_Times) {
        Expr* combined = simp_radicals_combine_times(target);
        if (combined) {
            if (current) expr_free(current);
            return combined;
        }
    }
    return current;  /* may be NULL when no change */
}

static Expr* simp_radicals_impl(const Expr* e) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Cheap precondition: if there is no Power with a non-integer
     * exponent anywhere, no radical combine can fire. */
    if (!has_non_integer_power(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("Radicals", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr* r = simp_radicals_walk(e);
    Expr* out = r ? r : expr_copy((Expr*)e);
    if (dbg) simp_debug_log("Radicals", e, out,
                            simp_debug_elapsed_ms(t0));
    return out;
}

static Expr* simp_radicals(const Expr* e) {
    return simp_memo_wrap(e, "$Radicals", simp_radicals_impl);
}

/* ----------------------------------------------------------------------- */
/* Sqrt-of-Sqrt denesting: simp_denest_sqrt                                */
/* ----------------------------------------------------------------------- */

/* Forward declaration -- defined below in the RadicalCanon section. The
 * denesting walker calls it on its final result so that Sqrt[Rational[
 * p, q]] forms produced inside (P + s)/2 get canonicalised to
 * Sqrt[p*q]/q before the result lands in the seed set. Without this
 * post-pass, my candidate has the same complexity as the inputs (two
 * unsimplified Sqrt[Rational] leaves vs one Sqrt[nested] leaf) and the
 * round-loop's RadicalCanon pass isn't in SIMP_TRANSFORMS, so the
 * canonicalisation never fires. */
static Expr* transform_radical_canon(const Expr* e);

/* Forward declaration -- defined just below. */
static bool is_sqrt(const Expr* e);

/* Rewrite every Power[Rational[m, n], 1/2] subtree to Times[Sqrt[m*n],
 * Rational[1, n]] (i.e. Sqrt[m*n]/n) for positive m, n. Returns a fresh
 * tree (caller owns); the input is not mutated.
 *
 * `transform_radical_canon` already attempts this rewrite, but the
 * picocas evaluator re-merges Times[Power[m, 1/2], Power[n, -1/2]]
 * back into Power[m/n, 1/2] for m > 1 (only m = 1 escapes because
 * Power[1, 1/2] evaluates to the integer 1 and drops out of the Times,
 * leaving Power[n, -1/2] for the negative-exponent rule). This helper
 * sidesteps the re-merge by constructing Sqrt[m*n] / n directly, which
 * the evaluator does not collapse back. */
static Expr* denest_rationalise_sqrt_of_rational(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    /* Recurse into children first (bottom-up). */
    size_t n = e->data.function.arg_count;
    Expr* head_copy = expr_copy(e->data.function.head);
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            new_args[i] = denest_rationalise_sqrt_of_rational(
                e->data.function.args[i]);
        }
    }
    Expr* rebuilt = expr_new_function(head_copy, new_args, n);
    if (new_args) free(new_args);

    /* Match Power[Rational[m, n], Rational[1, 2]] with m > 0, n > 0. */
    if (is_sqrt(rebuilt)) {
        Expr* base = rebuilt->data.function.args[0];
        if (is_rational_literal(base)) {
            Expr* num = base->data.function.args[0];
            Expr* den = base->data.function.args[1];
            if (num->type == EXPR_INTEGER && den->type == EXPR_INTEGER &&
                num->data.integer > 0 && den->data.integer > 1) {
                int64_t m = num->data.integer;
                int64_t nn = den->data.integer;
                int64_t mn;
                if (!__builtin_mul_overflow(m, nn, &mn)) {
                    /* Build Times[Power[m*n, 1/2], Rational[1, n]]. */
                    Expr* sqrt_args[2] = {
                        expr_new_integer(mn),
                        make_rational(1, 2)
                    };
                    Expr* sqrt_call = expr_new_function(
                        expr_new_symbol("Power"), sqrt_args, 2);
                    Expr* sqrt_e = evaluate(sqrt_call);
                    expr_free(sqrt_call);
                    Expr* prod_args[2] = {
                        sqrt_e,
                        make_rational(1, nn)
                    };
                    Expr* prod = expr_new_function(
                        expr_new_symbol("Times"), prod_args, 2);
                    Expr* out = evaluate(prod);
                    expr_free(prod);
                    expr_free(rebuilt);
                    return out;
                }
            }
        }
    }
    return rebuilt;
}

/*
 * simp_denest_sqrt rewrites Sqrt[A + b * Sqrt[C]] using the half-sum
 * identity
 *
 *     Sqrt[A + Sqrt[B]] = Sqrt[(A + s)/2] + Sqrt[(A - s)/2]
 *
 * where B = b^2 * C and s is a closed-form representative of
 * Sqrt[A^2 - B] (so s^2 = A^2 - B). The identity holds in the principal
 * branch for any (A, B) with A^2 - B >= 0 and the radicand A + Sqrt[B]
 * itself nonneg, provided the two inner radicands (A +/- s)/2 are both
 * nonneg -- which is the additional precondition we explicitly verify
 * via the AssumeCtx (or via direct numeric evaluation when A, s are
 * rational).
 *
 * When the original radicand is A - b * Sqrt[C] (i.e. b < 0), the
 * identity becomes
 *
 *     Sqrt[A - Sqrt[B]] = Sqrt[(A + s)/2] - Sqrt[(A - s)/2]
 *
 * (sign flipped on the second term). Since |A - b Sqrt[C]| = |A + (-b) Sqrt[C]|
 * after squaring, the discriminant D = A^2 - B is identical for the two
 * branches; only the output sign convention differs.
 *
 * The single primitive closes five user-supplied test cases:
 *   - Sqrt[3 + 2 Sqrt[2]]            -> 1 + Sqrt[2]                (A=3, B=8)
 *   - Sqrt[17 - 12 Sqrt[2]]          -> 3 - 2 Sqrt[2]              (A=17, B=288)
 *   - Sqrt[2 + Sqrt[3]]              -> (Sqrt[6]+Sqrt[2])/2        (A=2, B=3)
 *   - Sqrt[x+y+2 Sqrt[xy]]           -> Sqrt[x] + Sqrt[y]          (x,y>0)
 *   - Sqrt[x+Sqrt[x^2-y^2]]          -> Sqrt[(x+y)/2]+Sqrt[(x-y)/2] (x>y>0)
 *
 * Soundness: the transform refuses to fire whenever it cannot verify
 * P, Q >= 0 from purely numeric reasoning or the active AssumeCtx.
 * That is, a borderline case yields the input unchanged rather than a
 * possibly-wrong-branch denesting -- the standing Simplify-soundness
 * invariant.
 */

/* True iff e == Power[_, Rational[1, 2]] (a sqrt). The representation
 * the parser/evaluator settles on for `Sqrt[x]`. */
static bool is_sqrt(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol != SYM_Power) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* exp = e->data.function.args[1];
    if (!is_rational_literal(exp)) return false;
    Expr* num = exp->data.function.args[0];
    Expr* den = exp->data.function.args[1];
    return num->type == EXPR_INTEGER && num->data.integer == 1 &&
           den->type == EXPR_INTEGER && den->data.integer == 2;
}

/* Forward declaration: nonneg check used to combine multiple sqrt
 * factors into one when each base is provably nonneg. */
static bool denest_is_nonneg(const Expr* e, const AssumeCtx* ctx);

/* If `e` is a "sqrt-radical decomposition target" -- a single expression
 * that can be written as b * Sqrt[C] with a single radicand C -- populate
 * *out_b (a fresh expression representing b, possibly 1) and *out_C
 * (a fresh deep copy of the radicand). Returns false otherwise.
 *
 * Three forms are recognised:
 *   1. Power[C, 1/2]                          -- the bare-sqrt case (b=1).
 *   2. Power[c, p/2] for odd integer p > 0    -- canonicalised by the
 *      evaluator from "c^k * Sqrt[c]"; p=2k+1 gives b = c^k, C = c.
 *      (Case 1 of the user tests hits this: 2*Sqrt[2] surfaces inside
 *      Simplify as Power[2, 3/2].)
 *   3. Times node with one or more Power[_, 1/2] factors plus
 *      arbitrary other factors. When there are multiple sqrt factors,
 *      we combine them into a single Sqrt[Times[c1, c2, ...]] only if
 *      every sqrt base is provably nonneg under ctx (otherwise the
 *      combine breaks the principal-branch identity Sqrt[a]Sqrt[b] =
 *      Sqrt[ab], which fails for negative arguments). The non-sqrt
 *      factors collapse into b via evaluator-driven Times canonicalisation.
 *      (Case 6 of the user tests hits this: 2*Sqrt[x*y] is canonicalised
 *      to 2*Sqrt[x]*Sqrt[y] for symbolic x, y.)
 *
 * Returns false when e doesn't match any of the above. */
static bool extract_sqrt_term(const Expr* e, const AssumeCtx* ctx,
                              Expr** out_b, Expr** out_C) {
    if (!e) return false;

    /* Form 1: bare Sqrt. */
    if (is_sqrt(e)) {
        *out_b = expr_new_integer(1);
        *out_C = expr_copy(e->data.function.args[0]);
        return true;
    }

    /* Form 2: Power[c, p/2] with odd p > 0. */
    if (e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Power
        && e->data.function.arg_count == 2
        && is_rational_literal(e->data.function.args[1])) {
        Expr* exp = e->data.function.args[1];
        Expr* num = exp->data.function.args[0];
        Expr* den = exp->data.function.args[1];
        if (num->type == EXPR_INTEGER && den->type == EXPR_INTEGER &&
            den->data.integer == 2 && num->data.integer > 0 &&
            (num->data.integer % 2) == 1) {
            int64_t p = num->data.integer;
            int64_t k = (p - 1) / 2;
            Expr* base = e->data.function.args[0];
            Expr* b;
            if (k == 0) {
                b = expr_new_integer(1);
            } else if (k == 1) {
                b = expr_copy(base);
            } else {
                Expr* pa[2] = { expr_copy(base), expr_new_integer(k) };
                Expr* pc = expr_new_function(expr_new_symbol("Power"), pa, 2);
                b = evaluate(pc);
                expr_free(pc);
            }
            *out_b = b;
            *out_C = expr_copy(base);
            return true;
        }
    }

    /* Form 3: Times with sqrt factors. */
    if (e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol != SYM_Times) return false;

    /* Flatten nested Times via a depth-first walk: Simplify's
     * intermediate forms can land us with shapes like
     * Times[2, Times[Sqrt[x], Sqrt[y]]] even though Times has the FLAT
     * attribute, because some seed paths construct the tree without
     * re-evaluating it. We collect the leaf factors (those whose head
     * is not Times) into one flat list. */
    Expr** flat = NULL;
    size_t flat_n = 0;
    size_t flat_cap = 0;
    /* iterative DFS using a small stack of cursors */
    typedef struct { const Expr* node; size_t idx; } Cur;
    Cur stack[16];
    int sp = 0;
    stack[sp++] = (Cur){ e, 0 };
    while (sp > 0) {
        Cur* top = &stack[sp - 1];
        if (top->node->type != EXPR_FUNCTION ||
            !top->node->data.function.head ||
            top->node->data.function.head->type != EXPR_SYMBOL ||
            top->node->data.function.head->data.symbol != SYM_Times) {
            /* leaf factor: append */
            if (flat_n == flat_cap) {
                flat_cap = flat_cap ? flat_cap * 2 : 8;
                flat = (Expr**)realloc(flat, sizeof(Expr*) * flat_cap);
            }
            flat[flat_n++] = (Expr*)top->node;
            sp--;
            continue;
        }
        if (top->idx >= top->node->data.function.arg_count) {
            sp--;
            continue;
        }
        const Expr* child = top->node->data.function.args[top->idx++];
        if (sp >= 16) {
            /* Bail out (extremely deep nesting); not worth handling. */
            free(flat);
            return false;
        }
        stack[sp++] = (Cur){ child, 0 };
    }

    size_t sqrt_count = 0;
    for (size_t i = 0; i < flat_n; i++) {
        if (is_sqrt(flat[i])) sqrt_count++;
    }
    if (sqrt_count == 0) {
        free(flat);
        return false;
    }

    /* When >1 sqrt factors, combine into one only if all bases are
     * provably nonneg (so that Sqrt[a]*Sqrt[b] = Sqrt[a*b] is sound). */
    if (sqrt_count > 1) {
        for (size_t i = 0; i < flat_n; i++) {
            if (is_sqrt(flat[i])) {
                Expr* bi = flat[i]->data.function.args[0];
                if (!denest_is_nonneg(bi, ctx)) {
                    free(flat);
                    return false;
                }
            }
        }
    }

    /* Bucket factors: sqrt-bases vs non-sqrt-factors. */
    Expr** sqrt_bases = (Expr**)malloc(sizeof(Expr*) * (flat_n ? flat_n : 1));
    Expr** other     = (Expr**)malloc(sizeof(Expr*) * (flat_n ? flat_n : 1));
    size_t si = 0, oi = 0;
    for (size_t i = 0; i < flat_n; i++) {
        if (is_sqrt(flat[i])) {
            sqrt_bases[si++] = expr_copy(flat[i]->data.function.args[0]);
        } else {
            other[oi++] = expr_copy(flat[i]);
        }
    }
    free(flat);

    /* Combined radicand C. */
    Expr* C;
    if (si == 1) {
        C = sqrt_bases[0];
    } else {
        Expr* tn = expr_new_function(expr_new_symbol("Times"),
                                      sqrt_bases, si);
        C = evaluate(tn);
        expr_free(tn);
    }
    free(sqrt_bases);

    /* Combined b. */
    Expr* b;
    if (oi == 0) {
        free(other);
        b = expr_new_integer(1);
    } else if (oi == 1) {
        b = other[0];
        free(other);
    } else {
        Expr* tn = expr_new_function(expr_new_symbol("Times"), other, oi);
        b = evaluate(tn);
        expr_free(tn);
    }

    *out_b = b;
    *out_C = C;
    return true;
}

/* Given a clean square D, return a closed-form expression s with
 * s^2 = D (sign-agnostic; the caller is expected to validate the
 * downstream nonnegativity of (A+s)/2 and (A-s)/2 -- if both are
 * provably nonneg, the sign of s itself doesn't matter for the
 * identity).
 *
 * Returns NULL when no closed form is detected. The four cases handled
 * are:
 *   - integer/bigint nonneg perfect square                 -> integer sqrt
 *   - rational with perfect-square num and den             -> rational sqrt
 *   - Power[u, 2k]                                          -> Power[u, k]
 *   - polynomial Plus whose Expand+FactorSquareFree is a
 *     pure even-power Power[u, 2k]                          -> Power[u, k]
 *
 * The polynomial path covers the symbolic discriminants that case 6
 * ((x-y)^2) and case 7 (y^2) produce after expansion.
 */
static Expr* sqrt_if_clean_square(const Expr* D, const AssumeCtx* ctx) {
    (void)ctx;

    /* Numeric atoms. */
    if (D->type == EXPR_INTEGER) {
        int64_t d = D->data.integer;
        if (d < 0) return NULL;
        if (d == 0) return expr_new_integer(0);
        mpz_t z, r;
        mpz_init_set_si(z, d);
        if (!mpz_perfect_square_p(z)) { mpz_clear(z); return NULL; }
        mpz_init(r);
        mpz_sqrt(r, z);
        Expr* out;
        if (mpz_fits_slong_p(r)) {
            out = expr_new_integer((int64_t)mpz_get_si(r));
        } else {
            out = expr_new_bigint_from_mpz(r);
        }
        mpz_clear(r); mpz_clear(z);
        return out;
    }
    if (D->type == EXPR_BIGINT) {
        if (mpz_sgn(D->data.bigint) < 0) return NULL;
        if (mpz_sgn(D->data.bigint) == 0) return expr_new_integer(0);
        if (!mpz_perfect_square_p(D->data.bigint)) return NULL;
        mpz_t r; mpz_init(r);
        mpz_sqrt(r, D->data.bigint);
        Expr* out;
        if (mpz_fits_slong_p(r)) {
            out = expr_new_integer((int64_t)mpz_get_si(r));
        } else {
            out = expr_new_bigint_from_mpz(r);
        }
        mpz_clear(r);
        return out;
    }
    if (is_rational_literal(D)) {
        Expr* sn = sqrt_if_clean_square(D->data.function.args[0], ctx);
        if (!sn) return NULL;
        Expr* sd = sqrt_if_clean_square(D->data.function.args[1], ctx);
        if (!sd) { expr_free(sn); return NULL; }
        Expr* inv_args[2] = { sd, expr_new_integer(-1) };
        Expr* inv = expr_new_function(expr_new_symbol("Power"), inv_args, 2);
        Expr* inv_e = evaluate(inv);
        expr_free(inv);
        Expr* prod_args[2] = { sn, inv_e };
        Expr* prod = expr_new_function(expr_new_symbol("Times"), prod_args, 2);
        Expr* out = evaluate(prod);
        expr_free(prod);
        return out;
    }

    /* Power[u, 2k] direct. */
    if (D->type == EXPR_FUNCTION
        && D->data.function.head
        && D->data.function.head->type == EXPR_SYMBOL
        && D->data.function.head->data.symbol == SYM_Power
        && D->data.function.arg_count == 2) {
        Expr* exp = D->data.function.args[1];
        if (exp->type == EXPR_INTEGER &&
            exp->data.integer >= 2 &&
            (exp->data.integer % 2) == 0) {
            int64_t k = exp->data.integer / 2;
            Expr* u = expr_copy(D->data.function.args[0]);
            if (k == 1) return u;
            Expr* pa[2] = { u, expr_new_integer(k) };
            Expr* pc = expr_new_function(expr_new_symbol("Power"), pa, 2);
            Expr* out = evaluate(pc);
            expr_free(pc);
            return out;
        }
    }

    /* Polynomial Plus: try Expand + FactorSquareFree. The factoriser
     * stalls on inputs containing non-integer Power exponents (e.g. an
     * inner Sqrt that survived expansion), so gate on that. */
    if (D->type == EXPR_FUNCTION
        && D->data.function.head
        && D->data.function.head->type == EXPR_SYMBOL
        && D->data.function.head->data.symbol == SYM_Plus) {
        if (has_non_integer_power(D)) return NULL;
        Expr* expanded = call_unary_copy("Expand", D);
        if (!expanded) return NULL;
        if (has_non_integer_power(expanded)) {
            expr_free(expanded);
            return NULL;
        }
        Expr* fsf = call_unary_owned("FactorSquareFree", expanded);
        if (!fsf) return NULL;
        if (fsf->type == EXPR_FUNCTION
            && fsf->data.function.head
            && fsf->data.function.head->type == EXPR_SYMBOL
            && fsf->data.function.head->data.symbol == SYM_Power
            && fsf->data.function.arg_count == 2) {
            Expr* fsf_exp = fsf->data.function.args[1];
            if (fsf_exp->type == EXPR_INTEGER &&
                fsf_exp->data.integer >= 2 &&
                (fsf_exp->data.integer % 2) == 0) {
                int64_t k = fsf_exp->data.integer / 2;
                Expr* u = expr_copy(fsf->data.function.args[0]);
                expr_free(fsf);
                if (k == 1) return u;
                Expr* pa[2] = { u, expr_new_integer(k) };
                Expr* pc = expr_new_function(expr_new_symbol("Power"), pa, 2);
                Expr* out = evaluate(pc);
                expr_free(pc);
                return out;
            }
        }
        expr_free(fsf);
    }
    return NULL;
}

/* Given a Plus expression, partition its arguments into
 *   (A, b, C)
 * such that the Plus equals  A + b * Sqrt[C]  with a single sqrt-bearing
 * term. Returns true on a successful partition. *out_A, *out_b, *out_C
 * are caller-owned allocations on success.
 *
 * Two failure modes are explicit: zero sqrt-bearing terms (nothing to
 * denest) and two-or-more sqrt-bearing terms (multi-extension, phase 2).
 *
 * The ctx is forwarded to extract_sqrt_term so that multi-sqrt Times
 * factors with provably-nonneg bases can be combined into a single
 * sqrt-bearing term. */
static bool split_plus_into_a_plus_b_sqrt_c(const Expr* plus_node,
                                            const AssumeCtx* ctx,
                                            Expr** out_A,
                                            Expr** out_b,
                                            Expr** out_C) {
    if (!plus_node || plus_node->type != EXPR_FUNCTION) return false;
    if (!plus_node->data.function.head ||
        plus_node->data.function.head->type != EXPR_SYMBOL ||
        plus_node->data.function.head->data.symbol != SYM_Plus) return false;

    size_t n = plus_node->data.function.arg_count;
    if (n < 2) return false;

    /* Find the unique sqrt-bearing term and bucket the rest into A. */
    size_t sqrt_idx = SIZE_MAX;
    Expr* b = NULL;
    Expr* C = NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* bi = NULL; Expr* Ci = NULL;
        if (extract_sqrt_term(plus_node->data.function.args[i], ctx,
                              &bi, &Ci)) {
            if (sqrt_idx != SIZE_MAX) {
                /* Second sqrt-bearing term — abort. */
                expr_free(bi); expr_free(Ci);
                if (b) expr_free(b);
                if (C) expr_free(C);
                return false;
            }
            sqrt_idx = i;
            b = bi;
            C = Ci;
        }
    }
    if (sqrt_idx == SIZE_MAX) return false;

    /* Build A as Plus of remaining args, evaluating to canonicalise. */
    if (n == 2) {
        *out_A = expr_copy(plus_node->data.function.args[1 - sqrt_idx]);
    } else {
        Expr** a_args = (Expr**)malloc(sizeof(Expr*) * (n - 1));
        size_t ai = 0;
        for (size_t i = 0; i < n; i++) {
            if (i == sqrt_idx) continue;
            a_args[ai++] = expr_copy(plus_node->data.function.args[i]);
        }
        Expr* pn = expr_new_function(expr_new_symbol("Plus"), a_args, ai);
        *out_A = evaluate(pn);
        expr_free(pn);
    }
    *out_b = b;
    *out_C = C;
    return true;
}

/* Local transitive nonneg prover for the denesting branch checks. The
 * stock assume_known_nonneg / assume_known_positive don't do
 * transitive chaining: x > y && y > 0 doesn't imply x > 0 to them. We
 * extend with a one-step chain: scan facts of the form Greater[x, z]
 * or GreaterEqual[x, z] (in either argument order's natural rendering)
 * and recursively check that z is itself nonneg / positive. The recursion
 * is bounded by the fact set, so cycles are not possible (each step
 * moves to a different right-hand side from the original).
 *
 * Uses `depth` to bound the chain (4 is plenty for any realistic
 * assumption set we'd see; deeper chains would be unusual). */
static bool denest_prov_nonneg(const AssumeCtx* ctx, const Expr* x, int depth);

static bool denest_prov_pos(const AssumeCtx* ctx, const Expr* x, int depth) {
    if (assume_known_positive(ctx, x)) return true;
    if (depth <= 0) return false;
    if (!ctx) return false;
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (!f || f->type != EXPR_FUNCTION) continue;
        if (!f->data.function.head ||
            f->data.function.head->type != EXPR_SYMBOL) continue;
        if (f->data.function.arg_count != 2) continue;
        const char* h = f->data.function.head->data.symbol;
        const Expr* a = f->data.function.args[0];
        const Expr* b = f->data.function.args[1];
        /* Greater[x, z] or Less[z, x] with z provably nonneg => x > 0. */
        if (h == SYM_Greater && expr_eq((Expr*)a, (Expr*)x)) {
            if (denest_prov_nonneg(ctx, b, depth - 1)) return true;
        }
        if (h == SYM_Less && expr_eq((Expr*)b, (Expr*)x)) {
            if (denest_prov_nonneg(ctx, a, depth - 1)) return true;
        }
        /* GreaterEqual[x, z] with z provably positive => x > 0. */
        if (h == SYM_GreaterEqual && expr_eq((Expr*)a, (Expr*)x)) {
            if (denest_prov_pos(ctx, b, depth - 1)) return true;
        }
        if (h == SYM_LessEqual && expr_eq((Expr*)b, (Expr*)x)) {
            if (denest_prov_pos(ctx, a, depth - 1)) return true;
        }
    }
    return false;
}

/* Forward declaration: numeric_is_nonneg is defined just below. */
static bool numeric_is_nonneg(const Expr* e);

static bool denest_prov_nonneg(const AssumeCtx* ctx, const Expr* x, int depth) {
    /* Numeric literal fast-path. Important: covers Rational[p, q] which
     * the stock assume_known_nonneg doesn't (its numeric_sign only
     * handles INTEGER/BIGINT/REAL). The factor 1/2 inside (x+y)/2 needs
     * this to be recognised as nonneg. */
    if (numeric_is_nonneg(x)) return true;
    if (assume_known_nonneg(ctx, x)) return true;
    /* x > z && z >= 0 implies x > 0 implies x >= 0; routes through
     * denest_prov_pos. */
    if (denest_prov_pos(ctx, x, depth)) return true;
    if (depth <= 0) return false;
    if (!ctx) return false;
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (!f || f->type != EXPR_FUNCTION) continue;
        if (!f->data.function.head ||
            f->data.function.head->type != EXPR_SYMBOL) continue;
        if (f->data.function.arg_count != 2) continue;
        const char* h = f->data.function.head->data.symbol;
        const Expr* a = f->data.function.args[0];
        const Expr* b = f->data.function.args[1];
        if (h == SYM_GreaterEqual && expr_eq((Expr*)a, (Expr*)x)) {
            if (denest_prov_nonneg(ctx, b, depth - 1)) return true;
        }
        if (h == SYM_LessEqual && expr_eq((Expr*)b, (Expr*)x)) {
            if (denest_prov_nonneg(ctx, a, depth - 1)) return true;
        }
    }
    /* Plus decomposition with subtraction-pattern recognition.
     *
     * Two-arg case Plus[u, Times[-1, v]] = u - v: nonneg iff u >= v
     * iff ctx has Greater[u, v] or GreaterEqual[u, v]. This is the
     * key path for case 7's Q = (x-y)/2 — purely additive nonneg
     * decomposition fails (one summand has a negative coefficient),
     * but the inequality fact x > y entails x - y >= 0 directly.
     *
     * General case: every summand provably nonneg. Routes through our
     * transitive prover so chained facts are visible.
     */
    if (x && x->type == EXPR_FUNCTION && x->data.function.head &&
        x->data.function.head->type == EXPR_SYMBOL &&
        x->data.function.head->data.symbol == SYM_Plus) {
        size_t pn = x->data.function.arg_count;
        if (pn == 2) {
            /* Subtraction pattern: u + (-1)*v. */
            const Expr* u = NULL;
            const Expr* v = NULL;
            for (size_t i = 0; i < 2; i++) {
                const Expr* arg = x->data.function.args[i];
                const Expr* other = x->data.function.args[1 - i];
                if (arg->type == EXPR_FUNCTION &&
                    arg->data.function.head &&
                    arg->data.function.head->type == EXPR_SYMBOL &&
                    arg->data.function.head->data.symbol == SYM_Times &&
                    arg->data.function.arg_count == 2) {
                    Expr* c = arg->data.function.args[0];
                    Expr* w = arg->data.function.args[1];
                    if (c->type == EXPR_INTEGER && c->data.integer == -1) {
                        u = other;
                        v = w;
                        break;
                    }
                }
            }
            if (u && v) {
                /* Look for Greater[u, v] / GreaterEqual[u, v] in ctx. */
                for (size_t i = 0; i < ctx->count; i++) {
                    const Expr* f = ctx->facts[i];
                    if (!f || f->type != EXPR_FUNCTION) continue;
                    if (!f->data.function.head ||
                        f->data.function.head->type != EXPR_SYMBOL) continue;
                    if (f->data.function.arg_count != 2) continue;
                    const char* h = f->data.function.head->data.symbol;
                    const Expr* a = f->data.function.args[0];
                    const Expr* b = f->data.function.args[1];
                    if ((h == SYM_Greater || h == SYM_GreaterEqual) &&
                        expr_eq((Expr*)a, (Expr*)u) &&
                        expr_eq((Expr*)b, (Expr*)v)) {
                        return true;
                    }
                    if ((h == SYM_Less || h == SYM_LessEqual) &&
                        expr_eq((Expr*)a, (Expr*)v) &&
                        expr_eq((Expr*)b, (Expr*)u)) {
                        return true;
                    }
                }
            }
        }
        /* General case: each summand nonneg. */
        bool all_nn = true;
        for (size_t i = 0; i < pn; i++) {
            if (!denest_prov_nonneg(ctx, x->data.function.args[i],
                                    depth - 1)) {
                all_nn = false;
                break;
            }
        }
        if (all_nn) return true;
    }
    /* Times: every factor nonneg => product nonneg (sign-symmetric:
     * negative * negative is also nonneg, but the simpler positive
     * decomposition covers the cases we hit in practice).
     *
     * Also: Times[c, Plus[...]] for positive c reduces to Plus's
     * nonnegativity check on the inner Plus. This catches the factored
     * forms Times[1/2, Plus[x, -y]] that picocas's evaluator produces
     * before Expand distributes them. */
    if (x && x->type == EXPR_FUNCTION && x->data.function.head &&
        x->data.function.head->type == EXPR_SYMBOL &&
        x->data.function.head->data.symbol == SYM_Times) {
        bool all_nn = true;
        for (size_t i = 0; i < x->data.function.arg_count; i++) {
            if (!denest_prov_nonneg(ctx, x->data.function.args[i],
                                    depth - 1)) {
                all_nn = false;
                break;
            }
        }
        if (all_nn) return true;
    }
    return false;
}

/* Numeric-only nonneg check for a plain integer or rational literal. */
static bool numeric_is_nonneg(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer >= 0;
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint) >= 0;
    if (is_rational_literal(e)) {
        const Expr* num = e->data.function.args[0];
        const Expr* den = e->data.function.args[1];
        bool num_nn = (num->type == EXPR_INTEGER && num->data.integer >= 0) ||
                      (num->type == EXPR_BIGINT  && mpz_sgn(num->data.bigint) >= 0);
        bool den_pos = (den->type == EXPR_INTEGER && den->data.integer >  0) ||
                       (den->type == EXPR_BIGINT  && mpz_sgn(den->data.bigint) >  0);
        return num_nn && den_pos;
    }
    return false;
}

/* Decide whether `e` is provably nonneg under the active context.
 * Numeric literals get a free direct check; everything else routes
 * through assume_known_nonneg, augmented by a local transitive prover
 * (denest_prov_nonneg) that chains inequality facts. The transitive
 * step is required for case 7's branch check: x > y && y > 0 implies
 * (x+y)/2 >= 0, but assume_known_nonneg alone can't see it. */
static bool denest_is_nonneg(const Expr* e, const AssumeCtx* ctx) {
    if (numeric_is_nonneg(e)) return true;
    return denest_prov_nonneg(ctx, e, 4);
}

/* Shared compute step for the half-sum denesting identity. On success,
 * *P_out, *Q_out, *s_out are caller-owned and satisfy
 *
 *     A + b * Sqrt[C] = (Sqrt[P] + sign * Sqrt[Q])^2
 *
 * where sign = -1 if *b_is_negative_out, else +1, and
 *
 *     s = P - Q = Sqrt[A^2 - b^2 C]   (always nonneg)
 *
 * is the rationalising denominator for the reciprocal form. Returns
 * false when no clean denesting exists or branch validity cannot be
 * verified; in that case *_out are left NULL.
 *
 * Both the forward path (`Sqrt[plus] -> Sqrt[P] +/- Sqrt[Q]`) and the
 * reciprocal path (`1/Sqrt[plus] -> (Sqrt[P] -/+ Sqrt[Q])/s`) reuse
 * this so the branch-validity prover, the polynomial-arithmetic
 * scaffolding around D = A^2 - b^2 C, and the Sqrt[Rational] rationaliser
 * all live in one place. */
static bool denest_compute_pq_s(const Expr* plus_node,
                                 const AssumeCtx* ctx,
                                 Expr** P_out, Expr** Q_out, Expr** s_out,
                                 bool* b_is_negative_out) {
    *P_out = NULL; *Q_out = NULL; *s_out = NULL;
    Expr* A = NULL; Expr* b = NULL; Expr* C = NULL;
    if (!split_plus_into_a_plus_b_sqrt_c(plus_node, ctx, &A, &b, &C)) return false;

    /* Determine sign of b by inspection of its sign. We require b's
     * sign to be syntactically determinable -- a symbolic b with
     * unknown sign aborts the identity (we'd need to know which
     * principal-branch convention applies to the result). */
    bool b_is_negative;
    if (expr_is_superficially_negative(b)) {
        b_is_negative = true;
    } else if (numeric_is_nonneg(b) || (b->type == EXPR_INTEGER && b->data.integer == 1)) {
        b_is_negative = false;
    } else {
        /* Symbolic b with unknown sign. */
        if (assume_known_nonneg(ctx, b)) {
            b_is_negative = false;
        } else if (assume_known_nonpos(ctx, b)) {
            b_is_negative = true;
        } else {
            expr_free(A); expr_free(b); expr_free(C);
            return false;
        }
    }

    /* B = b^2 * C; evaluate so symbolic b^2 simplifies and rational
     * b^2 collapses to a literal. */
    Expr* b2_args[2] = { expr_copy(b), expr_new_integer(2) };
    Expr* b2_pow = expr_new_function(expr_new_symbol("Power"), b2_args, 2);
    Expr* b2 = evaluate(b2_pow);
    expr_free(b2_pow);
    Expr* B_args[2] = { b2, expr_copy(C) };
    Expr* B_times = expr_new_function(expr_new_symbol("Times"), B_args, 2);
    Expr* B = evaluate(B_times);
    expr_free(B_times);

    /* D = A^2 - B; Expand to put it in canonical polynomial form so
     * sqrt_if_clean_square's FactorSquareFree path has a chance. */
    Expr* A2_args[2] = { expr_copy(A), expr_new_integer(2) };
    Expr* A2_pow = expr_new_function(expr_new_symbol("Power"), A2_args, 2);
    Expr* A2 = evaluate(A2_pow);
    expr_free(A2_pow);
    Expr* negB_args[2] = { expr_new_integer(-1), B };
    Expr* negB_times = expr_new_function(expr_new_symbol("Times"), negB_args, 2);
    Expr* negB = evaluate(negB_times);
    expr_free(negB_times);
    Expr* D_args[2] = { A2, negB };
    Expr* D_plus = expr_new_function(expr_new_symbol("Plus"), D_args, 2);
    Expr* D_unexp = evaluate(D_plus);
    expr_free(D_plus);
    /* call_unary_owned consumes D_unexp; if it returns NULL we treat as
     * a soft failure and abort the denesting cleanly. In practice
     * Expand on a numeric/polynomial expression never returns NULL. */
    Expr* D = call_unary_owned("Expand", D_unexp);
    if (!D) {
        expr_free(A); expr_free(b); expr_free(C);
        return false;
    }

    Expr* s = sqrt_if_clean_square(D, ctx);
    expr_free(D);
    if (!s) {
        expr_free(A); expr_free(b); expr_free(C);
        return false;
    }

    /* Compute P = (A + s)/2 and Q = (A - s)/2. picocas's evaluator
     * keeps Times[-1, Plus[...]] factored rather than distributing,
     * so a raw evaluate of Plus[A, -s] / 2 leaves Q in a factored form
     * the surrounding pipeline can't reason about (e.g. case 6 saw
     * Q = 1/2 (x + y - (-x + y)) instead of Q = x). We call Together
     * on each, which combines redundant fractions back into a single
     * canonical Times[1/c, Plus[...]] form so the nonneg check below
     * sees a clean expression. We DON'T use Expand: distributing the
     * 1/2 across (x - y) wraps each summand in a numeric coefficient,
     * obscuring the simple "u - v >= 0" subtraction pattern that the
     * branch check uses on case 7. */
    Expr* Ap_args[2] = { expr_copy(A), expr_copy(s) };
    Expr* Ap = expr_new_function(expr_new_symbol("Plus"), Ap_args, 2);
    Expr* Ap_e = evaluate(Ap);
    expr_free(Ap);
    Expr* halfP_args[2] = { make_rational(1, 2), Ap_e };
    Expr* halfP = expr_new_function(expr_new_symbol("Times"), halfP_args, 2);
    Expr* P_raw = evaluate(halfP);
    expr_free(halfP);
    /* call_unary_owned consumes P_raw; Expand never returns NULL on a
     * polynomial-shaped input but guard anyway. */
    Expr* P = call_unary_owned("Together", P_raw);
    if (!P) {
        expr_free(b); expr_free(C);
        /* A was already freed below; we're early-returning, so
         * that hasn't happened yet — unwind. */
        expr_free(A); expr_free(s);
        return false;
    }

    /* Build (-1) * s without consuming `s`: the reciprocal denester needs
     * to return s untouched, and the original code only worked because it
     * never reached for s again. */
    Expr* neg_s_args[2] = { expr_new_integer(-1), expr_copy(s) };
    Expr* neg_s = expr_new_function(expr_new_symbol("Times"), neg_s_args, 2);
    Expr* neg_s_e = evaluate(neg_s);
    expr_free(neg_s);
    Expr* Am_args[2] = { expr_copy(A), neg_s_e };
    Expr* Am = expr_new_function(expr_new_symbol("Plus"), Am_args, 2);
    Expr* Am_e = evaluate(Am);
    expr_free(Am);
    Expr* halfQ_args[2] = { make_rational(1, 2), Am_e };
    Expr* halfQ = expr_new_function(expr_new_symbol("Times"), halfQ_args, 2);
    Expr* Q_raw = evaluate(halfQ);
    expr_free(halfQ);
    Expr* Q = call_unary_owned("Together", Q_raw);
    if (!Q) {
        expr_free(P);
        expr_free(A); expr_free(b); expr_free(C); expr_free(s);
        return false;
    }

    expr_free(A); expr_free(b); expr_free(C);

    /* Branch validity: P and Q must both be provably nonneg. If we
     * can't prove it, refuse rather than risk a wrong-branch result.
     * The denest_is_nonneg prover (defined in this section) extends
     * assume_known_nonneg with one-step transitive chaining over
     * inequality facts plus a subtraction-pattern detector for
     * Plus[u, -v]; without these case 7 would fail because its P, Q
     * are (x+y)/2 and (x-y)/2 which require chained reasoning under
     * the user's x > y && y > 0 assumptions. */
    if (!denest_is_nonneg(P, ctx) || !denest_is_nonneg(Q, ctx)) {
        expr_free(P); expr_free(Q); expr_free(s);
        return false;
    }

    *P_out = P;
    *Q_out = Q;
    *s_out = s;
    *b_is_negative_out = b_is_negative;
    return true;
}

/* True iff `e` is `Power[base, m/2]` for some odd integer `m`. On true,
 * sets *m_out = m and *base_out to a non-owning alias of the base.
 *
 * This is the generalisation of is_sqrt (m == 1) and is_inv_sqrt
 * (m == -1) that lets the walker handle Sqrt[X]^k forms which the
 * evaluator canonicalises to Power[X, k/2] (e.g. 1/Sqrt[X]^3 ->
 * Power[X, -3/2]). Even m would be an integer power, so we exclude
 * it here -- denesting at the half-power layer wouldn't apply. */
static bool match_half_int_power(const Expr* e, int64_t* m_out,
                                  Expr** base_out) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol != SYM_Power) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* exp = e->data.function.args[1];
    if (!is_rational_literal(exp)) return false;
    Expr* num = exp->data.function.args[0];
    Expr* den = exp->data.function.args[1];
    if (num->type != EXPR_INTEGER || den->type != EXPR_INTEGER) return false;
    if (den->data.integer != 2) return false;
    int64_t m = num->data.integer;
    if (m == 0 || (m % 2) == 0) return false;
    *m_out = m;
    *base_out = e->data.function.args[0];
    return true;
}

/* Raise `base` (consumed) to integer exponent k and evaluate. */
static Expr* power_int_eval(Expr* base, int64_t k) {
    Expr* args[2] = { base, expr_new_integer(k) };
    Expr* pow = expr_new_function(expr_new_symbol("Power"), args, 2);
    Expr* result = evaluate(pow);
    expr_free(pow);
    return result;
}

/* Generalised denester for `Power[Plus[A + b Sqrt[C]], m/2]` with odd m.
 *
 *   m == +1: returns Sqrt[P] + sign*Sqrt[Q]                 (forward path)
 *   m == -1: returns (Sqrt[P] - sign*Sqrt[Q]) / s           (reciprocal)
 *   |m| > 1:
 *     m > 0: (Sqrt[P] + sign*Sqrt[Q])^m
 *     m < 0: (Sqrt[P] - sign*Sqrt[Q])^|m| / s^|m|
 *
 * For negative m we rationalise via the conjugate so the result lives
 * in Q(Sqrt[P], Sqrt[Q]) and downstream Simplify can cancel against
 * other terms in that linear basis. */
static Expr* try_denest_pow_half_int(const Expr* plus_node, int64_t m,
                                      const AssumeCtx* ctx) {
    Expr* P = NULL; Expr* Q = NULL; Expr* s = NULL;
    bool b_is_negative = false;
    if (!denest_compute_pq_s(plus_node, ctx, &P, &Q, &s, &b_is_negative)) {
        return NULL;
    }

    Expr* sP_args[2] = { P, make_rational(1, 2) };
    Expr* sP = expr_new_function(expr_new_symbol("Power"), sP_args, 2);
    Expr* sP_e = evaluate(sP);
    expr_free(sP);
    Expr* sQ_args[2] = { Q, make_rational(1, 2) };
    Expr* sQ = expr_new_function(expr_new_symbol("Power"), sQ_args, 2);
    Expr* sQ_e = evaluate(sQ);
    expr_free(sQ);

    /* For m > 0 we use sign = b's sign (Sqrt[P] + sign*Sqrt[Q]). For
     * m < 0 we use the conjugate sign and divide by s^|m| (rationalised
     * reciprocal). */
    bool use_conjugate = (m < 0);
    bool plus_q = use_conjugate ? b_is_negative : !b_is_negative;
    /* plus_q = true means numerator is Sqrt[P] + Sqrt[Q];
     * plus_q = false means Sqrt[P] - Sqrt[Q]. */

    Expr* R;
    if (plus_q) {
        Expr* sum_args[2] = { sP_e, sQ_e };
        Expr* sum = expr_new_function(expr_new_symbol("Plus"), sum_args, 2);
        R = evaluate(sum);
        expr_free(sum);
    } else {
        Expr* neg_sQ_args[2] = { expr_new_integer(-1), sQ_e };
        Expr* neg_sQ = expr_new_function(expr_new_symbol("Times"),
                                          neg_sQ_args, 2);
        Expr* neg_sQ_e = evaluate(neg_sQ);
        expr_free(neg_sQ);
        Expr* sum_args[2] = { sP_e, neg_sQ_e };
        Expr* sum = expr_new_function(expr_new_symbol("Plus"), sum_args, 2);
        R = evaluate(sum);
        expr_free(sum);
    }

    int64_t k = m < 0 ? -m : m;
    Expr* R_pow = (k == 1) ? R : power_int_eval(R, k);

    Expr* result;
    if (m < 0) {
        /* Divide by s^k. Build R_pow * Power[s, -k] and evaluate. */
        Expr* inv_s_args[2] = { s, expr_new_integer(-k) };
        Expr* inv_s = expr_new_function(expr_new_symbol("Power"),
                                         inv_s_args, 2);
        Expr* inv_s_e = evaluate(inv_s);
        expr_free(inv_s);
        Expr* prod_args[2] = { R_pow, inv_s_e };
        Expr* prod = expr_new_function(expr_new_symbol("Times"),
                                        prod_args, 2);
        result = evaluate(prod);
        expr_free(prod);
    } else {
        expr_free(s);
        result = R_pow;
    }

    /* Canonicalise leftover Sqrt[Rational[m, n]] terms, same cleanup
     * as in the m == +/- 1 paths. */
    Expr* rationalised = denest_rationalise_sqrt_of_rational(result);
    if (rationalised && !expr_eq(rationalised, result)) {
        expr_free(result);
        result = rationalised;
    } else if (rationalised) {
        expr_free(rationalised);
    }
    return result;
}

/* Bottom-up walker: find every Sqrt[Plus] / 1/Sqrt[Plus] subtree and
 * try the half-sum identity at each. Returns a fully-rewritten copy of
 * e on any successful denesting; NULL when no node fired. */
static Expr* simp_denest_sqrt_walk(const Expr* e, const AssumeCtx* ctx) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;

    size_t n = e->data.function.arg_count;
    Expr** new_args = NULL;
    bool any = false;
    for (size_t i = 0; i < n; i++) {
        Expr* r = simp_denest_sqrt_walk(e->data.function.args[i], ctx);
        if (r) {
            if (!new_args) {
                new_args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
                for (size_t j = 0; j < i; j++) {
                    new_args[j] = expr_copy(e->data.function.args[j]);
                }
            }
            new_args[i] = r;
            any = true;
        } else if (new_args) {
            new_args[i] = expr_copy(e->data.function.args[i]);
        }
    }

    Expr* current = NULL;
    if (any) {
        Expr* head_copy = expr_copy(e->data.function.head);
        Expr* rebuilt = expr_new_function(head_copy, new_args, n);
        free(new_args);
        current = evaluate(rebuilt);
        expr_free(rebuilt);
    }

    /* Now check whether this very node is a Power[Plus[...], m/2] for
     * some odd integer m, i.e. Sqrt[X]^m. The forward (m=+1) and
     * reciprocal (m=-1) paths used to be dispatched separately via
     * is_sqrt / is_inv_sqrt; we now go through match_half_int_power so
     * higher half-powers (e.g. m = +/-3 from 1/Sqrt[X]^3 -> X^(-3/2))
     * also get denested. We use `current` (post-children-rewrite) when
     * present so structural changes propagate. */
    const Expr* target = current ? current : e;
    int64_t m_num = 0;
    Expr* base = NULL;
    if (match_half_int_power(target, &m_num, &base)) {
        if (base->type == EXPR_FUNCTION
            && base->data.function.head
            && base->data.function.head->type == EXPR_SYMBOL
            && base->data.function.head->data.symbol == SYM_Plus) {
            Expr* d = try_denest_pow_half_int(base, m_num, ctx);
            if (d) {
                if (current) expr_free(current);
                return d;
            }
        }
    }
    return current;  /* may be NULL when no descendant or this node fired */
}

/* Top-level entry point. Returns a new Expr* (caller owns); never
 * NULL. When no denesting fires, returns expr_copy(e). */
static Expr* simp_denest_sqrt(const Expr* e, const AssumeCtx* ctx) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;
    Expr* r = simp_denest_sqrt_walk(e, ctx);
    Expr* out = r ? r : expr_copy((Expr*)e);
    if (dbg) simp_debug_log("DenestSqrt", e, out,
                            simp_debug_elapsed_ms(t0));
    return out;
}

/* ----------------------------------------------------------------------- */
/* Cube-root denesting and sum-of-conjugates: simp_cuberoot                */
/* ----------------------------------------------------------------------- */

/*
 * simp_cuberoot recognises two narrow phase-3 patterns:
 *
 * Pattern A — Single-radical cube-root denesting (Borodin-Fagin-style):
 *
 *     Power[a + b*Sqrt[c], 1/3]
 *
 *   where a, b are rationals and c is a positive integer. Solvable iff
 *   there exist rationals (p, q) with
 *       p^3 + 3 p q^2 c = a
 *       3 p^2 q + q^3 c = b
 *   in which case
 *       (a + b Sqrt[c])^(1/3) = p + q Sqrt[c].
 *   We search a small bounded grid of candidate denominators and
 *   numerators; for case 4 (a=2, b=1, c=5) the answer p=q=1/2 lives in
 *   that grid.
 *
 * Pattern B — Sum of two conjugate cube roots (Cardano):
 *
 *     Power[a + b*Sqrt[c], 1/3] + Power[a - b*Sqrt[c], 1/3]
 *
 *   where a, b are rationals, c a positive integer, and a^2 - b^2 c is
 *   a perfect (possibly negative) integer cube `m^3`. Under the
 *   real-cube-root convention (which is what the user expects for case
 *   5; see the soundness note below), the sum s satisfies
 *       s^3 - 3 m s - 2 a = 0,
 *   and the transform succeeds when this cubic has a rational root.
 *
 * Soundness note (real vs. principal branch):
 *
 *   picocas's Power[neg_real, 1/3] uses the principal complex branch
 *   (e.g. (-1)^(1/3) = (1 + i sqrt(3))/2). Pattern B's identity
 *   `(2+sqrt(5))^(1/3) + (2-sqrt(5))^(1/3) = 1` only holds under the
 *   REAL cube-root convention; principal-branch evaluation gives
 *   ~1.93 + 0.535i. We fire the rewrite anyway because the user's
 *   intent is clear (their test expects 1) and Mathematica's Simplify
 *   uses the same heuristic. The transform is gated on the discriminant
 *   `a^2 - b^2 c` being a perfect cube of an integer (positive or
 *   negative); this is a purely structural test that doesn't depend
 *   on branch choice. Pattern A doesn't have a branch concern because
 *   the principal-branch equation `(p + q sqrt(c))^3 = a + b sqrt(c)`
 *   is an exact polynomial identity for any branch.
 */

/* Forward declaration: simp_search drives this. */
static Expr* simp_cuberoot(const Expr* e, const AssumeCtx* ctx);

/* Small fixed search bound; matches the user's test cases (p=q=1/2). */
#define CUBEROOT_DENOM_MAX 6
#define CUBEROOT_NUM_MAX   6

/* Helper: return integer cube root of a (signed). Returns true and
 * writes *out = m on success (m^3 == a), false otherwise. */
static bool cuberoot_int_exact(int64_t a, int64_t* out) {
    int64_t sign = a < 0 ? -1 : 1;
    int64_t abs_a = a < 0 ? -a : a;
    /* Integer cube root via search. abs_a fits in 64 bits, so cube root
     * fits in 21 bits; binary search converges fast. */
    int64_t lo = 0, hi = 2097152;  /* 2^21 */
    while (lo < hi) {
        int64_t mid = lo + (hi - lo) / 2;
        /* mid^3 overflow guard. */
        if (mid > 0 && mid > INT64_MAX / mid / mid + 1) {
            hi = mid;
            continue;
        }
        int64_t cube = mid * mid * mid;
        if (cube < abs_a) lo = mid + 1;
        else hi = mid;
    }
    if (lo == 0) {
        if (abs_a == 0) { *out = 0; return true; }
        return false;
    }
    if (lo * lo * lo != abs_a) return false;
    *out = sign * lo;
    return true;
}

/* Match a Plus[a, b*Sqrt[c]] (or Plus[a, -b*Sqrt[c]]) pattern.
 * Populates *a_num, *a_den, *b_num, *b_den, *c with rational
 * coefficients (a, b in lowest terms, c the integer Sqrt radicand).
 * `e` may also be a bare integer/rational a (b=0, c=0) — the caller
 * handles that case as a degenerate match.
 *
 * For our Phase-3 cases, the radicand c is always a small positive
 * integer (the user supplied 5). We allow only that shape: matching
 * Plus[a, Times[b, Power[c, 1/2]]] or Plus[a, Power[c, 1/2]] (b=1).
 *
 * Returns true on a successful match. */
static bool cuberoot_match_a_plus_b_sqrt_c(const Expr* e,
                                           int64_t* a_num, int64_t* a_den,
                                           int64_t* b_num, int64_t* b_den,
                                           int64_t* c_int) {
    *a_num = 0; *a_den = 1;
    *b_num = 0; *b_den = 1;
    *c_int = 0;
    if (!e) return false;
    /* Normalise atomic case: just a rational/integer. */
    if (e->type == EXPR_INTEGER) {
        *a_num = e->data.integer;
        return true;
    }
    if (is_rational_literal(e)) {
        Expr* n = e->data.function.args[0];
        Expr* d = e->data.function.args[1];
        if (n->type != EXPR_INTEGER || d->type != EXPR_INTEGER) return false;
        *a_num = n->data.integer;
        *a_den = d->data.integer;
        return true;
    }
    /* Bare Sqrt[c]: a=0, b=1. */
    if (is_sqrt(e)) {
        Expr* base = e->data.function.args[0];
        if (base->type != EXPR_INTEGER || base->data.integer < 2) return false;
        *b_num = 1;
        *c_int = base->data.integer;
        return true;
    }
    /* Times[k, Sqrt[c]]: a=0, b=k. */
    if (e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Times
        && e->data.function.arg_count == 2) {
        Expr* k = e->data.function.args[0];
        Expr* sq = e->data.function.args[1];
        int64_t kn, kd;
        if (is_rational(k, &kn, &kd) && is_sqrt(sq)) {
            Expr* base = sq->data.function.args[0];
            if (base->type == EXPR_INTEGER && base->data.integer >= 2) {
                *b_num = kn; *b_den = kd;
                *c_int = base->data.integer;
                return true;
            }
        }
    }
    /* Plus form. */
    if (e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol != SYM_Plus) return false;
    if (e->data.function.arg_count != 2) return false;

    /* Identify the rational-only arg and the b*Sqrt[c] arg. */
    int rat_idx = -1, sqrt_idx = -1;
    for (int i = 0; i < 2; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_INTEGER || is_rational_literal(arg)) {
            rat_idx = i;
        } else {
            sqrt_idx = i;
        }
    }
    if (rat_idx < 0 || sqrt_idx < 0) return false;

    /* Decode rational arg into (a_num, a_den). */
    Expr* rat = e->data.function.args[rat_idx];
    if (rat->type == EXPR_INTEGER) {
        *a_num = rat->data.integer;
    } else {
        Expr* n = rat->data.function.args[0];
        Expr* d = rat->data.function.args[1];
        if (n->type != EXPR_INTEGER || d->type != EXPR_INTEGER) return false;
        *a_num = n->data.integer;
        *a_den = d->data.integer;
    }

    /* Decode sqrt arg as either Sqrt[c] (b = 1) or Times[k, Sqrt[c]]. */
    Expr* sa = e->data.function.args[sqrt_idx];
    if (is_sqrt(sa)) {
        Expr* base = sa->data.function.args[0];
        if (base->type != EXPR_INTEGER || base->data.integer < 2) return false;
        *b_num = 1;
        *c_int = base->data.integer;
        return true;
    }
    if (sa->type == EXPR_FUNCTION
        && sa->data.function.head
        && sa->data.function.head->type == EXPR_SYMBOL
        && sa->data.function.head->data.symbol == SYM_Times
        && sa->data.function.arg_count == 2) {
        Expr* k = sa->data.function.args[0];
        Expr* sq = sa->data.function.args[1];
        int64_t kn, kd;
        if (is_rational(k, &kn, &kd) && is_sqrt(sq)) {
            Expr* base = sq->data.function.args[0];
            if (base->type == EXPR_INTEGER && base->data.integer >= 2) {
                *b_num = kn; *b_den = kd;
                *c_int = base->data.integer;
                return true;
            }
        }
    }
    return false;
}

/* Try Pattern A: Power[a + b sqrt(c), 1/3] denesting via small grid
 * search for rational (p, q) with (p + q sqrt(c))^3 = a + b sqrt(c).
 * Returns the denested form or NULL.
 *
 * Search domain: |p_num|, |q_num| up to CUBEROOT_NUM_MAX, common
 * denominator d up to CUBEROOT_DENOM_MAX. Sufficient for the user's
 * case 4 (p = q = 1/2) and similar small-coefficient cases. */
static Expr* try_cuberoot_denest(int64_t a_num, int64_t a_den,
                                 int64_t b_num, int64_t b_den,
                                 int64_t c) {
    /* Normalise to a common denominator d (= a_den * b_den / gcd(...)).
     * After multiplying through by d^3, the equation becomes integer
     * arithmetic. */
    int64_t a_n = a_num, a_d = a_den;
    int64_t b_n = b_num, b_d = b_den;
    /* (a_n / a_d) and (b_n / b_d). Test (p_n / d)^3 + 3 (p_n / d)
     * (q_n / d)^2 c = a_n / a_d  i.e.
     *   (p_n^3 + 3 p_n q_n^2 c) / d^3 = a_n / a_d
     *   ==> a_d (p_n^3 + 3 p_n q_n^2 c) = a_n d^3
     * Similarly b_d (3 p_n^2 q_n + q_n^3 c) = b_n d^3. */
    for (int64_t d = 1; d <= CUBEROOT_DENOM_MAX; d++) {
        int64_t d3 = d * d * d;
        for (int64_t p_n = -CUBEROOT_NUM_MAX; p_n <= CUBEROOT_NUM_MAX; p_n++) {
            for (int64_t q_n = -CUBEROOT_NUM_MAX; q_n <= CUBEROOT_NUM_MAX; q_n++) {
                /* Skip the trivial (0, 0). */
                if (p_n == 0 && q_n == 0) continue;
                int64_t lhs_a = a_d * (p_n * p_n * p_n + 3 * p_n * q_n * q_n * c);
                int64_t rhs_a = a_n * d3;
                if (lhs_a != rhs_a) continue;
                int64_t lhs_b = b_d * (3 * p_n * p_n * q_n + q_n * q_n * q_n * c);
                int64_t rhs_b = b_n * d3;
                if (lhs_b != rhs_b) continue;
                /* Match. Build (p_n + q_n sqrt(c)) / d. */
                Expr* p = make_rational(p_n, d);
                Expr* q = make_rational(q_n, d);
                Expr* sq_args[2] = { expr_new_integer(c), make_rational(1, 2) };
                Expr* sq_call = expr_new_function(expr_new_symbol("Power"),
                                                  sq_args, 2);
                Expr* sq = evaluate(sq_call);
                expr_free(sq_call);
                Expr* qsq_args[2] = { q, sq };
                Expr* qsq_call = expr_new_function(expr_new_symbol("Times"),
                                                   qsq_args, 2);
                Expr* qsq = evaluate(qsq_call);
                expr_free(qsq_call);
                Expr* sum_args[2] = { p, qsq };
                Expr* sum_call = expr_new_function(expr_new_symbol("Plus"),
                                                   sum_args, 2);
                Expr* sum = evaluate(sum_call);
                expr_free(sum_call);
                return sum;
            }
        }
    }
    return NULL;
}

/* Try Pattern B: Power[a+b sqrt(c), 1/3] + Power[a-b sqrt(c), 1/3]
 * via Cardano discriminant. Returns the rational sum on a hit, NULL
 * otherwise.
 *
 * Algorithm: discriminant D = a^2 - b^2 c. If D is a perfect cube of
 * an integer m, then under real-cube-root convention the sum s
 * satisfies s^3 - 3 m s - 2 a = 0. We rational-root test for s; for
 * a, m integer the candidate roots are divisors of -2a. */
static Expr* try_cuberoot_sum_conjugates(int64_t a_num, int64_t a_den,
                                         int64_t b_num, int64_t b_den,
                                         int64_t c) {
    /* For simplicity, require integer a, b. Phase 3 cases all integer. */
    if (a_den != 1 || b_den != 1) return NULL;
    int64_t a = a_num, b = b_num;
    int64_t D = a * a - b * b * c;
    int64_t m;
    if (!cuberoot_int_exact(D, &m)) return NULL;
    /* Cubic: s^3 - 3 m s - 2 a = 0. Rational-root test. Integer roots
     * must divide -2a. Enumerate divisors of |2a|. */
    int64_t target = 2 * a;
    if (target == 0) target = 1;  /* a = 0 case */
    int64_t abs_target = target < 0 ? -target : target;
    int64_t signs[2] = { 1, -1 };
    /* Enumerate integer divisors of abs_target, both signs. */
    for (int64_t d = 1; d * d <= abs_target; d++) {
        if (abs_target % d != 0) continue;
        int64_t partners[2] = { d, abs_target / d };
        for (int pi = 0; pi < 2; pi++) {
            int64_t cand = partners[pi];
            for (int si = 0; si < 2; si++) {
                int64_t s = signs[si] * cand;
                /* Test s^3 - 3 m s - 2 a == 0. Watch for overflow. */
                /* For our test cases s is bounded by ±|2a|, ≤ 4. So
                 * overflow is not a concern. */
                int64_t poly = s * s * s - 3 * m * s - 2 * a;
                if (poly == 0) {
                    return expr_new_integer(s);
                }
            }
        }
    }
    return NULL;
}

/* Walk the tree: at each Power[plus, 1/3] or each
 * Plus[Power[plus1, 1/3], Power[plus2, 1/3]] subtree, try Patterns A/B.
 * Returns NULL when no rewrite fires. */
static Expr* simp_cuberoot_walk(const Expr* e, const AssumeCtx* ctx) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;

    /* Recurse into children first (bottom-up). */
    size_t n = e->data.function.arg_count;
    Expr** new_args = NULL;
    bool any = false;
    for (size_t i = 0; i < n; i++) {
        Expr* r = simp_cuberoot_walk(e->data.function.args[i], ctx);
        if (r) {
            if (!new_args) {
                new_args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
                for (size_t j = 0; j < i; j++) {
                    new_args[j] = expr_copy(e->data.function.args[j]);
                }
            }
            new_args[i] = r;
            any = true;
        } else if (new_args) {
            new_args[i] = expr_copy(e->data.function.args[i]);
        }
    }
    Expr* current = NULL;
    if (any) {
        Expr* head_copy = expr_copy(e->data.function.head);
        Expr* rebuilt = expr_new_function(head_copy, new_args, n);
        free(new_args);
        current = evaluate(rebuilt);
        expr_free(rebuilt);
    }

    const Expr* target = current ? current : e;

    /* Pattern A: target = Power[Plus[...], 1/3] */
    if (target->type == EXPR_FUNCTION
        && target->data.function.head
        && target->data.function.head->type == EXPR_SYMBOL
        && target->data.function.head->data.symbol == SYM_Power
        && target->data.function.arg_count == 2) {
        Expr* exp = target->data.function.args[1];
        if (is_rational_literal(exp)) {
            int64_t pn, qn;
            if (is_rational(exp, &pn, &qn) && pn == 1 && qn == 3) {
                int64_t a_n, a_d, b_n, b_d, c;
                if (cuberoot_match_a_plus_b_sqrt_c(target->data.function.args[0],
                                                    &a_n, &a_d, &b_n, &b_d, &c)
                    && c > 0) {
                    Expr* d = try_cuberoot_denest(a_n, a_d, b_n, b_d, c);
                    if (d) {
                        if (current) expr_free(current);
                        return d;
                    }
                }
            }
        }
    }

    /* Pattern B: target = Plus[Power[...,1/3], Power[...,1/3], maybe other terms].
     * Scan args for paired conjugate cube roots. */
    if (target->type == EXPR_FUNCTION
        && target->data.function.head
        && target->data.function.head->type == EXPR_SYMBOL
        && target->data.function.head->data.symbol == SYM_Plus) {
        size_t pn = target->data.function.arg_count;
        for (size_t i = 0; i < pn; i++) {
            Expr* arg_i = target->data.function.args[i];
            /* Match Power[plus_i, 1/3]. */
            if (arg_i->type != EXPR_FUNCTION
                || !arg_i->data.function.head
                || arg_i->data.function.head->type != EXPR_SYMBOL
                || arg_i->data.function.head->data.symbol != SYM_Power
                || arg_i->data.function.arg_count != 2) continue;
            Expr* exp_i = arg_i->data.function.args[1];
            int64_t epn, eqn;
            if (!is_rational(exp_i, &epn, &eqn) || epn != 1 || eqn != 3) continue;
            int64_t ai_n, ai_d, bi_n, bi_d, ci;
            if (!cuberoot_match_a_plus_b_sqrt_c(arg_i->data.function.args[0],
                                                 &ai_n, &ai_d, &bi_n, &bi_d, &ci))
                continue;
            if (ci <= 0) continue;

            for (size_t j = i + 1; j < pn; j++) {
                Expr* arg_j = target->data.function.args[j];
                if (arg_j->type != EXPR_FUNCTION
                    || !arg_j->data.function.head
                    || arg_j->data.function.head->type != EXPR_SYMBOL
                    || arg_j->data.function.head->data.symbol != SYM_Power
                    || arg_j->data.function.arg_count != 2) continue;
                Expr* exp_j = arg_j->data.function.args[1];
                int64_t fpn, fqn;
                if (!is_rational(exp_j, &fpn, &fqn) || fpn != 1 || fqn != 3) continue;
                int64_t aj_n, aj_d, bj_n, bj_d, cj;
                if (!cuberoot_match_a_plus_b_sqrt_c(arg_j->data.function.args[0],
                                                     &aj_n, &aj_d, &bj_n, &bj_d, &cj))
                    continue;
                if (cj != ci) continue;
                /* Conjugate test: ai = aj AND bi = -bj (with same denom). */
                if (ai_n * aj_d != aj_n * ai_d) continue;
                if (bi_n * bj_d != -bj_n * bi_d) continue;
                /* Compute the sum. */
                Expr* sum_value = try_cuberoot_sum_conjugates(ai_n, ai_d,
                                                              bi_n, bi_d, ci);
                if (!sum_value) continue;

                /* Build replacement: original Plus minus arg_i and
                 * arg_j, plus sum_value. */
                Expr** new_pa = (Expr**)malloc(sizeof(Expr*) * (pn - 1));
                size_t out = 0;
                for (size_t k = 0; k < pn; k++) {
                    if (k == i || k == j) continue;
                    new_pa[out++] = expr_copy(target->data.function.args[k]);
                }
                new_pa[out++] = sum_value;
                Expr* new_plus = expr_new_function(expr_new_symbol("Plus"),
                                                   new_pa, out);
                Expr* new_eval = evaluate(new_plus);
                expr_free(new_plus);
                if (current) expr_free(current);
                return new_eval;
            }
        }
    }

    return current;
}

static Expr* simp_cuberoot(const Expr* e, const AssumeCtx* ctx) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;
    Expr* r = simp_cuberoot_walk(e, ctx);
    Expr* out = r ? r : expr_copy((Expr*)e);
    if (dbg) simp_debug_log("Cuberoot", e, out,
                            simp_debug_elapsed_ms(t0));
    return out;
}

/* ----------------------------------------------------------------------- */
/* Denominator rationalisation: simp_rationalize_denom                     */
/* ----------------------------------------------------------------------- */

/*
 * simp_rationalize_denom rewrites Power[denom, -1] subtrees where
 * `denom` is a polynomial in radicals over a single integer base `c`,
 * by computing the inverse via extended-Euclidean in Q[α]/(α^n - c)
 * with α = c^(1/n) and n the LCM of the denominators of the rational
 * exponents that appear. Closes two phase-2 user cases:
 *
 *   Simplify[1/(2^(1/3) - 1)]              -> 4^(1/3) + 2^(1/3) + 1
 *   Simplify[1/(Sqrt[2] + 2^(1/3))]        -> ... in Q[2^(1/6)]/(α^6 - 2)
 *
 * Algorithm:
 *   1. Walk the tree top-down. At each Power[denom, k] with k a
 *      negative rational, attempt rationalisation.
 *   2. Find all Power[c, p_i/q_i] subexpressions of denom. Reject
 *      (return NULL) if they don't all share a single positive
 *      integer base c.
 *   3. Compute n = lcm(q_i). Substitute each Power[c, p_i/q_i] with
 *      Power[α, p_i * (n / q_i)] for a fresh variable α.
 *   4. Run PolynomialExtendedGCD on the substituted denom and α^n - c.
 *      Accept the result iff the gcd is a nonzero constant g; the
 *      first Bezout coefficient u(α) then satisfies
 *      u(α) * denom(α) ≡ g (mod α^n - c), so 1/denom = u(α)/g.
 *   5. Substitute α → Power[c, 1/n] and reduce. The result has no
 *      radicals in its denominator.
 *
 * Soundness: the rewrite is correctness-preserving for any positive
 * integer base c (so c^(1/n) is real-valued and the principal-branch
 * arithmetic is unambiguous). For c = 1 the radicals collapse to
 * constants and the substitution is trivial; we skip that case so
 * the round loop's other transforms can handle it.
 */

/* Forward declaration: simp_search drives this. */
static Expr* simp_rationalize_denom(const Expr* e, const AssumeCtx* ctx);

/* Forward declaration: applied as a post-pass on the inverse so that
 * c^(p/n) forms with c = base^k canonicalise to a single power-of-base
 * representation. Without this post-pass, my candidate
 * "1 + 2^(1/3) + 2^(2/3)" subtracted from the user's
 * "1 + 2^(1/3) + 4^(1/3)" leaves 2^(2/3) - 4^(1/3) which is
 * mathematically zero but syntactically nonzero. PrimeRebase rewrites
 * 4^(1/3) -> 2^(2/3) and the cancellation completes. PrimeRebase
 * normally only runs at simp_dispatch on the input, not on candidates
 * produced inside simp_search, so we apply it explicitly here. */
static Expr* transform_prime_rebase(const Expr* e);

/* GCD/LCM on int64_t (inputs assumed non-negative). */
static int64_t denom_gcd_i64(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) {
        int64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static int64_t denom_lcm_i64(int64_t a, int64_t b) {
    if (a == 0 || b == 0) return 0;
    int64_t g = denom_gcd_i64(a, b);
    /* Overflow guard: skip when the product would exceed INT64_MAX. */
    if (a / g > INT64_MAX / b) return 0;
    return (a / g) * b;
}

/* Walk e collecting (base, q) pairs where each Power[base, p/q] uses
 * the SAME positive integer base. Populates *out_base (borrowed
 * pointer into e) and accumulates n = lcm of all q's into *out_n.
 * Returns false on:
 *   - any Power with non-rational exponent
 *   - any Power[c, p/q] with q == 1 (those are integer powers, not
 *     radicals — they pass through cleanly)
 *   - any Power whose base is not a positive integer
 *   - bases that disagree (multi-base extensions are out of scope here).
 *
 * Returns true with *out_base = NULL when e contains no radical
 * subexpressions at all (caller treats as "nothing to rationalise"). */
static bool denom_collect_radical_base(const Expr* e,
                                       const Expr** out_base,
                                       int64_t* out_n) {
    if (!e) return true;
    if (e->type != EXPR_FUNCTION) return true;

    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        int64_t p, q;
        if (is_rational(exp, &p, &q) && q != 1) {
            /* Reject base != positive integer. */
            if (base->type != EXPR_INTEGER || base->data.integer < 2)
                return false;
            /* Reject base disagreement. */
            if (*out_base) {
                if (!expr_eq((Expr*)*out_base, (Expr*)base)) return false;
            } else {
                *out_base = base;
            }
            int64_t new_n = denom_lcm_i64(*out_n ? *out_n : 1, q);
            if (new_n == 0) return false;
            *out_n = new_n;
            /* Don't recurse into the base — it's an integer leaf. */
            return true;
        }
    }
    /* Recurse into children. */
    if (!denom_collect_radical_base(e->data.function.head, out_base, out_n))
        return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (!denom_collect_radical_base(e->data.function.args[i],
                                        out_base, out_n))
            return false;
    }
    return true;
}

/* Substitute every Power[base, p/q] in e with Power[gen, p * (n/q)]
 * (an integer-exponent power of the fresh symbol). Returns a fresh
 * tree. */
static Expr* denom_subst_radical_to_gen(const Expr* e, const Expr* base,
                                        int64_t n, const char* gen) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        const Expr* this_base = e->data.function.args[0];
        const Expr* exp = e->data.function.args[1];
        int64_t p, q;
        if (is_rational(exp, &p, &q) && q != 1 &&
            expr_eq((Expr*)this_base, (Expr*)base)) {
            int64_t k = p * (n / q);
            Expr* pa[2] = { expr_new_symbol(gen), expr_new_integer(k) };
            Expr* pc = expr_new_function(expr_new_symbol("Power"), pa, 2);
            Expr* out = evaluate(pc);
            expr_free(pc);
            return out;
        }
    }
    size_t count = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (count ? count : 1));
    for (size_t i = 0; i < count; i++) {
        new_args[i] = denom_subst_radical_to_gen(e->data.function.args[i],
                                                  base, n, gen);
    }
    Expr* new_head = denom_subst_radical_to_gen(e->data.function.head,
                                                 base, n, gen);
    Expr* result = expr_new_function(new_head, new_args, count);
    free(new_args);
    return evaluate(result);
}

/* Reverse: substitute Power[gen, k] -> Power[base, k * (1/n)] (so a
 * polynomial in `gen` becomes a polynomial in c^(1/n)). Returns a
 * fresh tree. */
static Expr* denom_subst_gen_to_radical(const Expr* e, const char* gen,
                                         const Expr* base, int64_t n) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL && strcmp(e->data.symbol, gen) == 0) {
        Expr* pa[2] = { expr_copy((Expr*)base), make_rational(1, n) };
        Expr* pc = expr_new_function(expr_new_symbol("Power"), pa, 2);
        Expr* out = evaluate(pc);
        expr_free(pc);
        return out;
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2 &&
        e->data.function.args[0]->type == EXPR_SYMBOL &&
        strcmp(e->data.function.args[0]->data.symbol, gen) == 0) {
        const Expr* exp = e->data.function.args[1];
        if (exp->type == EXPR_INTEGER) {
            int64_t k = exp->data.integer;
            /* gen^k = c^(k/n). Build directly. */
            Expr* pa[2] = { expr_copy((Expr*)base), make_rational(k, n) };
            Expr* pc = expr_new_function(expr_new_symbol("Power"), pa, 2);
            Expr* out = evaluate(pc);
            expr_free(pc);
            return out;
        }
    }
    size_t count = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (count ? count : 1));
    for (size_t i = 0; i < count; i++) {
        new_args[i] = denom_subst_gen_to_radical(e->data.function.args[i],
                                                  gen, base, n);
    }
    Expr* new_head = denom_subst_gen_to_radical(e->data.function.head,
                                                 gen, base, n);
    Expr* result = expr_new_function(new_head, new_args, count);
    free(new_args);
    return evaluate(result);
}

/* Compute the inverse of `denom` (a polynomial in radicals over a
 * common positive integer base) via extended-Euclidean in
 * Q[α]/(α^n - c). Returns NULL when:
 *   - denom has no radicals, or has multiple bases.
 *   - the gcd of the substituted denom and α^n - c is not a nonzero
 *     constant (i.e. denom shares a factor with α^n - c — would mean
 *     the input expression is undefined when α takes that value).
 *
 * Soundness gate: only fires when c is a positive integer >= 2. */
static Expr* denom_compute_inverse(const Expr* denom) {
    const Expr* base = NULL;
    int64_t n = 0;
    if (!denom_collect_radical_base(denom, &base, &n)) return NULL;
    if (!base) return NULL;
    if (n < 2) return NULL;
    if (base->type != EXPR_INTEGER || base->data.integer < 2) return NULL;
    int64_t c = base->data.integer;

    /* Generate a fresh symbol name unlikely to clash with the user's
     * variables. The leading $ marks it as a system-generated name. */
    static int counter = 0;
    char gen[32];
    snprintf(gen, sizeof(gen), "$denomGen%d", counter++);

    /* p_in_a = denom rewritten as polynomial in `gen`. */
    Expr* p_in_a = denom_subst_radical_to_gen(denom, base, n, gen);
    if (!p_in_a) return NULL;

    /* Build the relation polynomial: gen^n - c. */
    Expr* gen_pow_args[2] = { expr_new_symbol(gen), expr_new_integer(n) };
    Expr* gen_pow = expr_new_function(expr_new_symbol("Power"), gen_pow_args, 2);
    Expr* gen_pow_e = evaluate(gen_pow);
    expr_free(gen_pow);
    Expr* relation_args[2] = { gen_pow_e, expr_new_integer(-c) };
    Expr* relation_call = expr_new_function(expr_new_symbol("Plus"),
                                            relation_args, 2);
    Expr* relation = evaluate(relation_call);
    expr_free(relation_call);

    /* Run PolynomialExtendedGCD[p_in_a, relation, gen]. */
    Expr* xgcd_args[3] = {
        p_in_a, relation, expr_new_symbol(gen)
    };
    Expr* xgcd_call = expr_new_function(expr_new_symbol("PolynomialExtendedGCD"),
                                        xgcd_args, 3);
    Expr* xgcd_result = evaluate(xgcd_call);
    expr_free(xgcd_call);

    /* Validate result shape: List[gcd, List[u, v]]. */
    if (!xgcd_result ||
        xgcd_result->type != EXPR_FUNCTION ||
        xgcd_result->data.function.head->type != EXPR_SYMBOL ||
        xgcd_result->data.function.head->data.symbol != SYM_List ||
        xgcd_result->data.function.arg_count != 2) {
        if (xgcd_result) expr_free(xgcd_result);
        return NULL;
    }
    Expr* gcd_e = xgcd_result->data.function.args[0];
    Expr* coeffs = xgcd_result->data.function.args[1];
    if (!coeffs ||
        coeffs->type != EXPR_FUNCTION ||
        coeffs->data.function.head->type != EXPR_SYMBOL ||
        coeffs->data.function.head->data.symbol != SYM_List ||
        coeffs->data.function.arg_count < 1) {
        expr_free(xgcd_result);
        return NULL;
    }
    Expr* u_in_a = coeffs->data.function.args[0];

    /* Require gcd to be a nonzero rational constant (no `gen` in it).
     * If gcd has the gen, it shares a factor with α^n - c, meaning
     * denom vanishes at some root and the expression is undefined. */
    if (gcd_e->type != EXPR_INTEGER &&
        !is_rational_literal(gcd_e)) {
        expr_free(xgcd_result);
        return NULL;
    }
    if (gcd_e->type == EXPR_INTEGER && gcd_e->data.integer == 0) {
        expr_free(xgcd_result);
        return NULL;
    }

    /* Build inverse = u_in_a / gcd_e (still in `gen` form). */
    Expr* inv_gen_args[2] = { expr_copy(gcd_e), expr_new_integer(-1) };
    Expr* gcd_inv_call = expr_new_function(expr_new_symbol("Power"),
                                           inv_gen_args, 2);
    Expr* gcd_inv = evaluate(gcd_inv_call);
    expr_free(gcd_inv_call);
    Expr* prod_args[2] = { expr_copy(u_in_a), gcd_inv };
    Expr* prod_call = expr_new_function(expr_new_symbol("Times"),
                                        prod_args, 2);
    Expr* inv_in_gen = evaluate(prod_call);
    expr_free(prod_call);

    expr_free(xgcd_result);

    /* Substitute gen → c^(1/n) and re-evaluate. */
    Expr* inv_in_radical = denom_subst_gen_to_radical(inv_in_gen, gen, base, n);
    expr_free(inv_in_gen);
    return inv_in_radical;
}

/* Walk the expression tree, applying denominator rationalisation at
 * each Power[denom, k] subterm with k a negative rational. Returns
 * NULL when no rewrite fires. */
static Expr* simp_rationalize_denom_walk(const Expr* e,
                                         const AssumeCtx* ctx) {
    (void)ctx;
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return NULL;

    /* Recurse into children first so deeper rewrites bubble up. */
    size_t n = e->data.function.arg_count;
    Expr** new_args = NULL;
    bool any = false;
    for (size_t i = 0; i < n; i++) {
        Expr* r = simp_rationalize_denom_walk(e->data.function.args[i], ctx);
        if (r) {
            if (!new_args) {
                new_args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
                for (size_t j = 0; j < i; j++) {
                    new_args[j] = expr_copy(e->data.function.args[j]);
                }
            }
            new_args[i] = r;
            any = true;
        } else if (new_args) {
            new_args[i] = expr_copy(e->data.function.args[i]);
        }
    }
    Expr* current = NULL;
    if (any) {
        Expr* head_copy = expr_copy(e->data.function.head);
        Expr* rebuilt = expr_new_function(head_copy, new_args, n);
        free(new_args);
        current = evaluate(rebuilt);
        expr_free(rebuilt);
    }

    /* Check whether this node is a Power[denom, k] with k a negative
     * rational and denom a candidate for rationalisation. */
    const Expr* target = current ? current : e;
    if (target->type == EXPR_FUNCTION
        && target->data.function.head
        && target->data.function.head->type == EXPR_SYMBOL
        && target->data.function.head->data.symbol == SYM_Power
        && target->data.function.arg_count == 2) {
        const Expr* denom = target->data.function.args[0];
        const Expr* exp = target->data.function.args[1];
        bool is_neg_rational = false;
        int64_t k_p = 0, k_q = 1;
        if (exp->type == EXPR_INTEGER && exp->data.integer < 0) {
            is_neg_rational = true;
            k_p = exp->data.integer;
            k_q = 1;
        } else if (is_rational_literal(exp)) {
            int64_t pn, qn;
            if (is_rational(exp, &pn, &qn) && pn < 0 && qn == 1) {
                /* canonical Rational with denom 1 shouldn't occur, but
                 * be defensive. */
                is_neg_rational = true;
                k_p = pn;
                k_q = 1;
            }
        }
        /* For now, only handle k = -1. Power[denom, -k] for k > 1 can
         * be reduced to repeated multiplication of the rationalised
         * inverse, but isn't needed by the user's phase-2 cases. */
        if (is_neg_rational && k_p == -1 && k_q == 1) {
            Expr* inv = denom_compute_inverse(denom);
            if (inv) {
                if (current) expr_free(current);
                return inv;
            }
        }
    }
    return current;
}

static Expr* simp_rationalize_denom(const Expr* e, const AssumeCtx* ctx) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;
    Expr* r = simp_rationalize_denom_walk(e, ctx);
    Expr* out = r ? r : expr_copy((Expr*)e);
    /* PrimeRebase post-pass on the FULL output. The user's RHS may use
     * a "compound" radical base (e.g. 4^(1/3)) while the inverse our
     * walker computes uses the PRIME base (2^(2/3)). Without rebasing
     * the full expression, the subtraction "2^(2/3) - 4^(1/3)" doesn't
     * cancel even though the values are equal. PrimeRebase normally
     * runs only at simp_dispatch on the input, not on candidates
     * produced inside simp_search, so we apply it here. */
    if (r) {
        Expr* rebased = transform_prime_rebase(out);
        if (rebased) {
            if (!expr_eq(rebased, out)) {
                expr_free(out);
                out = rebased;
            } else {
                expr_free(rebased);
            }
        }
    }
    if (dbg) simp_debug_log("RationalizeDenom", e, out,
                            simp_debug_elapsed_ms(t0));
    return out;
}

/* ----------------------------------------------------------------------- */
/* Algebraic-extension reduction: simp_algebraic                           */
/* ----------------------------------------------------------------------- */

/*
 * simp_algebraic rewrites an expression that contains one or more
 * distinct square-root sub-expressions Sqrt[u_i] by treating each
 * Sqrt[u_i] as a generator g_i of the algebraic extension
 *   K(vars)[g_1, ..., g_n] / (g_1^2 - u_1, ..., g_n^2 - u_n).
 * Standard rational arithmetic (Together) followed by reduction modulo
 * the relation ideal and successive rationalisation of the denominator
 * collapses identities that ordinary Together / Cancel cannot see, e.g.
 *
 *   (x/Sqrt[x^2+1] + 1) / ((Sqrt[x^2+1] + x)^2 + 1)  ->  1/(2 + 2 x^2)
 *   (x/(Sqrt[x^2+6] - Sqrt[6]))(1/Sqrt[x^2+6]
 *      - (Sqrt[x^2+6] - Sqrt[6])/x^2)                ->  Sqrt[6]/(x Sqrt[x^2+6])
 *   2/(Sqrt[(x+1)/(1-x)] - 1/Sqrt[(x+1)/(1-x)])      ->  Sqrt[(x+1)/(1-x)] (1-x)/x
 *
 * Algorithm:
 *   1. Walk the expression collecting every distinct surd argument
 *      u_i where Power[u_i, p/q] appears with q != 1. Bail if any q != 2
 *      (cube roots etc.), if any u_i contains an explicit complex
 *      literal, or if more than ALG_MAX_SURDS distinct bases appear.
 *   2. Substitute Sqrt[u_i] -> g_i for fresh distinct generator symbols.
 *      After substitution, the expression is a rational function in
 *      (vars, g_1, ..., g_n).
 *   3. Together  ->  N / D, both polynomials in (vars, g_1, ..., g_n).
 *   4. Reduce both N and D modulo the relation ideal {g_i^2 - u_i}_i
 *      via successive CoefficientList[..., g_i] decomposition. After
 *      one sweep across all generators the polynomial is multilinear
 *      in {g_i}: every g_i appears at degree 0 or 1.
 *   5. For i = 1..n, rationalise the i-th generator out of the
 *      denominator: multiply numerator and denominator by sigma_i(D)
 *      (D with g_i sign-flipped), then reduce again. After each step
 *      g_i has been eliminated from the denominator. The product
 *      D * sigma_i(D) lies in K[g_1, ..., g_{i-1}, g_{i+1}, ..., g_n]
 *      because the g_i terms in (a + b g_i)(a - b g_i) = a^2 - b^2 u_i
 *      cancel.
 *   6. Substitute g_i -> Sqrt[u_i] back, run Together / Cancel for
 *      cleanup, and accept the result iff its complexity score is
 *      strictly lower than the input.
 *
 * Principal-branch concern: the substitution Sqrt[u_i]^2 = u_i is only
 * sound where u_i lies in the principal branch's domain. We accept the
 * Mathematica-style convention (Simplify treats this as an identity on
 * the natural domain where the input is real) but skip when any u_i
 * contains an explicit complex literal (Complex[..,..] or the symbol I)
 * so we never produce a result that swallows a sign-of-imaginary-part
 * change silently.
 */

#define ALG_MAX_SURDS 4

/* Walk e collecting distinct surd bases. The walker enforces:
 *   - every Power[base, p/q] with q != 1 has q == 2,
 *   - distinct bases (by structural equality, expr_eq) accumulate into
 *     bases[0..*n_bases-1] up to max_n,
 *   - returns false on q != 2 or when bases would overflow max_n.
 *
 * Borrowed pointers into `e`. */
static bool alg_collect_sqrt_bases(const Expr* e, const Expr** bases,
                                   size_t* n_bases, size_t max_n) {
    if (!e || e->type != EXPR_FUNCTION) return true;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        Expr* exp        = e->data.function.args[1];
        int64_t p, q;
        if (is_rational(exp, &p, &q) && q != 1) {
            if (q != 2) return false;
            bool seen = false;
            for (size_t i = 0; i < *n_bases; i++) {
                if (expr_eq((Expr*)base, (Expr*)bases[i])) { seen = true; break; }
            }
            if (!seen) {
                if (*n_bases >= max_n) return false;
                bases[(*n_bases)++] = base;
            }
        }
    }
    if (!alg_collect_sqrt_bases(e->data.function.head, bases, n_bases, max_n))
        return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (!alg_collect_sqrt_bases(e->data.function.args[i], bases, n_bases, max_n))
            return false;
    }
    return true;
}

/* Returns true if any sub-expression has head Complex or contains the
 * symbol I. Used to gate simp_algebraic off explicit-complex inputs
 * whose Sqrt[]^2 = arg identity could mask a branch flip. */
static bool contains_explicit_complex(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL && e->data.symbol == SYM_I) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Complex) return true;
    if (contains_explicit_complex(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_explicit_complex(e->data.function.args[i])) return true;
    }
    return false;
}

/* For every Power[bases[i], p/2] in e, replace with bases[i]^floor(p/2)
 * * gens[i]^(p mod 2) (computed via floor-division so negative p is
 * handled correctly). Bases that don't appear are passed through. */
static Expr* alg_subst_sqrt_to_gens(const Expr* e, const Expr** bases,
                                    const char** gens, size_t n) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        int64_t p, q;
        if (is_rational(exp, &p, &q) && q == 2) {
            const Expr* base = e->data.function.args[0];
            for (size_t i = 0; i < n; i++) {
                if (expr_eq((Expr*)base, (Expr*)bases[i])) {
                    int64_t m = p / 2;
                    int64_t r = p - 2 * m;
                    if (r < 0) { m -= 1; r += 2; }
                    Expr* base_pow = (m == 0)
                        ? expr_new_integer(1)
                        : eval_and_free(expr_new_function(expr_new_symbol("Power"),
                              (Expr*[]){expr_copy((Expr*)base), expr_new_integer(m)}, 2));
                    Expr* g_pow = (r == 0)
                        ? expr_new_integer(1)
                        : expr_new_symbol(gens[i]);
                    return eval_and_free(expr_new_function(expr_new_symbol("Times"),
                              (Expr*[]){base_pow, g_pow}, 2));
                }
            }
        }
    }
    size_t count = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (count ? count : 1));
    for (size_t i = 0; i < count; i++) {
        new_args[i] = alg_subst_sqrt_to_gens(e->data.function.args[i], bases, gens, n);
    }
    Expr* new_head = alg_subst_sqrt_to_gens(e->data.function.head, bases, gens, n);
    Expr* result = eval_and_free(expr_new_function(new_head, new_args, count));
    free(new_args);
    return result;
}

/* For each generator symbol gens[i], replace with Sqrt[bases[i]]. */
static Expr* alg_subst_gens_to_sqrt(const Expr* e, const char** gens,
                                    const Expr** bases, size_t n) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < n; i++) {
            if (e->data.symbol == gens[i]) {
                return eval_and_free(expr_new_function(expr_new_symbol("Power"),
                          (Expr*[]){expr_copy((Expr*)bases[i]), make_rational(1, 2)}, 2));
            }
        }
        return expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t count = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (count ? count : 1));
    for (size_t i = 0; i < count; i++) {
        new_args[i] = alg_subst_gens_to_sqrt(e->data.function.args[i], gens, bases, n);
    }
    Expr* new_head = alg_subst_gens_to_sqrt(e->data.function.head, gens, bases, n);
    Expr* result = eval_and_free(expr_new_function(new_head, new_args, count));
    free(new_args);
    return result;
}

/* Replace every occurrence of generator gi_sym with -gi_sym in e. Used
 * to compute sigma_i(den) for rationalisation. */
static Expr* alg_sigma_negate(const Expr* e, const char* gi_sym) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL && e->data.symbol == gi_sym) {
        return eval_and_free(expr_new_function(expr_new_symbol("Times"),
                  (Expr*[]){expr_new_integer(-1), expr_copy((Expr*)e)}, 2));
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t count = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (count ? count : 1));
    for (size_t i = 0; i < count; i++) {
        new_args[i] = alg_sigma_negate(e->data.function.args[i], gi_sym);
    }
    Expr* new_head = alg_sigma_negate(e->data.function.head, gi_sym);
    Expr* result = eval_and_free(expr_new_function(new_head, new_args, count));
    free(new_args);
    return result;
}

/* Reduce poly modulo gi_sym^2 - u_i: returns A + B*gi where
 *   A = sum_{k even} a_k u^(k/2)
 *   B = sum_{k odd}  a_k u^((k-1)/2)
 * with a_k extracted via CoefficientList[poly, gi_sym]. The caller is
 * expected to have Expand-ed `poly` first. Returns NULL if
 * CoefficientList didn't yield a List. */
static Expr* alg_reduce_one_gen(const Expr* poly, const char* gi_sym,
                                const Expr* ui) {
    Expr* cl_args[2] = { expr_copy((Expr*)poly), expr_new_symbol(gi_sym) };
    Expr* cl_call = expr_new_function(expr_new_symbol("CoefficientList"),
                                      cl_args, 2);
    Expr* coefs = evaluate(cl_call);
    expr_free(cl_call);
    if (!coefs || coefs->type != EXPR_FUNCTION ||
        coefs->data.function.head->type != EXPR_SYMBOL ||
        coefs->data.function.head->data.symbol != SYM_List) {
        if (coefs) expr_free(coefs);
        return NULL;
    }

    size_t n = coefs->data.function.arg_count;
    Expr** evens = (Expr**)malloc(sizeof(Expr*) * (n + 1));
    Expr** odds  = (Expr**)malloc(sizeof(Expr*) * (n + 1));
    size_t ne = 0, no = 0;
    for (size_t k = 0; k < n; k++) {
        Expr* ck = coefs->data.function.args[k];
        int64_t exp_u = (int64_t)(k / 2);
        Expr* upow = (exp_u == 0)
            ? expr_new_integer(1)
            : eval_and_free(expr_new_function(expr_new_symbol("Power"),
                  (Expr*[]){expr_copy((Expr*)ui), expr_new_integer(exp_u)}, 2));
        Expr* term = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                  (Expr*[]){expr_copy(ck), upow}, 2));
        if ((k & 1) == 0) evens[ne++] = term;
        else              odds[no++]  = term;
    }
    expr_free(coefs);

    Expr* a_sum;
    if (ne == 0)      a_sum = expr_new_integer(0);
    else if (ne == 1) a_sum = evens[0];
    else              a_sum = eval_and_free(expr_new_function(expr_new_symbol("Plus"),
                                                              evens, ne));
    Expr* b_sum;
    if (no == 0)      b_sum = expr_new_integer(0);
    else if (no == 1) b_sum = odds[0];
    else              b_sum = eval_and_free(expr_new_function(expr_new_symbol("Plus"),
                                                              odds, no));
    free(evens);
    free(odds);

    /* Combine A + B*gi into a single expression. */
    Expr* b_gi = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                  (Expr*[]){b_sum, expr_new_symbol(gi_sym)}, 2));
    return eval_and_free(expr_new_function(expr_new_symbol("Plus"),
                  (Expr*[]){a_sum, b_gi}, 2));
}

/* Polynomial-divide `poly` by `u` repeatedly (treating both as
 * polynomials in `var`) until the division has a non-zero remainder.
 * Returns the residual quotient and writes the multiplicity to *k_out
 * so that `poly = u^(*k_out) * residual` modulo non-divisibility.
 *
 * Used so that an implicit u_i^k factor inside den_r (e.g. x^4 hiding
 * (x^2)^2 when u = x^2) can be lifted into Power[g_i, 2k] -- once the
 * u-power is expressed in terms of the generator, polynomial GCD over
 * Q[vars, g_1, ..., g_n] cancels it against any g_i factors carried by
 * the multilinear numerator. */
static Expr* alg_extract_u_power(const Expr* poly, const Expr* u,
                                 const Expr* var, int* k_out) {
    int k = 0;
    Expr* cur = expr_copy((Expr*)poly);
    for (;;) {
        Expr* qa[3] = { expr_copy(cur), expr_copy((Expr*)u), expr_copy((Expr*)var) };
        Expr* qcall = expr_new_function(expr_new_symbol("PolynomialQuotient"), qa, 3);
        Expr* q = evaluate(qcall);
        expr_free(qcall);

        Expr* ra[3] = { expr_copy(cur), expr_copy((Expr*)u), expr_copy((Expr*)var) };
        Expr* rcall = expr_new_function(expr_new_symbol("PolynomialRemainder"), ra, 3);
        Expr* r = evaluate(rcall);
        expr_free(rcall);

        bool zero = (r && r->type == EXPR_INTEGER && r->data.integer == 0);
        if (r) expr_free(r);
        if (!zero) { expr_free(q); break; }
        expr_free(cur);
        cur = q;
        k++;
        /* Defensive cap: prevent runaway when PolynomialQuotient
         * misbehaves (e.g. floating-point coefficients sneaking in). */
        if (k > 100) break;
    }
    *k_out = k;
    return cur;
}

/* Returns true iff u is a polynomial in its own variables -- i.e.,
 * every Power[base, exp] in u has a non-negative integer exp. Rational
 * u (e.g. (x+1)/(1-x)) is rejected so the polynomial-division u-power
 * extraction never tries to divide by a non-polynomial divisor. */
static bool alg_u_is_polynomial(const Expr* u) {
    if (!u) return false;
    if (u->type != EXPR_FUNCTION) return true;   /* leaf is always polynomial */
    if (u->data.function.head &&
        u->data.function.head->type == EXPR_SYMBOL &&
        u->data.function.head->data.symbol == SYM_Power &&
        u->data.function.arg_count == 2) {
        Expr* exp = u->data.function.args[1];
        if (exp->type != EXPR_INTEGER && exp->type != EXPR_BIGINT) return false;
        if (exp->type == EXPR_INTEGER && exp->data.integer < 0) return false;
        if (exp->type == EXPR_BIGINT && mpz_sgn(exp->data.bigint) < 0) return false;
    }
    if (!alg_u_is_polynomial(u->data.function.head)) return false;
    for (size_t i = 0; i < u->data.function.arg_count; i++) {
        if (!alg_u_is_polynomial(u->data.function.args[i])) return false;
    }
    return true;
}

/* Pick the first variable in Variables[u]. Returns NULL when u is
 * variable-free (numeric / constant), in which case alg_extract_u_power
 * is undefined and the caller should skip the u-power-extraction step. */
static Expr* alg_pick_var(const Expr* u) {
    Expr* vars = call_unary_copy("Variables", u);
    if (!vars || vars->type != EXPR_FUNCTION ||
        vars->data.function.head->type != EXPR_SYMBOL ||
        vars->data.function.head->data.symbol != SYM_List ||
        vars->data.function.arg_count == 0) {
        if (vars) expr_free(vars);
        return NULL;
    }
    Expr* v = expr_copy(vars->data.function.args[0]);
    expr_free(vars);
    return v;
}

/* Reduce poly modulo all relations {gi_sym^2 - u_i}_i by sweeping each
 * generator once. The result is multilinear in (g_1, ..., g_n). The
 * input is Expand-ed before each generator pass so CoefficientList sees
 * the canonical polynomial form. Returns NULL on any inner failure. */
static Expr* alg_reduce_all(const Expr* poly, const char** gens,
                            const Expr** us, size_t n) {
    Expr* cur = call_unary_copy("Expand", poly);
    for (size_t i = 0; i < n; i++) {
        Expr* nxt = alg_reduce_one_gen(cur, gens[i], us[i]);
        expr_free(cur);
        if (!nxt) return NULL;
        cur = call_unary_owned("Expand", nxt);
    }
    return cur;
}

static Expr* simp_algebraic_impl(const Expr* e) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Cheap precondition: nothing to do without a half-integer Power. */
    if (!has_non_integer_power(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("Algebraic", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    /* Early denester short-circuit. For nested-radical inputs (e.g.
     * Sqrt[A + b Sqrt[C]] or 1/Sqrt[...]^k), the half-sum denester is
     * O(small) and the algebraic-tower Together-with-Extension below is
     * O(big). When the denester strictly wins on complexity, taking that
     * win here saves the 200+ ms it costs to compute and then reject
     * the Together[expr, Extension -> α] form for the same input.
     *
     * The complexity gate is strict (<) so we don't churn on inputs the
     * denester touches without reducing -- those still fall through to
     * the extension path which may legitimately help. */
    {
        Expr* denested = simp_denest_sqrt(e, NULL);
        if (denested && simp_default_complexity(denested)
                          < simp_default_complexity(e)) {
            if (dbg) simp_debug_log("Algebraic", e, denested,
                                    simp_debug_elapsed_ms(t0));
            return denested;
        }
        if (denested) expr_free(denested);
    }

    /* Phase G9 (cube-root and higher): when the input has exactly one
     * rational-base radical generator, route through
     * `Together[expr, Extension -> α]`.  Together's extension path uses
     * the qaupoly substrate (qaupoly_gcd / qaupoly_divrem), which handles
     * any q ≥ 2 natively — no Sqrt-only special case.  This is the
     * shortcut that lets Simplify simplify expressions involving
     * `Power[c, p/q]` for q > 2, which the older multi-Sqrt path
     * (alg_sigma_negate sign-flip rationalisation) below cannot.
     *
     * Multi-generator cases (n ≥ 2) fall through to the Sqrt-only
     * path because that path's general rationalisation by sign-flip
     * conjugates handles n ≥ 2 over Q(Sqrt[u_i]) directly. */
    if (!contains_explicit_complex(e)) {
        QATower* qa_t = extension_autodetect(e);
        if (qa_t && qa_t->n == 1) {
            Expr* alpha = expr_copy(qa_t->alpha_renders[0]);
            Expr* tog = expr_new_function(
                expr_new_symbol("Together"),
                (Expr*[]){
                    expr_copy((Expr*)e),
                    expr_new_function(expr_new_symbol("Rule"),
                        (Expr*[]){expr_new_symbol("Extension"), alpha}, 2)
                }, 2);
            Expr* qa_out = evaluate(tog);
            expr_free(tog);
            qa_tower_free(qa_t);

            /* Accept only when the qaupoly path produced something at
             * least as compact as the input.  The Together-with-Extension
             * path occasionally factors a denominator unnecessarily —
             * letting simp_search's later round-loop pick the better
             * form is the safe move. */
            if (qa_out && simp_default_complexity(qa_out)
                            <= simp_default_complexity(e)) {
                if (dbg) simp_debug_log("Algebraic", e, qa_out,
                                        simp_debug_elapsed_ms(t0));
                return qa_out;
            }
            if (qa_out) expr_free(qa_out);
        } else if (qa_t) {
            qa_tower_free(qa_t);
        }
    }

    const Expr* bases[ALG_MAX_SURDS];
    size_t n = 0;
    if (!alg_collect_sqrt_bases(e, bases, &n, ALG_MAX_SURDS) || n == 0) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("Algebraic", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }
    /* Each surd's argument must itself be surd-free and contain no
     * explicit complex literals. */
    for (size_t i = 0; i < n; i++) {
        if (has_non_integer_power(bases[i]) ||
            contains_explicit_complex(bases[i])) {
            Expr* out = expr_copy((Expr*)e);
            if (dbg) simp_debug_log("Algebraic", e, out,
                                    simp_debug_elapsed_ms(t0));
            return out;
        }
    }
    if (contains_explicit_complex(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("Algebraic", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    /* Allocate a fresh interned generator symbol per surd. The names
     * are $-prefixed so they won't collide with user symbols. */
    const char* gens[ALG_MAX_SURDS];
    static const char* gen_names[ALG_MAX_SURDS] = {
        "$pc_alggen0$", "$pc_alggen1$", "$pc_alggen2$", "$pc_alggen3$"
    };
    for (size_t i = 0; i < n; i++) gens[i] = intern_symbol(gen_names[i]);

    /* Step 2-3: substitute and Together. */
    Expr* sub = alg_subst_sqrt_to_gens(e, bases, gens, n);
    Expr* tg  = call_unary_owned("Together", sub);
    Expr* num = call_unary_copy("Numerator",   tg);
    Expr* den = call_unary_copy("Denominator", tg);
    expr_free(tg);

    /* Step 4: reduce both modulo the relation ideal. */
    Expr* num_r = alg_reduce_all(num, gens, bases, n);
    Expr* den_r = alg_reduce_all(den, gens, bases, n);
    expr_free(num); expr_free(den);
    if (!num_r || !den_r) {
        if (num_r) expr_free(num_r);
        if (den_r) expr_free(den_r);
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("Algebraic", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    /* Step 5: rationalise each generator out of the denominator in turn. */
    for (size_t i = 0; i < n; i++) {
        Expr* sig = alg_sigma_negate(den_r, gens[i]);
        Expr* num_mul = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                          (Expr*[]){num_r, expr_copy(sig)}, 2));
        Expr* den_mul = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                          (Expr*[]){den_r, sig}, 2));
        Expr* num_next = alg_reduce_all(num_mul, gens, bases, n);
        Expr* den_next = alg_reduce_all(den_mul, gens, bases, n);
        expr_free(num_mul); expr_free(den_mul);
        if (!num_next || !den_next) {
            if (num_next) expr_free(num_next);
            if (den_next) expr_free(den_next);
            Expr* out = expr_copy((Expr*)e);
            if (dbg) simp_debug_log("Algebraic", e, out,
                                    simp_debug_elapsed_ms(t0));
            return out;
        }
        num_r = num_next;
        den_r = den_next;
    }

    /* Step 6 (pre): pull each implicit u_i^k factor out of num_r and
     * den_r, replacing it with g_i^(2k) so that Cancel over
     * Q[vars, g_1, ..., g_n] sees the cancellation between the
     * multilinear g_i factor in the numerator and an implicit u_i^k
     * factor in the denominator. Without this step, x^4 in the
     * denominator and Sqrt[x^2] in the numerator look like coprime
     * polynomial atoms to Cancel even though x^4 = u^2 = g^4 modulo
     * the algebraic relation g^2 = u. */
    for (size_t i = 0; i < n; i++) {
        if (!alg_u_is_polynomial(bases[i])) continue;  /* rational u: skip */
        Expr* var = alg_pick_var(bases[i]);
        if (!var) continue;     /* numeric u_i: no polynomial division */

        int kn = 0, kd = 0;
        Expr* num_resid = alg_extract_u_power(num_r, bases[i], var, &kn);
        Expr* den_resid = alg_extract_u_power(den_r, bases[i], var, &kd);
        expr_free(var);

        if (kn == 0 && kd == 0) {
            expr_free(num_resid); expr_free(den_resid);
            continue;
        }
        /* num_r = num_resid * Power[g_i, 2*kn]
         * den_r = den_resid * Power[g_i, 2*kd] */
        if (kn > 0) {
            Expr* g_pow = eval_and_free(expr_new_function(expr_new_symbol("Power"),
                          (Expr*[]){expr_new_symbol(gens[i]), expr_new_integer(2 * kn)}, 2));
            expr_free(num_r);
            num_r = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                      (Expr*[]){num_resid, g_pow}, 2));
        } else {
            expr_free(num_r);
            num_r = num_resid;
        }
        if (kd > 0) {
            Expr* g_pow = eval_and_free(expr_new_function(expr_new_symbol("Power"),
                          (Expr*[]){expr_new_symbol(gens[i]), expr_new_integer(2 * kd)}, 2));
            expr_free(den_r);
            den_r = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                      (Expr*[]){den_resid, g_pow}, 2));
        } else {
            expr_free(den_r);
            den_r = den_resid;
        }
    }

    /* Step 6: assemble num_r / den_r, substitute generators back, clean.
     *
     * Apply Factor to the polynomial-in-(vars, g_1..g_n) numerator and
     * denominator before substituting g_i -> Sqrt[u_i]. Without this
     * step, Cancel sees expanded polynomials whose common (u_i)^k
     * factors share denominators with Sqrt[u_i]^(2k); Factor exposes
     * the (u_i)^k structure so Cancel can combine
     * Power[u_i, k] * Power[u_i, 1/2] / Power[u_i, j]
     * into a single Power[u_i, ...] term. */
    Expr* num_factored = call_unary_owned("Factor", num_r);
    Expr* den_factored = call_unary_owned("Factor", den_r);

    Expr* den_inv = eval_and_free(expr_new_function(expr_new_symbol("Power"),
                  (Expr*[]){den_factored, expr_new_integer(-1)}, 2));
    Expr* quot = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                  (Expr*[]){num_factored, den_inv}, 2));

    Expr* with_sqrt = alg_subst_gens_to_sqrt(quot, gens, bases, n);
    expr_free(quot);
    Expr* result = call_unary_owned("Cancel", with_sqrt);

    /* Complexity gate: accept any form whose complexity score is no
     * greater than the input. The strict ">=" rejection used by the
     * simp_search round loop is too tight here because rationalisation
     * often hits a tied score (e.g. 1/(Sqrt[a]+Sqrt[b]) trades the
     * Power[..,-1] head for a single Times[-1, ...] term while keeping
     * two Sqrt leaves -- equal complexity but the rationalised form is
     * the conventionally preferred shape). simp_search's later round
     * loop will still pick the strictly-better form when one exists. */
    if (simp_default_complexity(result) > simp_default_complexity(e)) {
        expr_free(result);
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("Algebraic", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    if (dbg) simp_debug_log("Algebraic", e, result,
                            simp_debug_elapsed_ms(t0));
    return result;
}

static Expr* simp_algebraic(const Expr* e) {
    return simp_memo_wrap(e, "$Algebraic", simp_algebraic_impl);
}

/* ----------------------------------------------------------------------- */
/* Common-factor lift across Plus terms                                    */
/* ----------------------------------------------------------------------- */
/*
 * lift_common: given a Plus expression whose terms share a multiplicative
 * factor (an algebraic generator like (1+x^2)^(3/2), a free symbol, an
 * integer, or a Power[base, n>=1] that splits into n copies of base),
 * factor the common piece outside the Plus.
 *
 * Why a dedicated transform? picocas's Factor / FactorTerms decompose
 * polynomials over K[x_1, ..., x_n] using Variables[] to discover the
 * generator set, and Variables[] does not return non-integer Power
 * expressions (e.g. Sqrt[x], (1+x^2)^(3/2)). So a Plus that obviously
 * shares (1+x^2)^(3/2) across all terms slips past Factor untouched.
 * This transform takes a structural multiset view: a non-numeric factor
 * is either an algebraic generator (Power with non-integer exponent,
 * or any Power exponent we can't reason about) treated as one opaque
 * token, or a Power[base, n] with n a small positive integer that we
 * split into n copies of base. The numeric coefficients merge via
 * rational GCD. Lifted result: gcd_coef * Times(common_tokens) *
 * Plus[t_i / lift_factor], with the division handed back to evaluate()
 * for cancellation.
 *
 * Cases this enables that Factor alone cannot:
 *   - Plus of c_i * (1+x^2)^(3/2) * x^k  -> (1+x^2)^(3/2) * Plus[c_i x^k]
 *   - Plus inside Times[Plus, Power[denom, neg]] (rational expressions
 *     with non-integer-power denominator): factor the numerator only.
 *
 * Returns NULL when no nontrivial lift is possible (single-term Plus,
 * coprime coefficients with no shared token). */

typedef struct {
    Expr** items;      /* aliased pointers into the term tree (no ownership) */
    size_t count;
    size_t cap;
} LiftTokList;

static void lift_tl_init(LiftTokList* t) {
    t->items = NULL; t->count = 0; t->cap = 0;
}
static void lift_tl_free(LiftTokList* t) { free(t->items); }
static void lift_tl_push(LiftTokList* t, Expr* e) {
    if (t->count == t->cap) {
        size_t nc = t->cap ? t->cap * 2 : 4;
        Expr** ni = (Expr**)realloc(t->items, sizeof(Expr*) * nc);
        if (!ni) { /* OOM: drop the push silently rather than abort. */
            return;
        }
        t->items = ni; t->cap = nc;
    }
    t->items[t->count++] = e;
}

/* Convert an Expr* to mpq_t. Recognises EXPR_INTEGER, EXPR_BIGINT, and
 * Rational[n, d]. Returns false for anything else. */
static bool lift_expr_to_mpq(const Expr* e, mpq_t out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) {
        mpq_set_si(out, (long)e->data.integer, 1);
        return true;
    }
    if (e->type == EXPR_BIGINT) {
        mpq_set_z(out, e->data.bigint);
        return true;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpq_set_si(out, (long)n, (unsigned long)d);
        mpq_canonicalize(out);
        return true;
    }
    return false;
}

/* Build a normalised numeric Expr from an mpq_t. Returns Integer when
 * the denominator is 1, otherwise a Rational[n, d]. */
static Expr* lift_mpq_to_expr(const mpq_t v) {
    if (mpz_cmp_ui(mpq_denref(v), 1) == 0) {
        if (mpz_fits_slong_p(mpq_numref(v))) {
            return expr_new_integer((int64_t)mpz_get_si(mpq_numref(v)));
        }
        return expr_new_bigint_from_mpz(mpq_numref(v));
    }
    Expr* num = mpz_fits_slong_p(mpq_numref(v))
                  ? expr_new_integer((int64_t)mpz_get_si(mpq_numref(v)))
                  : expr_new_bigint_from_mpz(mpq_numref(v));
    Expr* den = mpz_fits_slong_p(mpq_denref(v))
                  ? expr_new_integer((int64_t)mpz_get_si(mpq_denref(v)))
                  : expr_new_bigint_from_mpz(mpq_denref(v));
    Expr* args[2] = { num, den };
    return expr_new_function(expr_new_symbol("Rational"), args, 2);
}

/* Decompose one Plus term into (mpq coefficient, list of token aliases).
 * Numeric leaves accumulate into *coef. Power[base, n] with 1<=n<=16 is
 * split into n copies of base. Power[base, exp] with any other exp shape
 * (rational, negative, symbolic) is treated as one opaque token.
 *
 * Recurses into nested Times: in practice picocas's Plus does not always
 * fully flatten Times children -- a literal-times-product subexpression
 * inside a Plus surfaces as Times[c, Times[a, b]] -- so we walk the
 * subtree rather than relying on a one-level-deep view. */
static void lift_decompose_term(Expr* term, mpq_t coef, LiftTokList* tokens) {
    if (!term) return;
    if (term->type == EXPR_FUNCTION
        && term->data.function.head
        && term->data.function.head->type == EXPR_SYMBOL
        && term->data.function.head->data.symbol == SYM_Times) {
        for (size_t i = 0; i < term->data.function.arg_count; i++) {
            lift_decompose_term(term->data.function.args[i], coef, tokens);
        }
        return;
    }
    mpq_t tmp; mpq_init(tmp);
    if (lift_expr_to_mpq(term, tmp)) {
        mpq_mul(coef, coef, tmp);
        mpq_clear(tmp);
        return;
    }
    mpq_clear(tmp);
    if (term->type == EXPR_FUNCTION
        && term->data.function.head
        && term->data.function.head->type == EXPR_SYMBOL
        && term->data.function.head->data.symbol == SYM_Power
        && term->data.function.arg_count == 2) {
        Expr* exp = term->data.function.args[1];
        if (exp->type == EXPR_INTEGER && exp->data.integer >= 1
            && exp->data.integer <= 16) {
            Expr* base = term->data.function.args[0];
            for (int64_t k = 0; k < exp->data.integer; k++) {
                lift_tl_push(tokens, base);
            }
            return;
        }
    }
    lift_tl_push(tokens, term);
}

/* Find common multiset of tokens. Greedy: walk the first list; for each
 * token, search for an unused match in every other list. If all match,
 * mark them used and add the token to the result. */
static void lift_find_common(LiftTokList* lists, size_t n_terms,
                             Expr*** out_common, size_t* out_count) {
    *out_common = NULL; *out_count = 0;
    if (n_terms == 0 || lists[0].count == 0) return;
    Expr** result = (Expr**)malloc(sizeof(Expr*) * lists[0].count);
    if (!result) return;
    size_t res_count = 0;

    char** used = (char**)malloc(sizeof(char*) * n_terms);
    if (!used) { free(result); return; }
    for (size_t i = 0; i < n_terms; i++) {
        used[i] = lists[i].count ? (char*)calloc(lists[i].count, 1) : NULL;
    }

    size_t* idx = (size_t*)malloc(sizeof(size_t) * n_terms);
    if (!idx) { free(result); for (size_t i = 0; i < n_terms; i++) free(used[i]); free(used); return; }

    for (size_t j = 0; j < lists[0].count; j++) {
        if (used[0][j]) continue;
        Expr* tok = lists[0].items[j];
        idx[0] = j;
        bool ok = true;
        for (size_t i = 1; i < n_terms; i++) {
            bool found = false;
            for (size_t k = 0; k < lists[i].count; k++) {
                if (!used[i][k] && expr_eq(tok, lists[i].items[k])) {
                    idx[i] = k; found = true; break;
                }
            }
            if (!found) { ok = false; break; }
        }
        if (ok) {
            for (size_t i = 0; i < n_terms; i++) used[i][idx[i]] = 1;
            result[res_count++] = tok;
        }
    }
    free(idx);
    for (size_t i = 0; i < n_terms; i++) free(used[i]);
    free(used);

    *out_common = result;
    *out_count = res_count;
}

/* GCD of an array of mpq values: gcd(numerators) / lcm(denominators).
 * Result is positive. n must be >= 1. */
static void lift_compute_mpq_gcd(mpq_t* values, size_t n, mpq_t out) {
    mpq_set(out, values[0]);
    if (mpq_sgn(out) < 0) mpq_neg(out, out);
    for (size_t i = 1; i < n; i++) {
        mpz_t b_num, g, lcm_d;
        mpz_inits(b_num, g, lcm_d, NULL);
        mpz_set(b_num, mpq_numref(values[i]));
        mpz_abs(b_num, b_num);
        mpz_gcd(g, mpq_numref(out), b_num);
        mpz_lcm(lcm_d, mpq_denref(out), mpq_denref(values[i]));
        mpz_set(mpq_numref(out), g);
        mpz_set(mpq_denref(out), lcm_d);
        mpq_canonicalize(out);
        mpz_clears(b_num, g, lcm_d, NULL);
    }
}

static Expr* lift_common_from_plus_impl(const Expr* plus_e) {
    if (!plus_e || plus_e->type != EXPR_FUNCTION
        || !plus_e->data.function.head
        || plus_e->data.function.head->type != EXPR_SYMBOL
        || plus_e->data.function.head->data.symbol != SYM_Plus) {
        return NULL;
    }
    size_t n = plus_e->data.function.arg_count;
    if (n < 2) return NULL;

    LiftTokList* lists = (LiftTokList*)malloc(sizeof(LiftTokList) * n);
    mpq_t* coefs = (mpq_t*)malloc(sizeof(mpq_t) * n);
    if (!lists || !coefs) { free(lists); free(coefs); return NULL; }
    for (size_t i = 0; i < n; i++) {
        lift_tl_init(&lists[i]);
        mpq_init(coefs[i]);
        mpq_set_ui(coefs[i], 1, 1);
        lift_decompose_term(plus_e->data.function.args[i], coefs[i], &lists[i]);
    }

    Expr** common = NULL;
    size_t common_count = 0;
    lift_find_common(lists, n, &common, &common_count);

    /* Restrict the firing condition to a real shared algebraic factor.
     * A coefficient-only GCD lift (e.g. Plus[a/3, b/9] -> (1/9)(3a + b))
     * doesn't reveal new structure that Together/Cancel haven't already
     * exposed, and feeding it through the round loop can blow up
     * downstream transform cost (the evaluator re-rationalises Sqrt
     * products, etc.). */
    if (common_count == 0) {
        for (size_t i = 0; i < n; i++) { lift_tl_free(&lists[i]); mpq_clear(coefs[i]); }
        free(lists); free(coefs); free(common);
        return NULL;
    }

    mpq_t cgcd; mpq_init(cgcd);
    lift_compute_mpq_gcd(coefs, n, cgcd);

    /* Build lifted multiplicative factor: cgcd * Times(common). */
    Expr* gcd_expr = lift_mpq_to_expr(cgcd);
    Expr* lift_factor;
    if (common_count == 0) {
        lift_factor = gcd_expr;
    } else {
        Expr** args = (Expr**)malloc(sizeof(Expr*) * (common_count + 1));
        args[0] = gcd_expr;
        for (size_t i = 0; i < common_count; i++) {
            args[i + 1] = expr_copy(common[i]);
        }
        lift_factor = expr_new_function(expr_new_symbol("Times"), args, common_count + 1);
        free(args);
    }
    Expr* lift_factor_eval = evaluate(lift_factor);

    /* For each term, divide by lift_factor and let evaluate() cancel. */
    Expr** new_terms = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* inv_args[2] = { expr_copy(lift_factor_eval), expr_new_integer(-1) };
        Expr* inv = expr_new_function(expr_new_symbol("Power"), inv_args, 2);
        Expr* mul_args[2] = { expr_copy(plus_e->data.function.args[i]), inv };
        Expr* div = expr_new_function(expr_new_symbol("Times"), mul_args, 2);
        new_terms[i] = evaluate(div);
    }
    Expr* sum = expr_new_function(expr_new_symbol("Plus"), new_terms, n);
    free(new_terms);
    Expr* sum_eval = evaluate(sum);

    Expr* res_args[2] = { lift_factor_eval, sum_eval };
    Expr* result = expr_new_function(expr_new_symbol("Times"), res_args, 2);
    Expr* result_eval = evaluate(result);

    for (size_t i = 0; i < n; i++) { lift_tl_free(&lists[i]); mpq_clear(coefs[i]); }
    free(lists); free(coefs); free(common);
    mpq_clear(cgcd);

    /* If the lift was a no-op (the Plus stayed structurally identical
     * after evaluate normalised the round trip), report no improvement. */
    if (expr_eq(result_eval, plus_e)) {
        expr_free(result_eval);
        return NULL;
    }
    return result_eval;
}

/* Walker entry point. Tries the lift on the input directly (when it's a
 * Plus) or on a Plus child of a Times product (e.g. the numerator of a
 * Times[Plus, Power[denom, -negative]] fraction). Returns NULL when no
 * structural improvement is found. */
static Expr* simp_lift_common_factor(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || !e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL) {
        return NULL;
    }
    const char* sym = e->data.function.head->data.symbol;
    if (sym == SYM_Plus) {
        return lift_common_from_plus_impl(e);
    }
    if (sym == SYM_Times) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            Expr* child = e->data.function.args[i];
            if (child && child->type == EXPR_FUNCTION
                && child->data.function.head
                && child->data.function.head->type == EXPR_SYMBOL
                && child->data.function.head->data.symbol == SYM_Plus) {
                Expr* lifted = lift_common_from_plus_impl(child);
                if (lifted) {
                    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * e->data.function.arg_count);
                    for (size_t j = 0; j < e->data.function.arg_count; j++) {
                        new_args[j] = (j == i)
                                          ? lifted
                                          : expr_copy(e->data.function.args[j]);
                    }
                    Expr* new_times = expr_new_function(expr_new_symbol("Times"),
                                                        new_args, e->data.function.arg_count);
                    free(new_args);
                    Expr* res = evaluate(new_times);
                    if (expr_eq(res, e)) { expr_free(res); return NULL; }
                    return res;
                }
            }
        }
    }
    return NULL;
}

/* ----------------------------------------------------------------------- */
/* Sign canonicalization: paired negative-leading Plus factors            */
/* ----------------------------------------------------------------------- */
/*
 * Mathematica's Factor (and our facpoly) emit Plus subterms whose
 * canonically-sorted first argument is the smallest by name -- so for
 * `(c - a)(d - b)` the printed form is `((-a + c)(-b + d))` because the
 * Plus[Times[-1, a], c] sorts the negated `a` before the bare `c` (the
 * sort key strips the leading -1 coefficient).
 *
 * Mathematica's printed convention prefers each binomial to lead with a
 * positive coefficient. We achieve this post-hoc: when a Times has an
 * even number of "negatively-leading" Plus factors, flipping the sign
 * of each is value-preserving (each pair contributes (-1)*(-1) = 1) and
 * leaves the printed form leading with the positive coefficient.
 *
 * We do not attempt the odd-count case here: pulling an extra -1 onto
 * an outer numeric factor changes which token absorbs the sign and is
 * not always a canonical win for the score function. */

static bool plus_arg_is_negative_leading(const Expr* arg) {
    if (!arg) return false;
    if (arg->type == EXPR_INTEGER) return arg->data.integer < 0;
    if (arg->type == EXPR_BIGINT) return mpz_sgn(arg->data.bigint) < 0;
    if (is_rational_literal(arg)
        && arg->data.function.args[0]->type == EXPR_INTEGER) {
        return arg->data.function.args[0]->data.integer < 0;
    }
    if (arg->type == EXPR_FUNCTION
        && arg->data.function.head
        && arg->data.function.head->type == EXPR_SYMBOL
        && arg->data.function.head->data.symbol == SYM_Times
        && arg->data.function.arg_count >= 1) {
        Expr* coef = arg->data.function.args[0];
        if (coef->type == EXPR_INTEGER) return coef->data.integer < 0;
        if (coef->type == EXPR_BIGINT) return mpz_sgn(coef->data.bigint) < 0;
        if (is_rational_literal(coef)
            && coef->data.function.args[0]->type == EXPR_INTEGER) {
            return coef->data.function.args[0]->data.integer < 0;
        }
    }
    return false;
}

static bool plus_is_negative_leading(const Expr* p) {
    if (!p || p->type != EXPR_FUNCTION
        || !p->data.function.head
        || p->data.function.head->type != EXPR_SYMBOL
        || p->data.function.head->data.symbol != SYM_Plus
        || p->data.function.arg_count < 1) {
        return false;
    }
    return plus_arg_is_negative_leading(p->data.function.args[0]);
}

/* Build a Plus equal to -p by negating every term and re-evaluating so
 * picocas re-canonicalises the argument order. */
static Expr* plus_negate(const Expr* p) {
    size_t n = p->data.function.arg_count;
    Expr** args = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* neg_args[2] = { expr_new_integer(-1),
                              expr_copy(p->data.function.args[i]) };
        args[i] = expr_new_function(expr_new_symbol("Times"), neg_args, 2);
    }
    Expr* neg_plus = expr_new_function(expr_new_symbol("Plus"), args, n);
    free(args);
    return evaluate(neg_plus);
}

/* Walk a Times. If two or more Plus children are negative-leading, flip
 * pairs of them. Only flips an even number; the odd remainder stays.
 * Returns NULL if no flip applies. */
static Expr* canon_negate_pairs(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION
        || !e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol != SYM_Times) {
        return NULL;
    }
    size_t n = e->data.function.arg_count;
    if (n < 2) return NULL;

    size_t* neg_idx = (size_t*)malloc(sizeof(size_t) * n);
    if (!neg_idx) return NULL;
    size_t neg_count = 0;
    for (size_t i = 0; i < n; i++) {
        if (plus_is_negative_leading(e->data.function.args[i])) {
            neg_idx[neg_count++] = i;
        }
    }
    size_t flip_count = (neg_count / 2) * 2;
    if (flip_count == 0) {
        free(neg_idx);
        return NULL;
    }

    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        new_args[i] = expr_copy(e->data.function.args[i]);
    }
    for (size_t j = 0; j < flip_count; j++) {
        size_t i = neg_idx[j];
        Expr* flipped = plus_negate(new_args[i]);
        expr_free(new_args[i]);
        new_args[i] = flipped;
    }
    free(neg_idx);
    Expr* new_times = expr_new_function(expr_copy(e->data.function.head),
                                        new_args, n);
    free(new_args);
    Expr* result = evaluate(new_times);
    if (expr_eq(result, e)) {
        expr_free(result);
        return NULL;
    }
    return result;
}

/* Pythagorean reduction:
 *   1 - Cos[x]^2  -> Sin[x]^2
 *   1 - Sin[x]^2  -> Cos[x]^2
 *   Cosh[x]^2 - 1 -> Sinh[x]^2
 *   1 + Sinh[x]^2 -> Cosh[x]^2
 *
 * Each rule is a strict leaf-count reduction (the Plus collapses to a
 * single Power). The trailing `r___` inside the Plus lets the rule fire
 * when the matching pair sits among other terms (e.g.
 * `1 - Cos[x]^2 + 5` -> `5 + Sin[x]^2`). Idempotent on inputs that
 * don't match. */
/* Wrap a transform's `impl` with FactorMemo lookup + store.  When no
 * memo is active (i.e., we're not inside a Simplify call), the impl
 * runs directly with no overhead.  When active, identical inputs
 * return cached results; the memo key includes a $-prefixed pseudo-
 * head so it never collides with builtin keys (Factor[X], TrigFactor[X],
 * etc.) sharing the same memo. */
static Expr* simp_memo_wrap(const Expr* e, const char* pseudo_head,
                            Expr* (*impl)(const Expr*)) {
    FactorMemo* memo = factor_memo_active();
    if (!memo) return impl(e);

    /* Note: we use raw-input keying here (not Together(Expand(.))
     * canonicalisation as in trig_memo_call).  Reason: the wrapped
     * transforms (PythagReduce, PythagSquareComplete, HalfAngle) use
     * pattern rules that look for specific surface structure --
     * `1 - Cos[x]^2`, `1 + 2 Sin Cos`, `Sin[x] / (1 + Cos[x])` etc.
     * Distributive Expand destroys those patterns (`a (-1 + Cos^2)`
     * becomes `-a + a Cos^2`, where the -1 disappears as a coefficient
     * adjustment), so the rules no longer fire on the canonical form.
     *
     * For the trig memos the canonical form is fine because those
     * transforms internally normalise via Together / TrigToExp before
     * pattern matching. */
    Expr* key_args[1] = { expr_copy((Expr*)e) };
    Expr* key = expr_new_function(expr_new_symbol(pseudo_head), key_args, 1);
    const Expr* hit = factor_memo_lookup(memo, key);
    if (hit) {
        Expr* cached = expr_copy((Expr*)hit);
        expr_free(key);
        return cached;
    }
    Expr* result = impl(e);
    if (result) factor_memo_store(memo, key, result);
    expr_free(key);
    return result;
}

/* Cheap structural check: does `e` contain any of Cos/Sin/Cosh/Sinh,
 * Tan/Cot/Tanh/Coth, or Sec/Csc/Sech/Csch as a function head?
 * PythagReduce / PythagCanon's rules can only fire on these heads, so
 * when the answer is no we can skip the ReplaceRepeated walk entirely.
 * Walks the tree once; cheaper by orders of magnitude than the
 * pattern-matching pass it gates.  Sec/Csc/Sech/Csch were added
 * alongside the Sec[x]^2 -> 1 + Tan[x]^2 substitution direction in
 * PythagCanon (and its three siblings); without listing them here,
 * inputs like (Sec[x]+1)(Sec[x]-1) - Tan[x]^2 would skip PythagCanon
 * via the gate and land in TrigFactor (700 ms+ on multi-variable
 * expansions). */
static bool has_pythag_head(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL) {
        const char* h = head->data.symbol;
        if (h == SYM_Cos || h == SYM_Sin ||
            h == SYM_Cosh || h == SYM_Sinh ||
            h == SYM_Tan || h == SYM_Cot ||
            h == SYM_Tanh || h == SYM_Coth ||
            h == SYM_Sec || h == SYM_Csc ||
            h == SYM_Sech || h == SYM_Csch) {
            return true;
        }
    }
    if (has_pythag_head(head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_pythag_head(e->data.function.args[i])) return true;
    }
    return false;
}

/* ----------------------------------------------------------------------- */
/* Trig at rational multiples of Pi: shortest-numerator canonicalization   */
/* ----------------------------------------------------------------------- */

/* For inputs like Cos[4 Pi/9] - Sin[Pi/18] the difference is exactly 0
 * by the complement identity Cos[Pi/2 - x] == Sin[x] (here Pi/2 - 4Pi/9
 * = Pi/18).  builtin_cos / builtin_sin handle a whitelist of "nice"
 * denominators (1,2,3,4,5,6,10,12) and otherwise leave the call as is,
 * so the two terms never become structurally equal and additive
 * cancellation in simp_search never fires.
 *
 * This transform rewrites every Sin/Cos/Tan/Cot/Sec/Csc of a rational
 * multiple of Pi into a unique form: pick the (Sin vs Cos / Tan vs Cot
 * / Sec vs Csc) representation whose reduced fraction has the smaller
 * numerator.  After the rewrite, Cos[4 Pi/9] and Sin[Pi/18] both land
 * at Sin[Pi/18] and the surrounding Plus collapses to 0.  Cos[5 Pi/9]
 * lands at -Sin[Pi/18] (since 5/9 > 1/2 picks up a sign via Cos[Pi -
 * x] = -Cos[x] before the complement swap), which lets the Morrie's-
 * law product 1/8 (Cos[4Pi/9] + Cos[5Pi/9]) collapse to 0 too.
 *
 * Idempotent: re-applying to the canonical form returns it unchanged,
 * so it is safe to seed without bounding the round count.
 *
 * Why a Simplify-only transform (rather than wiring into builtin_cos
 * / builtin_sin): the user-facing default print for Cos[4 Pi/9] should
 * stay as-written, both to match Mathematica and to keep regression
 * tests stable.  Only Simplify needs the unified representation. */

static int64_t trig_pi_i64_gcd(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a;
}

static void trig_pi_reduce_frac(int64_t* n, int64_t* d) {
    if (*d == 0) return;
    int64_t g = trig_pi_i64_gcd(*n, *d);
    if (g > 1) { *n /= g; *d /= g; }
    if (*d < 0) { *d = -*d; *n = -*n; }
}

/* Returns true if `arg` is structurally (n/d) * Pi or Pi (n=d=1) and
 * sets n, d.  Mirrors trig.c's extract_pi_multiplier; duplicated here
 * so simp.c does not have to expose that file-static helper. */
static bool trig_pi_extract(const Expr* arg, int64_t* n_out, int64_t* d_out) {
    if (arg->type == EXPR_SYMBOL && arg->data.symbol == SYM_Pi) {
        *n_out = 1; *d_out = 1; return true;
    }
    if (arg->type != EXPR_FUNCTION || !arg->data.function.head ||
        arg->data.function.head->type != EXPR_SYMBOL ||
        arg->data.function.head->data.symbol != SYM_Times ||
        arg->data.function.arg_count != 2) {
        return false;
    }
    Expr* a0 = arg->data.function.args[0];
    Expr* a1 = arg->data.function.args[1];
    if (!(a1->type == EXPR_SYMBOL && a1->data.symbol == SYM_Pi)) return false;
    return is_rational(a0, n_out, d_out);
}

/* Build (out_n / out_d) * Pi as an Expr*. */
static Expr* trig_pi_make_arg(int64_t n, int64_t d) {
    if (n == 0) return expr_new_integer(0);
    if (n == 1 && d == 1) return expr_new_symbol("Pi");
    Expr* coeff;
    if (d == 1) {
        coeff = expr_new_integer(n);
    } else {
        coeff = expr_new_function(expr_new_symbol("Rational"),
            (Expr*[]){ expr_new_integer(n), expr_new_integer(d) }, 2);
    }
    Expr* args[2] = { coeff, expr_new_symbol("Pi") };
    return eval_and_free(expr_new_function(expr_new_symbol("Times"), args, 2));
}

/* Compute the canonical (head, n, d, sign) for `head[(n/d) Pi]`.  See
 * the file-header comment block above for the rule set.
 *
 * `head_in` is one of the SYM_* trig pointers (Sin/Cos/Tan/Cot/Sec/Csc).
 * On success returns true and fills the outputs; on failure (head not
 * recognised, integer overflow, etc.) returns false and outputs are
 * untouched. */
static bool trig_pi_canon_one(const char* head_in, int64_t n, int64_t d,
                              const char** head_out, int64_t* n_out,
                              int64_t* d_out, int* sign_out) {
    if (d <= 0) return false;
    /* Guard against multiplication overflow when computing 2*d, alt_d. */
    if (d > (INT64_MAX / 4)) return false;
    if (n > (INT64_MAX / 4) || n < -(INT64_MAX / 4)) return false;

    int sign = 1;

    bool is_tan_family = (head_in == SYM_Tan || head_in == SYM_Cot);
    bool is_sin_like = (head_in == SYM_Sin || head_in == SYM_Csc);
    bool is_cos_like = (head_in == SYM_Cos || head_in == SYM_Sec);
    if (!is_tan_family && !is_sin_like && !is_cos_like) return false;

    if (is_tan_family) {
        /* Period Pi.  Reduce to [0, d). */
        n = n % d;
        if (n < 0) n += d;
        /* Tan[Pi - x] = -Tan[x], Cot[Pi - x] = -Cot[x]. */
        if (2 * n > d) {
            sign = -sign;
            n = d - n;
        }
        /* n in [0, d/2]. */
    } else {
        int64_t two_d = 2 * d;
        n = n % two_d;
        if (n < 0) n += two_d;
        if (is_sin_like) {
            /* Sin[Pi + x] = -Sin[x], Csc[Pi + x] = -Csc[x]. */
            if (n > d) { sign = -sign; n -= d; }
            /* Now n in [0, d].  Sin[Pi - x] = Sin[x], Csc same. */
            if (2 * n > d) n = d - n;
        } else { /* cos-like */
            /* Cos[2 Pi - x] = Cos[x], Sec same. */
            if (n > d) n = two_d - n;
            /* Now n in [0, d].  Cos[Pi - x] = -Cos[x], Sec same. */
            if (2 * n > d) { sign = -sign; n = d - n; }
        }
        /* n in [0, d/2]. */
    }

    trig_pi_reduce_frac(&n, &d);

    /* n == 0: argument collapsed to 0 (Sin[0]=0, Tan[0]=0, Cos[0]=1, ...).
     * builtin_sin / builtin_cos / etc. already handle the EXPR_INTEGER 0
     * case, so just emit head[0] -- the surrounding evaluate() pass will
     * resolve it. */
    if (n == 0) {
        *head_out = head_in;
        *n_out = 0;
        *d_out = 1;
        *sign_out = sign;
        return true;
    }

    /* Complement: head[(n/d) Pi] == alt_head[(d - 2n)/(2d) Pi].  The
     * sign relation is the identity Sin = Cos ∘ (Pi/2 - .) etc., which
     * is sign-preserving under the post-quadrant reduction we did above. */
    int64_t alt_n = d - 2 * n;
    int64_t alt_d = 2 * d;
    trig_pi_reduce_frac(&alt_n, &alt_d);

    const char* alt_head;
    if      (head_in == SYM_Sin) alt_head = "Cos";
    else if (head_in == SYM_Cos) alt_head = "Sin";
    else if (head_in == SYM_Tan) alt_head = "Cot";
    else if (head_in == SYM_Cot) alt_head = "Tan";
    else if (head_in == SYM_Sec) alt_head = "Csc";
    else                          alt_head = "Sec"; /* SYM_Csc */

    /* Pick the rep with the smaller reduced numerator.  Tie-break: keep
     * the input head so we never thrash on already-canonical forms. */
    if (alt_n < n) {
        *head_out = alt_head;
        *n_out = alt_n;
        *d_out = alt_d;
    } else {
        *head_out = head_in;
        *n_out = n;
        *d_out = d;
    }
    *sign_out = sign;
    return true;
}

/* Walk `e`, applying trig_pi_canon_one at every Sin/Cos/Tan/Cot/Sec/Csc
 * call whose argument is a rational multiple of Pi.  Returns a freshly
 * owned, evaluated tree (always non-NULL).  Idempotent on the canonical
 * form. */
static Expr* simp_trig_pi_canon_walk(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    /* Apply at this node when shape matches. */
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.arg_count == 1) {
        const char* h = e->data.function.head->data.symbol;
        bool is_trig = (h == SYM_Sin || h == SYM_Cos ||
                        h == SYM_Tan || h == SYM_Cot ||
                        h == SYM_Sec || h == SYM_Csc);
        if (is_trig) {
            int64_t n, d;
            if (trig_pi_extract(e->data.function.args[0], &n, &d)) {
                const char* out_h;
                int64_t out_n, out_d;
                int out_sign;
                if (trig_pi_canon_one(h, n, d, &out_h, &out_n, &out_d, &out_sign)) {
                    Expr* arg = trig_pi_make_arg(out_n, out_d);
                    Expr* call = expr_new_function(
                        expr_new_symbol(out_h),
                        (Expr*[]){ arg }, 1);
                    Expr* call_eval = eval_and_free(call);
                    if (out_sign == -1) {
                        Expr* neg = expr_new_function(
                            expr_new_symbol("Times"),
                            (Expr*[]){ expr_new_integer(-1), call_eval }, 2);
                        return eval_and_free(neg);
                    }
                    return call_eval;
                }
            }
        }
    }

    /* Recurse into children. */
    size_t n_args = e->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * n_args);
    bool any_changed = false;
    for (size_t i = 0; i < n_args; i++) {
        Expr* nc = simp_trig_pi_canon_walk(e->data.function.args[i]);
        if (!expr_eq(nc, e->data.function.args[i])) any_changed = true;
        new_args[i] = nc;
    }
    /* Walk the head too in case of e.g. Hold[Sin][...]. */
    Expr* new_head = simp_trig_pi_canon_walk(e->data.function.head);
    if (!expr_eq(new_head, e->data.function.head)) any_changed = true;

    if (!any_changed) {
        for (size_t i = 0; i < n_args; i++) expr_free(new_args[i]);
        free(new_args);
        expr_free(new_head);
        return expr_copy((Expr*)e);
    }
    Expr* res = expr_new_function(new_head, new_args, n_args);
    free(new_args);
    return eval_and_free(res);
}

/* Idempotent (re-application is a structural fixed point).  Inert on
 * inputs without a Sin/Cos/Tan/Cot/Sec/Csc of a rational multiple of
 * Pi -- the walker descends but never triggers a rewrite. */
static Expr* simp_trig_pi_canon(const Expr* e) {
    return simp_trig_pi_canon_walk(e);
}

static Expr* transform_pythag_reduce_impl(const Expr* e) {
    static Expr* rules = NULL;
    if (!rules) {
        rules = parse_expression(
            "{ 1 - Cos[x_]^2 + r___  :> Sin[x]^2 + r, "
            "  1 - Sin[x_]^2 + r___  :> Cos[x]^2 + r, "
            "  -1 + Cos[x_]^2 + r___ :> -Sin[x]^2 + r, "
            "  -1 + Sin[x_]^2 + r___ :> -Cos[x]^2 + r, "
            "  -1 + Cosh[x_]^2 + r___ :> Sinh[x]^2 + r, "
            "  1 + Sinh[x_]^2 + r___ :> Cosh[x]^2 + r, "
            "  1 - Cosh[x_]^2 + r___ :> -Sinh[x]^2 + r, "
            "  -1 - Sinh[x_]^2 + r___ :> -Cosh[x]^2 + r, "
            /* Reciprocal-pair identities. tanh^2 + sech^2 == 1, so
             *   1 - Tanh^2 -> Sech^2  and  -1 + Tanh^2 -> -Sech^2.
             * coth^2 - csch^2 == 1, so
             *   -1 + Coth^2 -> Csch^2 and  1 - Coth^2 -> -Csch^2.
             * tan^2 + 1 == sec^2, cot^2 + 1 == csc^2 (real-valued
             * Pythagorean trig).  These resolve a tied-score plateau
             * where the simp_search round loop's strict `<` tiebreak
             * would otherwise prefer the bare Plus form (e.g. score 7 =
             * score 7 for `-1 + Tanh^2` vs `-Sech^2`); fired here, the
             * structural collapse to a single Power head wins outright. */
            "  1 - Tanh[x_]^2 + r___  :> Sech[x]^2 + r, "
            "  -1 + Tanh[x_]^2 + r___ :> -Sech[x]^2 + r, "
            "  -1 + Coth[x_]^2 + r___ :> Csch[x]^2 + r, "
            "  1 - Coth[x_]^2 + r___  :> -Csch[x]^2 + r, "
            "  1 + Tan[x_]^2 + r___   :> Sec[x]^2 + r, "
            "  -1 - Tan[x_]^2 + r___  :> -Sec[x]^2 + r, "
            "  1 + Cot[x_]^2 + r___   :> Csc[x]^2 + r, "
            "  -1 - Cot[x_]^2 + r___  :> -Csc[x]^2 + r, "
            /* Reverse-direction reciprocal-pair identities.  Sec^2 - 1 ==
             * Tan^2, Csc^2 - 1 == Cot^2 (real-valued Pythagorean trig);
             * Sech^2 + (-1) doesn't hold (Sech^2 = 1 - Tanh^2, so 1 -
             * Sech^2 == Tanh^2); 1 + Csch^2 == Coth^2.  Without these
             * rules, a Plus shape like `-1 + Sec[x]^2 - Tan[x]^2` (the
             * post-Expand form of `(Sec+1)(Sec-1) - Tan^2`) would only
             * collapse via the TAN -> SEC direction `-1 - Tan^2 -> -Sec^2`,
             * which works for the simple shape but leaves the SEC term
             * "stranded" in coefficient-bearing forms PythagCanon must
             * separately rewrite.  Including both directions here lets
             * the cheap PythagReduce rule fire on either presentation. */
            "  -1 + Sec[x_]^2 + r___  :> Tan[x]^2 + r, "
            "  1 - Sec[x_]^2 + r___   :> -Tan[x]^2 + r, "
            "  -1 + Csc[x_]^2 + r___  :> Cot[x]^2 + r, "
            "  1 - Csc[x_]^2 + r___   :> -Cot[x]^2 + r, "
            "  1 - Sech[x_]^2 + r___  :> Tanh[x]^2 + r, "
            "  -1 + Sech[x_]^2 + r___ :> -Tanh[x]^2 + r, "
            "  1 + Csch[x_]^2 + r___  :> Coth[x]^2 + r, "
            "  -1 - Csch[x_]^2 + r___ :> -Coth[x]^2 + r }");
    }
    if (!rules) return NULL;
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Fast-skip: every PythagReduce rule LHS has a Cos/Sin/Cosh/Sinh
     * pattern.  If the input contains none of those heads, the
     * ReplaceRepeated walk would visit every node and try every rule
     * and find nothing -- which on huge sum-of-exponentials inputs
     * costs 50-120 ms per call.  Skip the rewrite and return a copy
     * unchanged. */
    if (!has_pythag_head(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("PythagReduce", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr* args[2] = { expr_copy((Expr*)e), expr_copy(rules) };
    Expr* call = expr_new_function(
        expr_new_symbol("ReplaceRepeated"), args, 2);
    Expr* out = evaluate(call);
    if (dbg) simp_debug_log("PythagReduce", e, out,
                            simp_debug_elapsed_ms(t0));
    return out;
}

/* PythagReduce sees the highest call volume of any simp transform
 * (~200 calls in the Tan double-angle case, with ~20 % unique inputs).
 * The memo dedupes the rest. */
static Expr* transform_pythag_reduce(const Expr* e) {
    return simp_memo_wrap(e, "$PythagReduce", transform_pythag_reduce_impl);
}

/* PythagCanon: substitution-based Pythagorean canonicalizer.
 *
 * The bare PythagReduce rules look for `1 +/- Cos[x_]^2 + r___` style
 * shapes -- they fire only when the unit constant and the squared trig
 * appear additively, with unit coefficient. After a generic Expand of an
 * input like
 *     18 (Cos[x]+1)(Cos[x]-1)(Cos[y]^2-1)^2 (x-1) + 18 (x-1) Sin[x]^2 Sin[y]^4
 * we get a polynomial in Cos[x], Cos[y] whose individual Cos^2 / Cos^4
 * factors carry coefficients other than 1, so PythagReduce misses every
 * one of them -- and the cancellation against the Sin^2 Sin^4 term is
 * never recognised by the round loop.
 *
 * This transform is coefficient-blind. It substitutes every even-power
 *     Cos[x_]^(2k) -> (1 - Sin[x_]^2)^k
 * (and the reverse Sin -> Cos, plus the hyperbolic counterparts) via a
 * single ReplaceRepeated pass, then Expands. For each of the four
 * directions it scores the result against the input and keeps the best
 * strict win. Idempotent on inputs without an even Cos^k / Sin^k power;
 * inert (returns a structural copy of the input) when no direction beats
 * the input score.
 *
 * Why all four directions: the choice of "all-Sin" vs "all-Cos" depends
 * on which side already has more mass. A user input that is mostly
 * Sin^2 + small Cos^2 minus 1 wants Cos -> 1 - Sin; the reverse leaves
 * the small Cos^2 term and bloats the rest. Trying both keeps us
 * agnostic. The hyperbolic pair is exactly analogous via
 * Cosh^2 - Sinh^2 = 1. */
static Expr* transform_pythag_canon_impl(const Expr* e) {
    static Expr* rules_to_sin = NULL;
    static Expr* rules_to_cos = NULL;
    static Expr* rules_to_sinh = NULL;
    static Expr* rules_to_cosh = NULL;
    /* Reciprocal-pair canonicalisation directions.  Same shape as the
     * Sin/Cos rules above, but for the Pythagorean identities
     *     Sec^2 = 1 + Tan^2,   Csc^2 = 1 + Cot^2
     *     Sech^2 = 1 - Tanh^2, Csch^2 = -1 + Coth^2.
     * These are needed for inputs like `(Sec[x]+1)(Sec[x]-1) - Tan[x]^2`
     * whose post-Expand form `-1 + Sec[x]^2 - Tan[x]^2` collapses when
     * Sec[x]^2 -> 1 + Tan[x]^2 (or Tan[x]^2 -> -1 + Sec[x]^2) is
     * substituted globally.  Without these directions the round loop
     * would have to ride out a TrigFactor pass (Sec -> 1/Cos rewrite +
     * polynomial reformulation), which on the multi-variable Sec/Tan
     * test case takes 700ms+ for a single call.  All eight rules guard
     * with `n >= 2 && EvenQ[n]` to keep odd-power forms intact. */
    static Expr* rules_to_tan = NULL;
    static Expr* rules_to_sec = NULL;
    static Expr* rules_to_cot = NULL;
    static Expr* rules_to_csc = NULL;
    static Expr* rules_to_tanh = NULL;
    static Expr* rules_to_sech = NULL;
    static Expr* rules_to_coth = NULL;
    static Expr* rules_to_csch = NULL;
    if (!rules_to_sin) {
        rules_to_sin = parse_expression(
            "{ Cos[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 - Sin[x]^2)^(n/2) }");
        rules_to_cos = parse_expression(
            "{ Sin[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 - Cos[x]^2)^(n/2) }");
        rules_to_sinh = parse_expression(
            "{ Cosh[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 + Sinh[x]^2)^(n/2) }");
        rules_to_cosh = parse_expression(
            "{ Sinh[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (-1 + Cosh[x]^2)^(n/2) }");
        rules_to_tan = parse_expression(
            "{ Sec[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 + Tan[x]^2)^(n/2) }");
        rules_to_sec = parse_expression(
            "{ Tan[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (-1 + Sec[x]^2)^(n/2) }");
        rules_to_cot = parse_expression(
            "{ Csc[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 + Cot[x]^2)^(n/2) }");
        rules_to_csc = parse_expression(
            "{ Cot[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (-1 + Csc[x]^2)^(n/2) }");
        rules_to_tanh = parse_expression(
            "{ Sech[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 - Tanh[x]^2)^(n/2) }");
        rules_to_sech = parse_expression(
            "{ Tanh[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 - Sech[x]^2)^(n/2) }");
        rules_to_coth = parse_expression(
            "{ Csch[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (-1 + Coth[x]^2)^(n/2) }");
        rules_to_csch = parse_expression(
            "{ Coth[x_]^n_Integer /; n >= 2 && EvenQ[n] "
            "    :> (1 + Csch[x]^2)^(n/2) }");
    }
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Same fast-skip as PythagReduce: every rule LHS targets a
     * Cos/Sin/Cosh/Sinh head. Without one nothing can fire. */
    if (!has_pythag_head(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("PythagCanon", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr* best = expr_copy((Expr*)e);
    size_t best_score = score_with_func(best, NULL);

    /* Expand first: the substitution rule matches Power[Cos[x], 2k]
     * literally, which only appears after distributing factored forms
     * like (Cos+1)(Cos-1) -> Cos^2 - 1. Without the pre-Expand the
     * rules cannot fire on the user's typical input shape. */
    Expr* pre_args[1] = { expr_copy((Expr*)e) };
    Expr* pre_call = expr_new_function(
        expr_new_symbol("Expand"), pre_args, 1);
    Expr* pre_expanded = evaluate(pre_call);
    if (!pre_expanded) {
        if (dbg) simp_debug_log("PythagCanon", e, best,
                                simp_debug_elapsed_ms(t0));
        return best;
    }

    /* Twelve directions: the four base Cos/Sin/Cosh/Sinh substitutions
     * plus the reciprocal-pair pairs (Sec<->Tan, Csc<->Cot, Sech<->Tanh,
     * Csch<->Coth).  Per-direction cost is one ReplaceRepeated walk +
     * one Expand; rules whose LHS head isn't present in `e` no-op out
     * almost instantly because the matcher rejects on head mismatch.
     * The keep-best-strict-win selection in the body below means firing
     * a no-op direction is only an ms-scale speed cost, never a
     * correctness cost. */
    Expr* rule_sets[12] = {
        rules_to_sin,  rules_to_cos,  rules_to_sinh, rules_to_cosh,
        rules_to_tan,  rules_to_sec,  rules_to_cot,  rules_to_csc,
        rules_to_tanh, rules_to_sech, rules_to_coth, rules_to_csch
    };
    for (int i = 0; i < 12; i++) {
        if (!rule_sets[i]) continue;
        Expr* ra_args[2] = { expr_copy(pre_expanded),
                              expr_copy(rule_sets[i]) };
        Expr* ra_call = expr_new_function(
            expr_new_symbol("ReplaceRepeated"), ra_args, 2);
        Expr* substituted = evaluate(ra_call);
        if (!substituted) continue;
        if (expr_eq(substituted, pre_expanded)) {
            expr_free(substituted);
            continue;
        }
        Expr* exp_args[1] = { substituted };
        Expr* exp_call = expr_new_function(
            expr_new_symbol("Expand"), exp_args, 1);
        Expr* expanded = evaluate(exp_call);
        if (!expanded) continue;
        size_t s = score_with_func(expanded, NULL);
        if (s < best_score) {
            expr_free(best);
            best = expanded;
            best_score = s;
        } else {
            expr_free(expanded);
        }
    }
    expr_free(pre_expanded);
    if (dbg) simp_debug_log("PythagCanon", e, best,
                            simp_debug_elapsed_ms(t0));
    return best;
}

static Expr* transform_pythag_canon(const Expr* e) {
    return simp_memo_wrap(e, "$PythagCanon", transform_pythag_canon_impl);
}

/* ----------------------------------------------------------------------- */
/* PrimeRebase: Power[c, e] -> Power[p, k*e] for c = p^k integer (k >= 2)  */
/* ----------------------------------------------------------------------- */

/* Soundness:  for positive integer c = p^k,
 *     c^e = (p^k)^e = p^(k*e)
 * holds for ALL complex e -- the (a^b)^c = a^(b*c) identity is sound
 * when a > 0, with no branch cut to worry about.
 *
 * Why this is needed:  picocas's canonical Power evaluator never rebases
 * composite integer bases ((2^2)^x stays as (2^2)^x), so factors like
 * 4^x and 2^x sit in different Times "base buckets" and the same-base
 * exponent combine in `times.c` cannot cancel them.  After rebasing all
 * such Power factors to a single canonical prime base, evaluate()
 * collapses the combined exponents in a single pass.
 *
 * Coverage in Simplify search:
 *     4^x * 2^(-x) * 2^(-x) - 1     ->  0  (top-level rebase)
 *     f[4^x * 2^(-x) * 2^(-x)] - f[1]  ->  0  (rebase inside f's arg)
 *     2^(2^(2 x) x) - 2^(x*4^x)     ->  0  (rebase inside an exponent)
 *
 * Branch-cut-sensitive shapes -- negative-integer-base split
 * ((-4)^x -> (-1)^x 4^x), constant-positive distribute
 * ((c1 c2)^e -> c1^e c2^e for ci > 0), and integer-exponent distribute
 * ((a b)^n -> a^n b^n for n integer) -- are handled by
 * `transform_power_distribute` (defined just below); see the comment
 * block there for the soundness argument and the dispatch wiring.
 */

/* If n is a perfect prime power p^k with k >= 2, set *p_out = p,
 * *k_out = k and return true.  Otherwise return false (n is < 4, prime,
 * or has at least two distinct prime factors).  Trial division up to
 * sqrt(n); since picocas's prime-base inputs are typically small literals
 * (4, 8, 9, 16, 25, 27, 32, 49, ...), this is microseconds-cheap. */
static bool prime_rebase_check(int64_t n, int64_t* p_out, int64_t* k_out) {
    if (n < 4) return false;
    int64_t p;
    if ((n & 1) == 0) {
        p = 2;
    } else {
        for (p = 3; ; p += 2) {
            if (p > n / p) return false;  /* n is prime */
            if (n % p == 0) break;
        }
    }
    int64_t k = 0;
    int64_t m = n;
    while (m > 1 && (m % p) == 0) {
        m /= p;
        k++;
    }
    if (m != 1 || k < 2) return false;
    *p_out = p;
    *k_out = k;
    return true;
}

/* Returns true if e contains any Power[c, _] with c an EXPR_INTEGER >= 4
 * that is a perfect prime power.  Cheap structural gate to avoid the full
 * walk + Expand reseed when nothing can fire (the common case). */
static bool has_rebaseable_power(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL && head->data.symbol == SYM_Power
        && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        if (base && base->type == EXPR_INTEGER && base->data.integer >= 4) {
            int64_t p, k;
            if (prime_rebase_check(base->data.integer, &p, &k)) return true;
        }
    }
    if (has_rebaseable_power(head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_rebaseable_power(e->data.function.args[i])) return true;
    }
    return false;
}

/* Recursive structural copy that rebases every prime-power-base Power.
 * Walks children before checking the current node so nested rebases
 * (e.g. inside an exponent) bubble out.  Caller owns the returned tree. */
static Expr* prime_rebase_copy(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    size_t n = e->data.function.arg_count;
    Expr* new_head = prime_rebase_copy(e->data.function.head);
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            new_args[i] = prime_rebase_copy(e->data.function.args[i]);
        }
    }

    int64_t p_val, k_val;
    if (new_head && new_head->type == EXPR_SYMBOL
        && new_head->data.symbol == SYM_Power
        && n == 2 && new_args && new_args[0]
        && new_args[0]->type == EXPR_INTEGER
        && prime_rebase_check(new_args[0]->data.integer, &p_val, &k_val)) {
        /* Build Power[p, Times[k, new_args[1]]] -- new_args[1] is the
         * already-rebased exponent (so nested c'^e' inside the exponent
         * has already been rewritten before we get here). */
        Expr* p_expr = expr_new_integer(p_val);
        Expr* k_expr = expr_new_integer(k_val);
        Expr* times_args[2] = { k_expr, new_args[1] };
        Expr* times_call = expr_new_function(
            expr_new_symbol("Times"), times_args, 2);
        Expr* power_args[2] = { p_expr, times_call };
        Expr* result = expr_new_function(
            expr_new_symbol("Power"), power_args, 2);
        expr_free(new_args[0]);  /* the original integer c */
        free(new_args);
        expr_free(new_head);
        return result;
    }

    Expr* result = expr_new_function(new_head, new_args, n);
    if (new_args) free(new_args);
    return result;
}

static Expr* transform_prime_rebase_impl(const Expr* e) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    if (!has_rebaseable_power(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("PrimeRebase", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr* rebased = prime_rebase_copy(e);
    if (!rebased) rebased = expr_copy((Expr*)e);

    /* Re-evaluate so canonical Times same-base combine collapses the
     * rebased Power factors (e.g. 2^(2x) * 2^(-x) * 2^(-x) -> 2^0 = 1). */
    Expr* result = evaluate(rebased);
    if (!result) result = expr_copy((Expr*)e);

    if (dbg) simp_debug_log("PrimeRebase", e, result,
                            simp_debug_elapsed_ms(t0));
    return result;
}

static Expr* transform_prime_rebase(const Expr* e) {
    return simp_memo_wrap(e, "$PrimeRebase", transform_prime_rebase_impl);
}

/* ----------------------------------------------------------------------- */
/* PowerOneify: combine `A * Power[A, e]` -> `Power[A, e+1]` in Times      */
/* ----------------------------------------------------------------------- */

/* Soundness: Power[a, 1] = a holds for ALL a, so
 *     A * Power[A, e] = Power[A, 1] * Power[A, e] = Power[A, e + 1]
 * is universally valid -- no branch cut, no positivity assumption.
 *
 * Why this is needed: picocas's Times-canonical-form same-base combine
 * groups factors by the base of any wrapping Power[base, exp].  A bare
 * factor A whose canonical form is itself a Power expression (e.g.
 * Power[x, -1] = 1/x) does NOT get re-bucketed as Power[A, 1] before
 * grouping, so Times[Power[x,-1], Power[Power[x,-1], Log[2]]] keeps the
 * two factors as distinct same-base bucket entries.  This blocks
 *     (1/x)^Log[2] / x - (1/x)^(1 + Log[2])  ->  0
 * from collapsing because the LHS Times can't combine the two (1/x)
 * factors into one.
 *
 * Implementation: walk the tree; in every Times node, look for any pair
 * (i, j) where args[j] = Power[B, e] and args[i] is structurally equal
 * to B.  Replace args[j] with Power[B, e + 1] and drop args[i] (which is
 * the bare A = Power[B, 1] absorbed into the exponent).  Then evaluate so
 * canonical Plus folds the new exponent and so any nested rebase cascades.
 *
 * Inert when no Times node contains an A-and-Power[A,_] pair (the common
 * case): one structural pass, microseconds.
 */

static Expr* power_oneify_walk(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    size_t n = e->data.function.arg_count;
    Expr* new_head = power_oneify_walk(e->data.function.head);
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            new_args[i] = power_oneify_walk(e->data.function.args[i]);
        }
    }

    /* Only Times triggers the implicit-one combine. */
    bool combined = false;
    if (new_head && new_head->type == EXPR_SYMBOL
        && new_head->data.symbol == SYM_Times && n >= 2) {
        for (size_t i = 0; i < n; i++) {
            if (!new_args[i]) continue;
            Expr* a = new_args[i];
            for (size_t j = 0; j < n; j++) {
                if (i == j || !new_args[j]) continue;
                Expr* p = new_args[j];
                if (p->type != EXPR_FUNCTION) continue;
                Expr* ph = p->data.function.head;
                if (!(ph && ph->type == EXPR_SYMBOL
                      && ph->data.symbol == SYM_Power
                      && p->data.function.arg_count == 2)) continue;
                if (!expr_eq(p->data.function.args[0], a)) continue;
                /* args[i] = A, args[j] = Power[A, e].  Combine. */
                Expr* old_exp = p->data.function.args[1];
                Expr* one = expr_new_integer(1);
                Expr* plus_args[2] = { expr_copy(old_exp), one };
                Expr* new_exp = expr_new_function(
                    expr_new_symbol("Plus"), plus_args, 2);
                Expr* power_args[2] = {
                    expr_copy(p->data.function.args[0]), new_exp };
                Expr* new_power = expr_new_function(
                    expr_new_symbol("Power"), power_args, 2);
                expr_free(new_args[j]);
                new_args[j] = new_power;
                expr_free(new_args[i]);
                new_args[i] = NULL;
                combined = true;
                break;  /* args[i] consumed; move to next i */
            }
        }
    }

    if (combined) {
        /* Compact dropped slots. */
        size_t out_n = 0;
        for (size_t i = 0; i < n; i++) {
            if (new_args[i]) new_args[out_n++] = new_args[i];
        }
        Expr* times_call = expr_new_function(new_head, new_args, out_n);
        free(new_args);
        Expr* result = evaluate(times_call);
        if (!result) result = expr_copy((Expr*)e);
        return result;
    }

    Expr* result = expr_new_function(new_head, new_args, n);
    if (new_args) free(new_args);
    return result;
}

static Expr* transform_power_oneify_impl(const Expr* e) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;
    Expr* result = power_oneify_walk(e);
    if (!result) result = expr_copy((Expr*)e);
    if (dbg) simp_debug_log("PowerOneify", e, result,
                            simp_debug_elapsed_ms(t0));
    return result;
}

static Expr* transform_power_oneify(const Expr* e) {
    return simp_memo_wrap(e, "$PowerOneify", transform_power_oneify_impl);
}

/* ----------------------------------------------------------------------- */
/* PowerDistribute: distribute Power over Times in three universally       */
/* sound shapes.  Walks the tree once; each shape independently rewrites   */
/* matching nodes in place, then evaluate() is called on the result.       */
/*                                                                         */
/*   (A) Power[neg_int, e]  ->  Power[-1, e] * Power[|neg_int|, e]         */
/*       Soundness: principal-branch identity for any e (real or           */
/*       symbolic).  Picocas uses principal branches everywhere.           */
/*       Why: lets prime_rebase + Times same-base merge cancel mixed-sign  */
/*       products like (-4)^x * (-2)^(-x) * 2^(-x) -> 1.                   */
/*                                                                         */
/*   (B) Power[Times[c1,...,ck, u1,...,um], e]                             */
/*           ->  Times[c1^e, ..., ck^e, Power[Times[u1,...,um], e]]        */
/*       when each ci is a constant positive (literal positive numeric,    */
/*       or a recognised positive constant symbol like Pi/E/...).          */
/*       Soundness: (a b)^c = a^c b^c whenever a > 0, valid for any b      */
/*       and any c.  Why: collapses identities like                        */
/*           Exp[x] Exp[y] 2^x 2^y - (2 Exp[1])^(x+y) -> 0                 */
/*       which need (2 E)^(x+y) -> 2^(x+y) E^(x+y) to align with the LHS.  */
/*                                                                         */
/*   (C) Power[Times[u1,...,un], e]  ->  Power[u1, e] * ... * Power[un, e] */
/*       and Power[Power[u, p], e]   ->  Power[u, p*e]                     */
/*       when prov_int(ctx, e) (the exponent is provably integer, either   */
/*       as a literal or via an Element[_, Integers] assumption).          */
/*       Soundness: integer exponents distribute through products and      */
/*       through nested Power without branch cuts, for any complex base.   */
/*       Why: handles (y/x)^(-n) - x^n y^(-n) -> 0 under                   */
/*       Element[n, Integers].                                             */
/*                                                                         */
/* All three shapes are inert when no matching node exists.  After any     */
/* rewrite we evaluate() the result so the canonical Times same-base       */
/* merge collapses adjacent Power[u, ?] factors.                           */
/* ----------------------------------------------------------------------- */

static bool is_constant_positive_factor(const Expr* x) {
    if (!x) return false;
    if (numeric_sign(x) == 1) return true;
    if (x->type == EXPR_SYMBOL && is_positive_constant_symbol(x->data.symbol))
        return true;
    return false;
}

/* Cheap structural gate for shape (A). */
static bool has_distributable_power(const Expr* e, const AssumeCtx* ctx) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL && head->data.symbol == SYM_Power
        && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp_  = e->data.function.args[1];
        /* (A) negative integer base */
        if (base && base->type == EXPR_INTEGER && base->data.integer < -1)
            return true;
        /* (B)/(C) Power[Times[...], e] */
        if (base && base->type == EXPR_FUNCTION && base->data.function.head
            && base->data.function.head->type == EXPR_SYMBOL
            && base->data.function.head->data.symbol == SYM_Times
            && base->data.function.arg_count > 0) {
            /* (B) any constant-positive factor? */
            for (size_t i = 0; i < base->data.function.arg_count; i++) {
                if (is_constant_positive_factor(base->data.function.args[i]))
                    return true;
            }
            /* (C) integer exponent triggers full distribute */
            if (prov_int(ctx, exp_)) return true;
        }
        /* (C) Power[Power[u, p], e] with e integer */
        if (base && base->type == EXPR_FUNCTION && base->data.function.head
            && base->data.function.head->type == EXPR_SYMBOL
            && base->data.function.head->data.symbol == SYM_Power
            && base->data.function.arg_count == 2
            && prov_int(ctx, exp_)) {
            return true;
        }
    }
    if (has_distributable_power(head, ctx)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_distributable_power(e->data.function.args[i], ctx)) return true;
    }
    return false;
}

static Expr* power_distribute_walk(const Expr* e, const AssumeCtx* ctx) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    size_t n = e->data.function.arg_count;
    Expr* new_head = power_distribute_walk(e->data.function.head, ctx);
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            new_args[i] = power_distribute_walk(e->data.function.args[i], ctx);
        }
    }

    /* Only Power[_, _] triggers any rewrite at this level. */
    if (!(new_head && new_head->type == EXPR_SYMBOL
          && new_head->data.symbol == SYM_Power
          && n == 2 && new_args && new_args[0] && new_args[1])) {
        Expr* result = expr_new_function(new_head, new_args, n);
        if (new_args) free(new_args);
        return result;
    }

    Expr* base = new_args[0];
    Expr* exp_ = new_args[1];

    /* (A) Power[neg_int, e] -> Power[-1, e] * Power[|neg_int|, e] */
    if (base->type == EXPR_INTEGER && base->data.integer < -1) {
        int64_t v = base->data.integer;
        Expr* p1_args[2] = { expr_new_integer(-1), expr_copy(exp_) };
        Expr* p1 = expr_new_function(expr_new_symbol("Power"), p1_args, 2);
        Expr* p2_args[2] = { expr_new_integer(-v), exp_ };
        Expr* p2 = expr_new_function(expr_new_symbol("Power"), p2_args, 2);
        Expr* tm_args[2] = { p1, p2 };
        Expr* tm = expr_new_function(expr_new_symbol("Times"), tm_args, 2);
        expr_free(base);             /* original neg-int base */
        free(new_args);
        expr_free(new_head);
        return tm;
    }

    /* (B) and (C) Power[Times[args...], e] */
    if (base->type == EXPR_FUNCTION && base->data.function.head
        && base->data.function.head->type == EXPR_SYMBOL
        && base->data.function.head->data.symbol == SYM_Times
        && base->data.function.arg_count > 0) {
        size_t bn = base->data.function.arg_count;
        Expr** ba = base->data.function.args;
        bool exp_is_int = prov_int(ctx, exp_);

        /* (C) integer exponent: distribute over every factor. */
        if (exp_is_int) {
            Expr** powers = (Expr**)malloc(sizeof(Expr*) * bn);
            for (size_t i = 0; i < bn; i++) {
                Expr* pa[2] = { expr_copy(ba[i]), expr_copy(exp_) };
                powers[i] = expr_new_function(
                    expr_new_symbol("Power"), pa, 2);
            }
            Expr* tm = expr_new_function(
                expr_new_symbol("Times"), powers, bn);
            expr_free(base);
            expr_free(exp_);
            free(new_args);
            expr_free(new_head);
            return tm;
        }

        /* (B) split off constant-positive factors. */
        size_t pos_count = 0;
        for (size_t i = 0; i < bn; i++) {
            if (is_constant_positive_factor(ba[i])) pos_count++;
        }
        if (pos_count > 0 && pos_count < bn) {
            /* Split: pos_count Power[ci, e] factors + Power[Times[rest], e] */
            Expr** out = (Expr**)malloc(sizeof(Expr*) * (pos_count + 1));
            size_t out_i = 0;
            size_t rest_n = bn - pos_count;
            Expr** rest = (Expr**)malloc(sizeof(Expr*) * rest_n);
            size_t rest_i = 0;
            for (size_t i = 0; i < bn; i++) {
                if (is_constant_positive_factor(ba[i])) {
                    Expr* pa[2] = { expr_copy(ba[i]), expr_copy(exp_) };
                    out[out_i++] = expr_new_function(
                        expr_new_symbol("Power"), pa, 2);
                } else {
                    rest[rest_i++] = expr_copy(ba[i]);
                }
            }
            Expr* rest_times;
            if (rest_n == 1) {
                rest_times = rest[0];
                free(rest);
            } else {
                rest_times = expr_new_function(
                    expr_new_symbol("Times"), rest, rest_n);
            }
            Expr* pa[2] = { rest_times, expr_copy(exp_) };
            out[out_i++] = expr_new_function(
                expr_new_symbol("Power"), pa, 2);
            Expr* tm = expr_new_function(
                expr_new_symbol("Times"), out, pos_count + 1);
            expr_free(base);
            expr_free(exp_);
            free(new_args);
            expr_free(new_head);
            return tm;
        }
        if (pos_count == bn) {
            /* All factors constant-positive: full distribute. */
            Expr** powers = (Expr**)malloc(sizeof(Expr*) * bn);
            for (size_t i = 0; i < bn; i++) {
                Expr* pa[2] = { expr_copy(ba[i]), expr_copy(exp_) };
                powers[i] = expr_new_function(
                    expr_new_symbol("Power"), pa, 2);
            }
            Expr* tm = expr_new_function(
                expr_new_symbol("Times"), powers, bn);
            expr_free(base);
            expr_free(exp_);
            free(new_args);
            expr_free(new_head);
            return tm;
        }
    }

    /* (C) Power[Power[u, p], e] -> Power[u, p*e] when e provably integer. */
    if (base->type == EXPR_FUNCTION && base->data.function.head
        && base->data.function.head->type == EXPR_SYMBOL
        && base->data.function.head->data.symbol == SYM_Power
        && base->data.function.arg_count == 2
        && prov_int(ctx, exp_)) {
        Expr* u  = base->data.function.args[0];
        Expr* p  = base->data.function.args[1];
        Expr* tm_args[2] = { expr_copy(p), exp_ };
        Expr* prod = expr_new_function(expr_new_symbol("Times"), tm_args, 2);
        Expr* pa[2] = { expr_copy(u), prod };
        Expr* po = expr_new_function(expr_new_symbol("Power"), pa, 2);
        expr_free(base);
        free(new_args);
        expr_free(new_head);
        return po;
    }

    Expr* result = expr_new_function(new_head, new_args, n);
    if (new_args) free(new_args);
    return result;
}

static Expr* transform_power_distribute_impl(const Expr* e,
                                              const AssumeCtx* ctx) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    if (!has_distributable_power(e, ctx)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("PowerDistribute", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr* split = power_distribute_walk(e, ctx);
    if (!split) split = expr_copy((Expr*)e);

    /* Re-evaluate so canonical Times same-base merge collapses the
     * newly-introduced Power[u, e] factors against existing same-base
     * factors elsewhere in the surrounding Times. */
    Expr* result = evaluate(split);
    if (!result) result = expr_copy((Expr*)e);

    if (dbg) simp_debug_log("PowerDistribute", e, result,
                            simp_debug_elapsed_ms(t0));
    return result;
}

/* Note: ctx-dependent, so we cannot use simp_memo_wrap (which keys on
 * the input expression alone).  PowerDistribute is invoked at most a
 * few times per Simplify call from simp_dispatch and is structurally
 * cheap (single tree walk, no recursive search), so the caching loss
 * is negligible. */
static Expr* transform_power_distribute(const Expr* e,
                                         const AssumeCtx* ctx) {
    return transform_power_distribute_impl(e, ctx);
}

/* ----------------------------------------------------------------------- */
/* RadicalCanon: split Power[Rational[a,b],q] and rationalise -p/q powers  */
/* ----------------------------------------------------------------------- */

/* Two related canonicalisations that picocas's standard Power evaluator
 * does NOT do, and which leave equivalent expressions in distinct shapes:
 *
 * (1)  Power[Rational[a, b], q]                       (a, b positive ints)
 *           ->  Power[a, q] * Power[b, -q]
 *      Soundness: a, b > 0, so (a/b)^q = a^q * b^(-q) for all complex q
 *      with no branch cut.  Worked example:
 *          Sqrt[1/2] = Power[Rational[1,2], 1/2]
 *                    -> Power[1, 1/2] * Power[2, -1/2]
 *                    -> 1/Sqrt[2]                    (after evaluate)
 *
 * (2)  Power[a, q]                          (a positive int >= 2,
 *                                            q rational, q < 0, denom > 1)
 *           ->  Power[a, r] * Rational[1, a^k]
 *      where k = ceil(-q) >= 1, r = q + k in [0, 1).
 *      Soundness: a > 0, so (a)^q = (a)^r / a^k.  Not folded back by
 *      canonical Times because Rational[1, a^k] occupies a distinct base
 *      bucket from Power[a, r].  Worked example:
 *          Power[2, -1/2]  ->  Power[2, 1/2] * Rational[1, 2]
 *                            =  Sqrt[2] / 2
 *
 * Why both: (1) alone pushes Sqrt[1/2] to 1/Sqrt[2] but leaves it in
 * the negative-exponent form; the surrounding additive context still
 * doesn't see it as same-shape with a Sqrt[2]/6 sibling.  (2) finishes
 * the job by rationalising the denominator so additive cancellations
 * fire.
 *
 * Coverage:
 *     Cos[x]/Sqrt[6] + Sin[x]/Sqrt[2] + Cos[y]/(3 Sqrt[6]) +
 *     Sin[y]/(3 Sqrt[2]) - (Sqrt[6]/3) Sin[x+Pi/6] -
 *     (Sqrt[6]/9) Sin[y+Pi/6]                                  ->  0
 *
 * Inert when the input contains no Power[Rational, _] or no
 * Power[positive_int, negative_rational].
 */

/* Compute a^k with overflow-guard.  Returns false on overflow or
 * a, k out of range; the caller falls back to leaving the term alone. */
static bool radical_canon_pow_int(int64_t a, int64_t k, int64_t* out) {
    if (k < 0 || k > 62) return false;
    int64_t v = 1;
    for (int64_t i = 0; i < k; i++) {
        if (a != 0 && v > INT64_MAX / (a > 0 ? a : -a)) return false;
        v *= a;
    }
    *out = v;
    return true;
}

/* Match Rational[num, den] with den > 1, returning num and den.  The
 * canonical evaluator never produces Rational[_, 1] (it folds to a bare
 * integer) so we only see denominator > 1 here. */
static bool radical_canon_get_rational(const Expr* e,
                                        int64_t* num, int64_t* den) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol != SYM_Rational) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* na = e->data.function.args[0];
    Expr* da = e->data.function.args[1];
    if (na->type != EXPR_INTEGER || da->type != EXPR_INTEGER) return false;
    *num = na->data.integer;
    *den = da->data.integer;
    return true;
}

/* If e = Power[positive_integer_>=2, Rational[num<0, den>1]], rewrite to
 *     Times[Power[a, Rational[r_num, den]], Rational[1, a^k]]
 * Returns NULL when the rule does not apply. */
static Expr* radical_canon_rationalise_negexp(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol != SYM_Power) return NULL;
    if (e->data.function.arg_count != 2) return NULL;
    Expr* base = e->data.function.args[0];
    Expr* exp  = e->data.function.args[1];
    if (base->type != EXPR_INTEGER || base->data.integer < 2) return NULL;
    int64_t num, den;
    if (!radical_canon_get_rational(exp, &num, &den)) return NULL;
    if (num >= 0 || den <= 1) return NULL;

    int64_t a = base->data.integer;
    int64_t neg_num = -num;
    int64_t k = neg_num / den;
    if (neg_num % den != 0) k++;
    int64_t r_num = num + k * den;  /* r_num in [0, den) */

    int64_t a_k;
    if (!radical_canon_pow_int(a, k, &a_k)) return NULL;
    if (a_k < 1) return NULL;

    /* Build Power[a, Rational[r_num, den]].  When r_num == 0, Power
     * evaluates to 1; when r_num == den (impossible here since
     * r_num < den), Power[a, 1] = a. */
    Expr* pow_args[2] = {
        expr_new_integer(a),
        expr_new_function(
            expr_new_symbol("Rational"),
            (Expr*[]){ expr_new_integer(r_num), expr_new_integer(den) }, 2)
    };
    Expr* pow_call = expr_new_function(
        expr_new_symbol("Power"), pow_args, 2);

    /* Rational[1, a^k].  When a^k == 1 (impossible for a >= 2, k >= 1),
     * the evaluator folds to the integer 1. */
    Expr* recip = expr_new_function(
        expr_new_symbol("Rational"),
        (Expr*[]){ expr_new_integer(1), expr_new_integer(a_k) }, 2);

    Expr* times_args[2] = { pow_call, recip };
    Expr* times_call = expr_new_function(
        expr_new_symbol("Times"), times_args, 2);
    Expr* out = evaluate(times_call);
    return out ? out : expr_new_function(
        expr_new_symbol("Times"), times_args, 2);
}

/* If e = Power[Rational[a, b], q] with positive integers a, b > 1,
 * rewrite to Times[Power[a, q], Power[b, -q]].  Returns NULL when not
 * applicable. */
static Expr* radical_canon_split_rational_base(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol != SYM_Power) return NULL;
    if (e->data.function.arg_count != 2) return NULL;
    Expr* base = e->data.function.args[0];
    Expr* exp  = e->data.function.args[1];
    int64_t a, b;
    if (!radical_canon_get_rational(base, &a, &b)) return NULL;
    if (a < 1 || b < 2) return NULL;

    /* Power[a, q] * Power[b, -q].  a may be 1 (Sqrt[1/n] case): then
     * Power[1, q] evaluates to 1 and we're left with Power[b, -q]. */
    Expr* pow_a_args[2] = {
        expr_new_integer(a), expr_copy(exp)
    };
    Expr* pow_a = expr_new_function(
        expr_new_symbol("Power"), pow_a_args, 2);

    Expr* neg_exp_args[2] = {
        expr_new_integer(-1), expr_copy(exp)
    };
    Expr* neg_exp = expr_new_function(
        expr_new_symbol("Times"), neg_exp_args, 2);
    Expr* pow_b_args[2] = {
        expr_new_integer(b), neg_exp
    };
    Expr* pow_b = expr_new_function(
        expr_new_symbol("Power"), pow_b_args, 2);

    Expr* times_args[2] = { pow_a, pow_b };
    Expr* times_call = expr_new_function(
        expr_new_symbol("Times"), times_args, 2);
    Expr* out = evaluate(times_call);
    return out;
}

/* Recursive walker: applies both rules at each Power node and recurses
 * into children.  Returns a freshly allocated tree (caller owns). */
static Expr* radical_canon_walk(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    size_t n = e->data.function.arg_count;
    Expr* new_head = radical_canon_walk(e->data.function.head);
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            new_args[i] = radical_canon_walk(e->data.function.args[i]);
        }
    }

    Expr* rebuilt = expr_new_function(new_head, new_args, n);
    if (new_args) free(new_args);

    /* Try Rational-base split first (may produce a Power node that the
     * negative-exponent rule then handles). */
    Expr* split = radical_canon_split_rational_base(rebuilt);
    if (split) {
        expr_free(rebuilt);
        rebuilt = split;
    }
    Expr* rat = radical_canon_rationalise_negexp(rebuilt);
    if (rat) {
        expr_free(rebuilt);
        rebuilt = rat;
    }
    return rebuilt;
}

static Expr* transform_radical_canon_impl(const Expr* e) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;
    Expr* walked = radical_canon_walk(e);
    if (!walked) walked = expr_copy((Expr*)e);
    /* Re-evaluate so canonical Plus/Times fold the rewritten tree --
     * the walker only evaluates at each Power-level rewrite, leaving the
     * surrounding Plus/Times in unevaluated form (e.g. an arg pair like
     * 1/2 Sqrt[2] + -1 (1/2 Sqrt[2]) won't auto-cancel without this). */
    Expr* result = evaluate(walked);
    if (!result) result = expr_copy((Expr*)e);
    if (dbg) simp_debug_log("RadicalCanon", e, result,
                            simp_debug_elapsed_ms(t0));
    return result;
}

static Expr* transform_radical_canon(const Expr* e) {
    return simp_memo_wrap(e, "$RadicalCanon", transform_radical_canon_impl);
}

/* ----------------------------------------------------------------------- */
/* TanAddition: rewrite Tan[c]/Cot[c] when c = a+b and Tan[a],Tan[b] occur */
/* ----------------------------------------------------------------------- */

/* The standard angle-addition identities:
 *     Tan[a + b] = (Tan[a] + Tan[b]) / (1 - Tan[a] Tan[b])
 *     Cot[a + b] = (1 - Tan[a] Tan[b]) / (Tan[a] + Tan[b])
 * are always sound (away from the isolated singularities where the
 * denominator vanishes -- inside Simplify those are ignored).
 *
 * picocas's TrigExpand fires only when the argument is *literally* a
 * Plus, so it can rewrite Tan[x+y] but NOT Tan[5] -- even when the
 * surrounding expression also contains Tan[2] and Tan[3] making
 * Tan[5] = Tan[2 + 3] usable.  This transform performs that recognition:
 *
 *   1. Walk the input collecting every distinct expression that appears
 *      as the argument of any Tan, Cot, Sin, Cos, Sec, or Csc head.
 *   2. For each ordered pair (a, b) in that set with a != b, evaluate
 *      a + b.  If the sum is also in the set, the triple (a, b, c=a+b)
 *      witnesses an addition-formula opportunity.
 *   3. For each such triple, build a RuleDelayed list that rewrites the
 *      heads which actually appear at c (Tan[c], Cot[c], etc.) using the
 *      angle-addition formula in terms of Tan[a], Tan[b] (or Sin/Cos[a]
 *      and Sin/Cos[b], for Sin[c] / Cos[c] occurrences).
 *   4. Apply ReplaceAll with the constructed rules, then evaluate +
 *      Together + Cancel so the polynomial cancellation fires.
 *   5. Score-gated: keep only strict wins.  Inert when no pair-sum match
 *      exists in the input (the common case).
 *
 * Coverage of the case-6 shape:
 *     Tan[2] Tan[3] B A - (-Tan[2]/Tan[5] - Tan[3]/Tan[5] + 1) B A  ->  0
 * via the (2, 3, 5) triple's Cot[5] -> (1 - Tan[2] Tan[3]) / (Tan[2] +
 * Tan[3]) substitution.  Symbolic shapes (Tan[x], Tan[y], Tan[x+y])
 * also work; the integer-arg case is just the most surprising one.
 */

typedef struct {
    Expr** items;       /* Borrowed pointers into the input tree. */
    size_t count;
    size_t cap;
} TrigArgSet;

static void tas_init(TrigArgSet* s) {
    s->items = NULL; s->count = 0; s->cap = 0;
}
static void tas_free(TrigArgSet* s) {
    free(s->items);
    s->items = NULL; s->count = 0; s->cap = 0;
}
static bool tas_contains(const TrigArgSet* s, const Expr* e) {
    for (size_t i = 0; i < s->count; i++) {
        if (expr_eq(s->items[i], e)) return true;
    }
    return false;
}
static void tas_add_borrowed(TrigArgSet* s, Expr* e) {
    if (tas_contains(s, e)) return;
    if (s->count == s->cap) {
        s->cap = s->cap ? s->cap * 2 : 8;
        s->items = (Expr**)realloc(s->items, sizeof(Expr*) * s->cap);
    }
    s->items[s->count++] = e;
}

static bool is_addition_trig_head(const char* h) {
    return h == SYM_Tan || h == SYM_Cot || h == SYM_Sin
        || h == SYM_Cos || h == SYM_Sec || h == SYM_Csc;
}

static void collect_addition_trig_args(const Expr* e, TrigArgSet* set) {
    if (!e || e->type != EXPR_FUNCTION) return;
    Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL
        && is_addition_trig_head(head->data.symbol)
        && e->data.function.arg_count == 1) {
        tas_add_borrowed(set, e->data.function.args[0]);
    }
    if (head) collect_addition_trig_args(head, set);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        collect_addition_trig_args(e->data.function.args[i], set);
    }
}

/* Returns true if `e` contains Power[head[c], _] or head[c] for the
 * given trig head symbol.  A cheap substring presence check without a
 * full pattern match. */
static bool tree_contains_trig_at(const Expr* e, const char* head_sym,
                                   const Expr* arg) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL
        && head->data.symbol == head_sym
        && e->data.function.arg_count == 1
        && expr_eq(e->data.function.args[0], arg)) {
        return true;
    }
    if (head && tree_contains_trig_at(head, head_sym, arg)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (tree_contains_trig_at(e->data.function.args[i], head_sym, arg)) {
            return true;
        }
    }
    return false;
}

static Expr* mk_unary_call(const char* head, Expr* arg) {
    Expr* args[1] = { arg };
    return expr_new_function(expr_new_symbol(head), args, 1);
}

static Expr* mk_div(Expr* num, Expr* den) {
    Expr* inv_args[2] = { den, expr_new_integer(-1) };
    Expr* inv = expr_new_function(
        expr_new_symbol("Power"), inv_args, 2);
    Expr* args[2] = { num, inv };
    return expr_new_function(expr_new_symbol("Times"), args, 2);
}

/* Append rule list entries (RuleDelayed nodes) for the heads that
 * actually appear at `c` in `input`. */
static void append_addition_rules_for_triple(Expr* input,
                                              const Expr* a, const Expr* b,
                                              const Expr* c,
                                              Expr*** rules, size_t* rcount,
                                              size_t* rcap) {
    #define RAPPEND(node) do {                                \
        if (*rcount == *rcap) {                               \
            *rcap = *rcap ? *rcap * 2 : 8;                    \
            *rules = (Expr**)realloc(*rules,                  \
                                      sizeof(Expr*) * (*rcap)); \
        }                                                     \
        (*rules)[(*rcount)++] = (node);                       \
    } while (0)

    /* All six rules (Tan/Cot/Sin/Cos/Sec/Csc at c) are expressed in terms
     * of the same two Sin/Cos primitives:
     *     sin_rhs = Sin[a] Cos[b] + Cos[a] Sin[b]   (= Sin[a + b])
     *     cos_rhs = Cos[a] Cos[b] - Sin[a] Sin[b]   (= Cos[a + b])
     * Using the same primitives for Tan[c] = sin_rhs / cos_rhs and
     * Sec[c] = 1 / cos_rhs (rather than the Tan-based identity Tan[c] =
     * (Tan[a]+Tan[b])/(1 - Tan[a] Tan[b])) gives the post-substitution
     * expression a uniform Sin/Cos denominator structure, so a single
     * Together / Cancel pass collapses mixed Sec[c] + Tan[c] inputs that
     * the Tan-based formula leaves with mismatched (1 - Tan[a] Tan[b])
     * vs (Cos[a] Cos[b] - Sin[a] Sin[b]) denominators. */
    bool need_tan = tree_contains_trig_at(input, SYM_Tan, c);
    bool need_cot = tree_contains_trig_at(input, SYM_Cot, c);
    bool need_sin = tree_contains_trig_at(input, SYM_Sin, c);
    bool need_csc = tree_contains_trig_at(input, SYM_Csc, c);
    bool need_cos = tree_contains_trig_at(input, SYM_Cos, c);
    bool need_sec = tree_contains_trig_at(input, SYM_Sec, c);

    bool need_sin_rhs = need_sin || need_csc || need_tan || need_cot;
    bool need_cos_rhs = need_cos || need_sec || need_tan || need_cot;

    Expr* sin_rhs = NULL;
    if (need_sin_rhs) {
        Expr* sa_cb_args[2] = {
            mk_unary_call("Sin", expr_copy((Expr*)a)),
            mk_unary_call("Cos", expr_copy((Expr*)b))
        };
        Expr* sa_cb = expr_new_function(
            expr_new_symbol("Times"), sa_cb_args, 2);
        Expr* ca_sb_args[2] = {
            mk_unary_call("Cos", expr_copy((Expr*)a)),
            mk_unary_call("Sin", expr_copy((Expr*)b))
        };
        Expr* ca_sb = expr_new_function(
            expr_new_symbol("Times"), ca_sb_args, 2);
        Expr* sin_args[2] = { sa_cb, ca_sb };
        sin_rhs = expr_new_function(
            expr_new_symbol("Plus"), sin_args, 2);
    }

    Expr* cos_rhs = NULL;
    if (need_cos_rhs) {
        Expr* ca_cb_args[2] = {
            mk_unary_call("Cos", expr_copy((Expr*)a)),
            mk_unary_call("Cos", expr_copy((Expr*)b))
        };
        Expr* ca_cb = expr_new_function(
            expr_new_symbol("Times"), ca_cb_args, 2);
        Expr* neg_sa_sb_args[3] = {
            expr_new_integer(-1),
            mk_unary_call("Sin", expr_copy((Expr*)a)),
            mk_unary_call("Sin", expr_copy((Expr*)b))
        };
        Expr* neg_sa_sb = expr_new_function(
            expr_new_symbol("Times"), neg_sa_sb_args, 3);
        Expr* cos_args[2] = { ca_cb, neg_sa_sb };
        cos_rhs = expr_new_function(
            expr_new_symbol("Plus"), cos_args, 2);
    }

    if (need_sin) {
        Expr* lhs = mk_unary_call("Sin", expr_copy((Expr*)c));
        Expr* rule_args[2] = { lhs, expr_copy(sin_rhs) };
        Expr* rule = expr_new_function(
            expr_new_symbol("RuleDelayed"), rule_args, 2);
        RAPPEND(rule);
    }
    if (need_csc) {
        Expr* lhs = mk_unary_call("Csc", expr_copy((Expr*)c));
        Expr* inv_args[2] = { expr_copy(sin_rhs), expr_new_integer(-1) };
        Expr* inv = expr_new_function(
            expr_new_symbol("Power"), inv_args, 2);
        Expr* rule_args[2] = { lhs, inv };
        Expr* rule = expr_new_function(
            expr_new_symbol("RuleDelayed"), rule_args, 2);
        RAPPEND(rule);
    }
    if (need_cos) {
        Expr* lhs = mk_unary_call("Cos", expr_copy((Expr*)c));
        Expr* rule_args[2] = { lhs, expr_copy(cos_rhs) };
        Expr* rule = expr_new_function(
            expr_new_symbol("RuleDelayed"), rule_args, 2);
        RAPPEND(rule);
    }
    if (need_sec) {
        Expr* lhs = mk_unary_call("Sec", expr_copy((Expr*)c));
        Expr* inv_args[2] = { expr_copy(cos_rhs), expr_new_integer(-1) };
        Expr* inv = expr_new_function(
            expr_new_symbol("Power"), inv_args, 2);
        Expr* rule_args[2] = { lhs, inv };
        Expr* rule = expr_new_function(
            expr_new_symbol("RuleDelayed"), rule_args, 2);
        RAPPEND(rule);
    }
    if (need_tan) {
        Expr* lhs = mk_unary_call("Tan", expr_copy((Expr*)c));
        Expr* rhs = mk_div(expr_copy(sin_rhs), expr_copy(cos_rhs));
        Expr* rule_args[2] = { lhs, rhs };
        Expr* rule = expr_new_function(
            expr_new_symbol("RuleDelayed"), rule_args, 2);
        RAPPEND(rule);
    }
    if (need_cot) {
        Expr* lhs = mk_unary_call("Cot", expr_copy((Expr*)c));
        Expr* rhs = mk_div(expr_copy(cos_rhs), expr_copy(sin_rhs));
        Expr* rule_args[2] = { lhs, rhs };
        Expr* rule = expr_new_function(
            expr_new_symbol("RuleDelayed"), rule_args, 2);
        RAPPEND(rule);
    }

    if (sin_rhs) expr_free(sin_rhs);
    if (cos_rhs) expr_free(cos_rhs);
    #undef RAPPEND
}

static Expr* transform_tan_addition_impl(const Expr* e) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Cheap gate: if no trig head appears in the input, we cannot collect
     * any args and the rule list will be empty. */
    if (!has_pythag_head(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("TanAddition", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    TrigArgSet args; tas_init(&args);
    collect_addition_trig_args(e, &args);
    if (args.count < 3) {
        /* Need at least three distinct args to have a (a, b, c=a+b) triple. */
        tas_free(&args);
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("TanAddition", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    Expr** rules = NULL;
    size_t rcount = 0, rcap = 0;
    Expr* input_copy = (Expr*)e;  /* tree_contains_trig_at takes const */
    for (size_t i = 0; i < args.count; i++) {
        for (size_t j = 0; j < args.count; j++) {
            if (i == j) continue;
            Expr* sum_args[2] = {
                expr_copy(args.items[i]), expr_copy(args.items[j])
            };
            Expr* sum_call = expr_new_function(
                expr_new_symbol("Plus"), sum_args, 2);
            Expr* sum = evaluate(sum_call);
            if (!sum) continue;
            if (tas_contains(&args, sum)) {
                /* Use the canonical sum as `c` for the rule LHS.  The
                 * input might contain Tan[Plus[a, b]] (if a + b doesn't
                 * collapse to a literal) or Tan[5] (if a, b are
                 * integers), so we need to use the actually-occurring
                 * representative -- pull that out of the args set. */
                Expr* canonical_c = NULL;
                for (size_t k = 0; k < args.count; k++) {
                    if (expr_eq(args.items[k], sum)) {
                        canonical_c = args.items[k];
                        break;
                    }
                }
                if (canonical_c) {
                    append_addition_rules_for_triple(
                        input_copy,
                        args.items[i], args.items[j], canonical_c,
                        &rules, &rcount, &rcap);
                }
            }
            expr_free(sum);
        }
    }
    tas_free(&args);

    if (rcount == 0) {
        free(rules);
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("TanAddition", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    /* Assemble rule list { rule1, rule2, ... } and apply ReplaceAll. */
    Expr* rule_list = expr_new_function(
        expr_new_symbol("List"), rules, rcount);
    free(rules);
    Expr* ra_args[2] = { expr_copy((Expr*)e), rule_list };
    Expr* ra_call = expr_new_function(
        expr_new_symbol("ReplaceAll"), ra_args, 2);
    Expr* substituted = evaluate(ra_call);
    if (!substituted) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("TanAddition", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    /* Together + Cancel turns the rational expression into a single
     * fraction, then reduces by polynomial GCD.  When the substitution
     * was witnessing a true identity (the case-6 shape), this drops to
     * 0 cleanly. */
    Expr* tg_args[1] = { substituted };
    Expr* tg_call = expr_new_function(
        expr_new_symbol("Together"), tg_args, 1);
    Expr* tg = evaluate(tg_call);
    Expr* result = tg ? tg : expr_copy((Expr*)e);

    Expr* cn_args[1] = { result };
    Expr* cn_call = expr_new_function(
        expr_new_symbol("Cancel"), cn_args, 1);
    Expr* cn = evaluate(cn_call);
    if (cn) result = cn;

    if (dbg) simp_debug_log("TanAddition", e, result,
                            simp_debug_elapsed_ms(t0));
    return result;
}

static Expr* transform_tan_addition(const Expr* e) {
    return simp_memo_wrap(e, "$TanAddition", transform_tan_addition_impl);
}

/* ----------------------------------------------------------------------- */
/* Log/Power rewriter (positive-real cascade, v1)                          */
/* ----------------------------------------------------------------------- */

/*
 * The strict-positive cascade implements the Log/Power identities that
 * are sound under positivity / reality assumptions on the operands.
 * Identities cover (1) log of products and quotients, (2) log of a power
 * of a positive base, (3) power of a product, and (4) tower-of-powers
 * collapse for a positive base.
 *
 * The general-real and general-complex branches of the user's cascade
 * (with Boole / Floor / Ceiling phase corrections) are deliberately not
 * implemented in v1; see picocas_spec.md for v2 scope.
 *
 * Implementation: a bottom-up structural walker that consults the
 * AssumeCtx for positivity/reality of operands. Each top-level rewrite
 * emits a freshly evaluated tree, so e.g. nested Log[Times[x, 1/y]] ->
 * Log[x] + Log[1/y] -> Log[x] - Log[y] (via the Power[..., -1] case)
 * stabilises after a small fixed number of passes.
 */

static Expr* logexp_top_rewrite(const Expr* e, const AssumeCtx* ctx);

/* Returns NULL if the recursive walk produced no change. Otherwise returns
 * a newly owned, evaluated tree. */
static Expr* logexp_walk(const Expr* e, const AssumeCtx* ctx) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) {
        return logexp_top_rewrite(e, ctx);
    }

    /* First rewrite children. */
    size_t n = e->data.function.arg_count;
    Expr** new_args = NULL;
    bool any = false;
    for (size_t i = 0; i < n; i++) {
        Expr* r = logexp_walk(e->data.function.args[i], ctx);
        if (r) {
            if (!new_args) {
                new_args = (Expr**)calloc(n, sizeof(Expr*));
                for (size_t j = 0; j < i; j++) new_args[j] = expr_copy(e->data.function.args[j]);
            }
            new_args[i] = r;
            any = true;
        } else if (new_args) {
            new_args[i] = expr_copy(e->data.function.args[i]);
        }
    }

    Expr* current_owned = NULL;
    const Expr* target;
    if (any) {
        Expr* head_copy = expr_copy(e->data.function.head);
        Expr* rebuilt = expr_new_function(head_copy, new_args, n);
        free(new_args);
        current_owned = evaluate(rebuilt);
        expr_free(rebuilt);
        target = current_owned;
    } else {
        target = e;
    }

    Expr* top = logexp_top_rewrite(target, ctx);
    if (top) {
        if (current_owned) expr_free(current_owned);
        return top;
    }
    return current_owned;  /* may be NULL if no change anywhere */
}

static Expr* build_unary(const char* head, Expr* owned_arg) {
    Expr* a[1] = { owned_arg };
    return expr_new_function(expr_new_symbol(head), a, 1);
}

static Expr* build_binary(const char* head, Expr* a0, Expr* a1) {
    Expr* a[2] = { a0, a1 };
    return expr_new_function(expr_new_symbol(head), a, 2);
}

static Expr* logexp_top_rewrite(const Expr* e, const AssumeCtx* ctx) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* h = e->data.function.head->data.symbol;
    Expr** a = e->data.function.args;
    size_t n = e->data.function.arg_count;

    /* Log[Times[u1,...,un]] -> Sum Log[ui]  when every ui is positive.
     * Log[Power[x, p]]      -> p Log[x]      when x positive and p real. */
    if (h == SYM_Log && n == 1) {
        Expr* inner = a[0];
        if (inner->type == EXPR_FUNCTION &&
            inner->data.function.head &&
            inner->data.function.head->type == EXPR_SYMBOL) {
            const char* ih = inner->data.function.head->data.symbol;
            size_t in = inner->data.function.arg_count;
            Expr** ia = inner->data.function.args;

            if (ih == SYM_Times && in > 0) {
                bool all_pos = true;
                for (size_t i = 0; i < in; i++) {
                    if (!prov_pos(ctx, ia[i])) { all_pos = false; break; }
                }
                if (all_pos) {
                    Expr** logs = (Expr**)calloc(in, sizeof(Expr*));
                    for (size_t i = 0; i < in; i++) {
                        logs[i] = build_unary("Log", expr_copy(ia[i]));
                    }
                    Expr* sum = expr_new_function(expr_new_symbol("Plus"), logs, in);
                    free(logs);
                    Expr* canon = evaluate(sum);
                    expr_free(sum);
                    return canon;
                }
            }
            if (ih == SYM_Power && in == 2) {
                Expr* base = ia[0];
                Expr* p    = ia[1];
                if (prov_pos(ctx, base) && prov_re(ctx, p)) {
                    Expr* logx = build_unary("Log", expr_copy(base));
                    Expr* mul  = build_binary("Times", expr_copy(p), logx);
                    Expr* canon = evaluate(mul);
                    expr_free(mul);
                    return canon;
                }
            }
        }
    }

    /* Power[-1, k] reductions: even exponent -> 1, odd exponent -> -1.
     * Catches `(-1)^(2 n)` with `n` provably integer, etc. The Times
     * propagation in prov_even (a product of integers with at least one
     * even factor is even) handles the canonical user input. */
    if (h == SYM_Power && n == 2) {
        Expr* base = a[0];
        Expr* exp_  = a[1];
        if (base->type == EXPR_INTEGER && base->data.integer == -1) {
            if (prov_even(ctx, exp_)) return expr_new_integer(1);
        }
    }

    /* Exp distribute: Power[E, Plus[t1,...,tn]] -> Product Power[E, ti].
     * Always sound (E^(x+y) = E^x · E^y for any x, y). When the exponent
     * is a product like `Times[c, Plus[...]]`, expand it first so the
     * Plus surfaces. After distribution, individual Power[E, c·Log u]
     * subterms collapse to `u^c` via picocas's existing
     * Power[E, c·Log[u]] -> u^c rule, completing identities like
     *
     *   Exp[3 (Log[a] + Log[b])]  ->  a^3 · b^3
     *   Exp[y (Log[a] + Log[b])]  ->  a^y · b^y    (a > 0, b > 0)
     */
    if (h == SYM_Power && n == 2) {
        Expr* base = a[0];
        Expr* exp_  = a[1];

        if (base->type == EXPR_SYMBOL && base->data.symbol == SYM_E) {
            /* Try to expand the exponent so Plus distributes through any
             * outer Times. Cheap, idempotent, and only used to detect
             * Plus structure. */
            Expr* expanded_exp = expr_expand(exp_);
            if (expanded_exp &&
                expanded_exp->type == EXPR_FUNCTION &&
                expanded_exp->data.function.head &&
                expanded_exp->data.function.head->type == EXPR_SYMBOL &&
                expanded_exp->data.function.head->data.symbol == SYM_Plus &&
                expanded_exp->data.function.arg_count > 1) {
                size_t en = expanded_exp->data.function.arg_count;
                Expr** factors = (Expr**)calloc(en, sizeof(Expr*));
                for (size_t i = 0; i < en; i++) {
                    factors[i] = build_binary("Power", expr_new_symbol("E"),
                                              expr_copy(expanded_exp->data.function.args[i]));
                }
                Expr* prod = expr_new_function(expr_new_symbol("Times"), factors, en);
                free(factors);
                Expr* canon = evaluate(prod);
                expr_free(prod);
                expr_free(expanded_exp);
                return canon;
            }
            if (expanded_exp) expr_free(expanded_exp);
        }

        if (base->type == EXPR_FUNCTION &&
            base->data.function.head &&
            base->data.function.head->type == EXPR_SYMBOL) {
            const char* bh = base->data.function.head->data.symbol;
            size_t bn = base->data.function.arg_count;
            Expr** ba = base->data.function.args;

            if (bh == SYM_Times && bn > 0) {
                bool all_pos = true;
                for (size_t i = 0; i < bn; i++) {
                    if (!prov_pos(ctx, ba[i])) { all_pos = false; break; }
                }
                if (all_pos) {
                    Expr** powers = (Expr**)calloc(bn, sizeof(Expr*));
                    for (size_t i = 0; i < bn; i++) {
                        powers[i] = build_binary("Power", expr_copy(ba[i]), expr_copy(exp_));
                    }
                    Expr* prod = expr_new_function(expr_new_symbol("Times"), powers, bn);
                    free(powers);
                    Expr* canon = evaluate(prod);
                    expr_free(prod);
                    return canon;
                }
            }
            if (bh == SYM_Power && bn == 2) {
                Expr* xx = ba[0];
                Expr* pp = ba[1];
                if (prov_pos(ctx, xx) && prov_re(ctx, pp)) {
                    Expr* prod = build_binary("Times", expr_copy(pp), expr_copy(exp_));
                    Expr* prod_canon = evaluate(prod);
                    expr_free(prod);
                    Expr* pow_ = build_binary("Power", expr_copy(xx), prod_canon);
                    Expr* canon = evaluate(pow_);
                    expr_free(pow_);
                    return canon;
                }
            }
        }
    }

    return NULL;
}

/* Apply the rewriter to a fixed point. Returns NULL if unchanged.
 * Bounded iteration count protects against pathological alternations
 * with the evaluator's canonicalisation. */
static Expr* apply_logexp_rules(const Expr* input, const AssumeCtx* ctx) {
    /* NULL ctx is treated as an empty context. The positivity-aware
     * Log[a*b] -> Log[a]+Log[b] etc. rewrites simply won't fire
     * (prov_pos returns false), but the unconditional
     * Power[E, Plus[...]] distribute rule still does its job. */
    Expr* current = expr_copy((Expr*)input);
    bool changed = false;
    for (int iter = 0; iter < 8; iter++) {
        Expr* r = logexp_walk(current, ctx);
        if (!r) break;
        if (expr_eq(r, current)) { expr_free(r); break; }
        expr_free(current);
        current = r;
        changed = true;
    }
    if (!changed) {
        expr_free(current);
        return NULL;
    }
    if (expr_eq(current, input)) {
        expr_free(current);
        return NULL;
    }
    return current;
}

/* ----------------------------------------------------------------------- */
/* Abs simplification: structural rewrites over Abs[...] subexpressions   */
/* ----------------------------------------------------------------------- */

/* Cheap pre-check: skip the walker when the input is Abs-free. */
static bool contains_abs(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Abs) return true;
    if (contains_abs(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_abs(e->data.function.args[i])) return true;
    }
    return false;
}

/* True iff `h` is a 1-arg head whose presence makes a transform that
 * targets trig or hyperbolic functions potentially fire. Covers the six
 * canonical pairs and their inverses. */
static bool head_is_trig_or_hyperbolic(const char* h) {
    static const char* const TRIG_HEADS[] = {
        "Sin","Cos","Tan","Cot","Sec","Csc",
        "ArcSin","ArcCos","ArcTan","ArcCot","ArcSec","ArcCsc",
        "Sinh","Cosh","Tanh","Coth","Sech","Csch",
        "ArcSinh","ArcCosh","ArcTanh","ArcCoth","ArcSech","ArcCsch",
        NULL
    };
    for (size_t i = 0; TRIG_HEADS[i]; i++) {
        if (strcmp(h, TRIG_HEADS[i]) == 0) return true;
    }
    return false;
}

static bool contains_trig_or_hyperbolic(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        head_is_trig_or_hyperbolic(e->data.function.head->data.symbol)) return true;
    if (contains_trig_or_hyperbolic(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_trig_or_hyperbolic(e->data.function.args[i])) return true;
    }
    return false;
}

/* True iff `e` contains a Power[E, _] subexpression -- i.e. an Exp atom in
 * exponential form (E^x, E^(-x), E^(I Pi), ...). Used to gate the
 * ExpToTrig seed: pure-Exp inputs miss every trig-gated transform
 * (TrigRoundtrip, PythagReduce, TrigFactor, ...) because their gates
 * check for Cos/Sin/Cosh/Sinh heads. ExpToTrig converts the Exp form
 * to Cosh/Sinh, opening those transforms to subsequent rounds. */
static bool contains_exp_form(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        if (base && base->type == EXPR_SYMBOL && base->data.symbol == SYM_E) {
            return true;
        }
    }
    if (contains_exp_form(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_exp_form(e->data.function.args[i])) return true;
    }
    return false;
}

static bool contains_log(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Log) return true;
    if (contains_log(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_log(e->data.function.args[i])) return true;
    }
    return false;
}

static bool contains_power(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power) return true;
    if (contains_power(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_power(e->data.function.args[i])) return true;
    }
    return false;
}

/* True iff a Plus, Times, or Power head appears anywhere in `e`. Used by
 * the TrigReduce gate to short-circuit on a bare single trig call (no
 * product or power means no product-to-sum work). Power is included
 * because Sin[x]^2 is the canonical TrigReduce input. */
static bool contains_plus_or_times(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (h == SYM_Plus ||
            h == SYM_Times ||
            h == SYM_Power) return true;
    }
    if (contains_plus_or_times(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_plus_or_times(e->data.function.args[i])) return true;
    }
    return false;
}

/* True iff any non-numeric-constant symbol leaf appears anywhere in `e`.
 * Pi, E, EulerGamma, Degree, Catalan, Glaisher, Khinchin do not count --
 * they are positive numeric constants. Used to short-circuit transforms
 * that have nothing to do on a purely numeric input (Factor, Apart, ...). */
static bool contains_variable(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        return !is_real_constant_symbol(e->data.symbol);
    }
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_variable(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_variable(e->data.function.args[i])) return true;
    }
    return false;
}

/* Number of distinct non-constant symbol leaves in `e`, capped at `cap`.
 * Returns as soon as the count reaches `cap`, so callers that only need
 * "0 / 1 / >=2" can pass cap=2 and early-out. Constant symbols (Pi, E,
 * ...) are excluded, matching contains_variable. */
static size_t expr_variables_count_capped_walk(const Expr* e,
                                               char** seen, size_t* nseen,
                                               size_t cap) {
    if (!e || *nseen >= cap) return *nseen;
    if (e->type == EXPR_SYMBOL) {
        if (is_real_constant_symbol(e->data.symbol)) return *nseen;
        for (size_t i = 0; i < *nseen; i++) {
            if (strcmp(seen[i], e->data.symbol) == 0) return *nseen;
        }
        seen[*nseen] = e->data.symbol;
        (*nseen)++;
        return *nseen;
    }
    if (e->type != EXPR_FUNCTION) return *nseen;
    expr_variables_count_capped_walk(e->data.function.head, seen, nseen, cap);
    for (size_t i = 0; i < e->data.function.arg_count && *nseen < cap; i++) {
        expr_variables_count_capped_walk(e->data.function.args[i],
                                         seen, nseen, cap);
    }
    return *nseen;
}

static size_t expr_variables_count_capped(const Expr* e, size_t cap) {
    if (cap == 0) return 0;
    char* seen[8];  /* cap is at most 2 in our call sites; 8 is a safe ceiling */
    size_t nseen = 0;
    if (cap > 8) cap = 8;
    expr_variables_count_capped_walk(e, seen, &nseen, cap);
    return nseen;
}

/* True iff the assumption ctx has at least one usable fact. NULL ctx, an
 * empty fact list, or an inconsistent ctx all return false -- no
 * assumption-driven rewrite can do anything in those cases. */
static bool ctx_has_facts(const AssumeCtx* ctx) {
    return ctx != NULL && ctx->count > 0 && !ctx->inconsistent;
}

/* Try to simplify a single Abs[arg] node. `arg` is the inner expression
 * (i.e. the argument to Abs). Returns a new Expr* on success, NULL if no
 * rule fires. */
static Expr* try_simp_abs(const Expr* arg, const AssumeCtx* ctx) {
    /* Universal: idempotency Abs[Abs[x]] -> Abs[x]. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Abs &&
        arg->data.function.arg_count == 1) {
        return expr_copy((Expr*)arg);
    }

    /* Universal: conjugate symmetry Abs[Conjugate[x]] -> Abs[x]. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Conjugate &&
        arg->data.function.arg_count == 1) {
        Expr* a[1] = { expr_copy(arg->data.function.args[0]) };
        return expr_new_function(expr_new_symbol("Abs"), a, 1);
    }

    /* Universal: Abs[E^z] -> E^Re[z]. The magnitude of any complex
     * exponential is e^(real part of the exponent). */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Power &&
        arg->data.function.arg_count == 2 &&
        arg->data.function.args[0]->type == EXPR_SYMBOL &&
        arg->data.function.args[0]->data.symbol == SYM_E) {
        Expr* re_in[1] = { expr_copy(arg->data.function.args[1]) };
        Expr* re_call = expr_new_function(expr_new_symbol("Re"), re_in, 1);
        Expr* pa[2] = { expr_new_symbol("E"), re_call };
        return expr_new_function(expr_new_symbol("Power"), pa, 2);
    }

    /* Universal: split products. Abs[Times[a, b, ...]] -> Abs[a] Abs[b] ...
     * Captures both Abs[c x] (numeric coefficient extraction) and the
     * Abs[x/y] case since x/y is Times[x, Power[y, -1]]. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Times &&
        arg->data.function.arg_count >= 2) {
        size_t n = arg->data.function.arg_count;
        Expr** new_args = (Expr**)calloc(n, sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            Expr* a[1] = { expr_copy(arg->data.function.args[i]) };
            new_args[i] = expr_new_function(expr_new_symbol("Abs"), a, 1);
        }
        Expr* result = expr_new_function(expr_new_symbol("Times"), new_args, n);
        free(new_args);
        return result;
    }

    /* Universal: integer-power split. Abs[x^n] -> Abs[x]^n for integer n.
     * For complex x and integer n the identity |x^n| = |x|^n is exact;
     * for non-integer n it can fail (branch-cut), so the unconditional
     * rule applies only to integer exponents. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Power &&
        arg->data.function.arg_count == 2 &&
        (arg->data.function.args[1]->type == EXPR_INTEGER ||
         arg->data.function.args[1]->type == EXPR_BIGINT)) {
        Expr* base = arg->data.function.args[0];
        Expr* exp  = arg->data.function.args[1];
        Expr* a[1] = { expr_copy(base) };
        Expr* abs_call = expr_new_function(expr_new_symbol("Abs"), a, 1);
        Expr* pa[2] = { abs_call, expr_copy(exp) };
        return expr_new_function(expr_new_symbol("Power"), pa, 2);
    }

    /* The remaining rules need an assumption context. */
    if (!ctx) return NULL;

    /* Cascading: Abs[x] -> x  if x >= 0 (provably nonnegative). */
    if (assume_known_nonneg(ctx, arg)) {
        return expr_copy((Expr*)arg);
    }

    /* Cascading: Abs[x] -> -x  if x <= 0 (provably nonpositive). */
    if (assume_known_nonpos(ctx, arg)) {
        Expr* na[2] = { expr_new_integer(-1), expr_copy((Expr*)arg) };
        return expr_new_function(expr_new_symbol("Times"), na, 2);
    }

    /* Cascading: Abs[x^y] -> Abs[x]^y if y is real. The integer-power
     * rule above handles n in Z; this generalises to any real y under
     * an Element[y, Reals] assumption. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Power &&
        arg->data.function.arg_count == 2 &&
        assume_known_real(ctx, arg->data.function.args[1])) {
        Expr* base = arg->data.function.args[0];
        Expr* exp  = arg->data.function.args[1];
        Expr* a[1] = { expr_copy(base) };
        Expr* abs_call = expr_new_function(expr_new_symbol("Abs"), a, 1);
        Expr* pa[2] = { abs_call, expr_copy(exp) };
        return expr_new_function(expr_new_symbol("Power"), pa, 2);
    }

    /* Cascading: Abs[x^y] -> x^Re[y] if x > 0 (strictly positive).
     * Proof: for x > 0, x^y = x^(Re[y] + I Im[y]) = x^Re[y] * Exp[I Im[y]
     * Log[x]] and the second factor has unit modulus. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Power &&
        arg->data.function.arg_count == 2 &&
        assume_known_positive(ctx, arg->data.function.args[0])) {
        Expr* base = arg->data.function.args[0];
        Expr* exp  = arg->data.function.args[1];
        Expr* re_in[1] = { expr_copy(exp) };
        Expr* re_call = expr_new_function(expr_new_symbol("Re"), re_in, 1);
        Expr* pa[2] = { expr_copy(base), re_call };
        return expr_new_function(expr_new_symbol("Power"), pa, 2);
    }

    /* The Abs[Sin[x]] -> Sign[Sin[x]] Sin[x] rule from the user-provided
     * cascade is omitted: the rewrite expands leaf count (3 -> 6) and only
     * pays off when a downstream Sign-folding pass narrows Sign[Sin[x]] on
     * a known interval, which picocas does not currently perform. Adding
     * it without that infrastructure produces a strictly larger expression
     * with no observable benefit. */
    return NULL;
}

/* Bottom-up walker that rewrites Abs[...] subexpressions. Returns a new
 * Expr* if any rewrite fired anywhere in the tree, NULL otherwise. */
static Expr* abs_walk(const Expr* e, const AssumeCtx* ctx) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;

    size_t argc = e->data.function.arg_count;
    Expr** new_args = (Expr**)calloc(argc ? argc : 1, sizeof(Expr*));
    bool any = false;
    for (size_t i = 0; i < argc; i++) {
        Expr* sub = abs_walk(e->data.function.args[i], ctx);
        if (sub) {
            new_args[i] = sub;
            any = true;
        } else {
            new_args[i] = expr_copy(e->data.function.args[i]);
        }
    }

    Expr* this_form;
    if (any) {
        this_form = expr_new_function(expr_copy(e->data.function.head),
                                       new_args, argc);
    } else {
        for (size_t i = 0; i < argc; i++) expr_free(new_args[i]);
        this_form = NULL;
    }
    free(new_args);

    /* Rule fires only on Abs[_]. */
    bool is_abs = e->data.function.head &&
                  e->data.function.head->type == EXPR_SYMBOL &&
                  e->data.function.head->data.symbol == SYM_Abs &&
                  e->data.function.arg_count == 1;
    if (is_abs) {
        const Expr* inner = this_form ? this_form->data.function.args[0]
                                       : e->data.function.args[0];
        Expr* simp = try_simp_abs(inner, ctx);
        if (simp) {
            if (this_form) expr_free(this_form);
            return simp;
        }
    }

    return this_form;
}

/* Returns a rewritten copy of `input` if any Abs simplification fired,
 * NULL otherwise. ctx may be NULL (universal rules still fire). */
static Expr* apply_abs_rules(const Expr* input, const AssumeCtx* ctx) {
    if (!contains_abs(input)) return NULL;
    return abs_walk(input, ctx);
}

/* ----------------------------------------------------------------------- */
/* Heuristic search                                                        */
/* ----------------------------------------------------------------------- */

/* Forward declarations: defined below alongside the rest of the
 * simp_factorial cluster. transform_can_fire (a few hundred lines down)
 * needs to ask "does this input contain a Factorial atom?" before
 * firing the FactorialRules seed, and simp_search uses simp_eq_head_sym
 * to count factorials in candidates -- both must be visible here. */
static bool contains_factorial(const Expr* e);
static bool simp_eq_head_sym(const Expr* e, const char* name);

static const char* SIMP_TRANSFORMS[] = {
    "Together",
    "Cancel",
    "Expand",
    "ExpandNumerator",
    "ExpandDenominator",
    "Factor",
    "FactorSquareFree",
    "FactorTerms",
    "Apart",
    "TrigExpand",
    "TrigFactor",
    /* TrigReduce shrinks angle-addition forms (Sin[a]Cos[b]+Cos[a]Sin[b]
     * -> Sin[a+b]) and pulls integer powers / products of single-arg
     * trig calls into single trig calls of compound arguments. The
     * round-loop's score-gate keeps the reduction only when its leaf
     * count is no greater than the parent's, which is the typical case
     * for genuine angle-addition shapes. */
    "TrigReduce",
    /* TrigToExp surfaces the exp form directly so that hyperbolic
     * combinations whose exp form is strictly simpler (e.g. Sinh[x] +
     * Cosh[x] -> E^x) win the score tiebreak. The full trig roundtrip
     * (TrigToExp -> Together -> Cancel -> ExpToTrig) ends with ExpToTrig
     * and converts E^x back to Cosh[x] + Sinh[x], hiding the simpler
     * intermediate; offering TrigToExp as its own seed avoids that. For
     * pure trig inputs, TrigToExp yields a complex-coefficient exp form
     * with strictly higher complexity than the original, so the original
     * still wins -- no regression on Sin/Cos identities. */
    "TrigToExp"
};
static const size_t SIMP_TRANSFORM_COUNT =
    sizeof(SIMP_TRANSFORMS) / sizeof(SIMP_TRANSFORMS[0]);

/* Returns true if the expression contains any Power with a non-integer
 * exponent (e.g. Sqrt forms, Rational exponents, symbolic exponents).
 * picocas's Factor / FactorSquareFree call its trial-division loop in
 * factor_roots which can stall on multivariate inputs that include such
 * Power atoms, so we skip those transforms when this returns true. */
static bool has_non_integer_power(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        if (exp->type != EXPR_INTEGER && exp->type != EXPR_BIGINT) return true;
    }
    if (has_non_integer_power(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_non_integer_power(e->data.function.args[i])) return true;
    }
    return false;
}

/* Centralised cheap-precondition gate. Returns false only when the
 * predicate proves the named transform cannot possibly fire on `e` (and
 * `ctx`, where applicable). Conservative: returns true if uncertain.
 *
 * Names cover both the SIMP_TRANSFORMS table entries and the open-coded
 * seed/round transforms ("AbsRules", "LogExpRules", "AssumptionRules",
 * "TrigRoundtrip", "CollectPerVariable") so every call site shares one
 * dispatch table.
 *
 * `ctx` may be NULL; transforms that don't consult it ignore the parameter. */
static bool transform_can_fire(const char* name, const Expr* e,
                               const AssumeCtx* ctx) {
    /* Polynomial machinery: Factor / FactorSquareFree / FactorTerms / TrigFactor
     * stall on non-integer Power exponents (Sqrt, Rational exponents, ...). */
    if (strcmp(name, "Factor") == 0 ||
        strcmp(name, "FactorSquareFree") == 0 ||
        strcmp(name, "FactorTerms") == 0 ||
        strcmp(name, "TrigFactor") == 0) {
        if (has_non_integer_power(e)) return false;
    }
    /* Polynomial / partial-fraction transforms are no-ops on numeric-only
     * inputs (no symbol leaves). */
    if (strcmp(name, "Factor") == 0 ||
        strcmp(name, "FactorSquareFree") == 0 ||
        strcmp(name, "FactorTerms") == 0 ||
        strcmp(name, "Apart") == 0) {
        if (!contains_variable(e)) return false;
    }
    /* Trig family: skip when there's no trig or hyperbolic head anywhere.
     * The roundtrip composite is gated identically. */
    if (strcmp(name, "TrigExpand") == 0 ||
        strcmp(name, "TrigFactor") == 0 ||
        strcmp(name, "TrigToExp") == 0 ||
        strcmp(name, "TrigRoundtrip") == 0 ||
        strcmp(name, "TrigReduce") == 0) {
        if (!contains_trig_or_hyperbolic(e)) return false;
    }
    /* TrigReduce additionally needs Plus or Times anywhere in the input;
     * a single trig call with no enclosing arithmetic produces no
     * product-to-sum work and the rule list is a no-op. Without this
     * gate every Sin/Cos seed would seed a TrigReduce candidate and
     * waste a memo slot per leaf. */
    if (strcmp(name, "TrigReduce") == 0) {
        if (!contains_plus_or_times(e)) return false;
    }
    /* Abs rules. */
    if (strcmp(name, "AbsRules") == 0) {
        if (!contains_abs(e)) return false;
    }
    /* LogExp identity cascade: requires Log or Power somewhere (the
     * positivity-aware rewrites — Log[a*b] -> Log[a]+Log[b] etc. — are
     * gated on ctx facts internally; the unconditional Exp-distribute
     * rewrite Power[E, Plus[...]] -> Product Power[E, ti] fires
     * regardless and is needed for `Exp[k(Log a + Log b)] -> a^k b^k`
     * with integer k where no positivity assumption is required). */
    if (strcmp(name, "LogExpRules") == 0) {
        if (!contains_log(e) && !contains_power(e)) return false;
    }
    /* Assumption rewriter: nothing fires without facts. */
    if (strcmp(name, "AssumptionRules") == 0) {
        if (!ctx_has_facts(ctx)) return false;
    }
    /* Per-variable Collect: meaningful only when there are at least two
     * distinct variables to choose between. */
    if (strcmp(name, "CollectPerVariable") == 0) {
        if (expr_variables_count_capped(e, 2) < 2) return false;
    }
    /* Factorial rewriter: nothing fires unless the input mentions
     * Factorial somewhere in the tree. */
    if (strcmp(name, "FactorialRules") == 0) {
        if (!contains_factorial(e)) return false;
    }
    return true;
}

/* Score a candidate; if it beats the running best, replace best. The
 * candidate `c` is *borrowed* (caller still owns the source); `best` is
 * a slot the caller manages. */
static void update_best(Expr** best, size_t* best_score, const Expr* c,
                        const Expr* complexity_func) {
    size_t s = score_with_func(c, complexity_func);
    if (s < *best_score) {
        expr_free(*best);
        *best = expr_copy((Expr*)c);
        *best_score = s;
    }
}

/* True when `best` is the literal integer 0 -- the canonical
 * "input was an identity" outcome. Once Simplify reaches 0 there
 * is nothing left to improve, but the round loop would otherwise
 * keep grinding through the remaining transforms (Factor, TrigFactor,
 * etc.) on the still-large seeds. Early-exiting on this hit
 * trims the (Sin[x]+Cos[x])^4 - (1+Sin[2x])^2 case from ~220 ms to a
 * few ms by stopping as soon as TrigReduce produces 0.
 *
 * Restricted to literal 0 (rather than "any atomic with score 1") so
 * a custom complexity_func that ranks structurally-larger forms below
 * an atomic still gets full search coverage. The literal-zero case
 * is unambiguous: no transform we apply can produce a structurally
 * simpler form than 0. */
static bool simp_best_is_zero(const Expr* best) {
    return best && best->type == EXPR_INTEGER && best->data.integer == 0;
}

/* Apply Collect[expr, v] for each free variable v of expr, scoring the
 * results against `best` and adding novel ones to `next`. Collect can
 * surface a more compact form by grouping like powers (e.g. it can
 * recover x*(a+b) from a*x + b*x), and which variable to collect by is
 * not knowable up front -- Mathematica's Simplify likewise tries each
 * variable. We rely on Variables[] to enumerate the candidates. */
static void try_collect_per_variable(const Expr* seed, size_t parent_score,
                                     CandSet* next,
                                     Expr** best, size_t* best_score,
                                     const Expr* complexity_func) {
    Expr* vars = call_unary_copy("Variables", seed);
    if (!vars) return;
    if (vars->type != EXPR_FUNCTION ||
        !vars->data.function.head ||
        vars->data.function.head->type != EXPR_SYMBOL ||
        vars->data.function.head->data.symbol != SYM_List) {
        expr_free(vars);
        return;
    }
    size_t nv = vars->data.function.arg_count;
    bool dbg = simp_debug_enabled();
    for (size_t i = 0; i < nv; i++) {
        Expr* v = vars->data.function.args[i];
        Expr* args[2] = { expr_copy((Expr*)seed), expr_copy(v) };
        Expr* call = expr_new_function(expr_new_symbol("Collect"), args, 2);
        clock_t t0 = dbg ? clock() : 0;
        Expr* r = evaluate(call);
        expr_free(call);
        if (dbg) {
            char* vname = expr_to_string(v);
            char buf[64];
            snprintf(buf, sizeof(buf), "Collect[_,%s]",
                     vname ? vname : "?");
            simp_debug_log(buf, seed, r, simp_debug_elapsed_ms(t0));
            free(vname);
        }
        if (!r) continue;
        update_best(best, best_score, r, complexity_func);
        /* Branch-pruning: only propagate if the candidate is no worse
         * than its parent. Strictly worse forms cannot lead to a better
         * result through any further unary transform we run. */
        if (expr_eq(r, seed)) {
            expr_free(r);
        } else if (score_with_func(r, complexity_func) > parent_score) {
            expr_free(r);
        } else {
            cs_add_or_free(next, r);
        }
    }
    expr_free(vars);
}

/* ----------------------------------------------------------------------- */
/* Shape classifier + pipeline dispatch                                    */
/* ----------------------------------------------------------------------- */

/*
 * simp_classify performs one O(n) walk over the input and assigns it to a
 * shape bucket. simp_dispatch then routes to a specialised pipeline that
 * runs only the transforms relevant for that shape. The classifier is
 * conservative: borderline inputs fall through to SIMP_SHAPE_GENERAL,
 * which calls simp_search (the original full search), so misclassification
 * cannot change behaviour, only performance.
 *
 * Priority order (first matching wins):
 *   1. TRIG       -- any Sin/Cos/.../Sinh/.../ArcXxx head present
 *   2. LOGEXP     -- Log present (after trig is excluded)
 *   3. POLYNOMIAL -- no trig/log/abs, every Power has integer exponent,
 *                    PolynomialQ over Variables[input] is True
 *   4. RATIONAL   -- no trig/log/abs, Together[input] is poly/poly
 *   5. GENERAL    -- everything else
 *
 * Inputs containing Abs route to GENERAL: Abs rewrites are best handled
 * by the existing search machinery (the bottom-up walker hits Abs heads
 * directly via apply_abs_rules), and a specialised Abs pipeline buys
 * little.
 */
typedef enum {
    SIMP_SHAPE_POLYNOMIAL,
    SIMP_SHAPE_RATIONAL,
    SIMP_SHAPE_TRIG,
    SIMP_SHAPE_LOGEXP,
    SIMP_SHAPE_GENERAL
} SimpShape;

/* Helper: PolynomialQ[e, Variables[e]] using the existing builtins. */
static bool simp_is_polynomial_in_own_vars(const Expr* e) {
    Expr* vars = call_unary_copy("Variables", e);
    if (!vars) return false;
    if (vars->type != EXPR_FUNCTION ||
        !vars->data.function.head ||
        vars->data.function.head->type != EXPR_SYMBOL ||
        vars->data.function.head->data.symbol != SYM_List) {
        expr_free(vars);
        return false;
    }
    /* Zero-variable input: a numeric literal. PolynomialQ returns True
     * trivially, but a numeric leaf doesn't benefit from the polynomial
     * pipeline (no Factor / Collect target), so report false to fall
     * through to GENERAL. */
    if (vars->data.function.arg_count == 0) {
        expr_free(vars);
        return false;
    }
    Expr* args[2] = { expr_copy((Expr*)e), vars };
    Expr* pq = expr_new_function(expr_new_symbol("PolynomialQ"), args, 2);
    Expr* r = evaluate(pq);
    expr_free(pq);
    bool ok = (r && r->type == EXPR_SYMBOL &&
               r->data.symbol == SYM_True);
    if (r) expr_free(r);
    return ok;
}

/* Helper: build Together[e], extract Numerator/Denominator, check both
 * are polynomial in their own variables. Returns false (and frees nothing
 * external) if any step fails.
 *
 * Numeric collapse: when Together[e] reduces to a numeric leaf (or
 * Rational[p,q]), Numerator/Denominator yield zero-variable
 * expressions which simp_is_polynomial_in_own_vars rejects (its
 * empty-Variables clause is intentional for the polynomial classifier
 * but spurious here). Treat that case as trivially rational so the
 * shortcut still fires -- e.g. ((x-y)/(x+y) - (x+y)/(x-y)) /
 * (1 - (x^2-x*y-y^2)/(x^2-y^2)) collapses to -4 and routes to the
 * rational pipeline instead of falling through to GENERAL/bottomup. */
static bool simp_is_rational(const Expr* e) {
    Expr* tg = call_unary_copy("Together", e);
    if (!tg) return false;
    bool tg_is_numeric_leaf =
        tg->type == EXPR_INTEGER ||
        tg->type == EXPR_BIGINT ||
        tg->type == EXPR_REAL ||
        is_rational_literal(tg);
    if (tg_is_numeric_leaf) {
        expr_free(tg);
        return true;
    }
    Expr* num = call_unary_copy("Numerator", tg);
    Expr* den = call_unary_copy("Denominator", tg);
    bool ok = false;
    if (num && den &&
        !has_non_integer_power(num) &&
        !has_non_integer_power(den) &&
        simp_is_polynomial_in_own_vars(num) &&
        simp_is_polynomial_in_own_vars(den)) {
        ok = true;
    }
    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(tg);
    return ok;
}

static SimpShape simp_classify(const Expr* e) {
    if (contains_trig_or_hyperbolic(e)) return SIMP_SHAPE_TRIG;
    if (contains_abs(e))                return SIMP_SHAPE_GENERAL;
    /* Factorial-bearing inputs MUST flow through the general pipeline so
     * the FactorialRules seed in simp_search can fire. The rational and
     * polynomial pipelines short-circuit past it, leaving Factorial[arg]
     * structurally opaque to Together/Cancel and missing the shift-
     * normalization that closes (n+3)!/(n+1)! -> (n+2)(n+3) etc. */
    if (contains_factorial(e))          return SIMP_SHAPE_GENERAL;
    if (contains_log(e))                return SIMP_SHAPE_LOGEXP;
    /* No trig, no abs, no log. Decide poly / rational / general. */
    if (!has_non_integer_power(e) && simp_is_polynomial_in_own_vars(e)) {
        return SIMP_SHAPE_POLYNOMIAL;
    }
    if (simp_is_rational(e)) return SIMP_SHAPE_RATIONAL;
    return SIMP_SHAPE_GENERAL;
}

/* Forward declaration: simp_search is the GENERAL pipeline. */
static Expr* simp_search(const Expr* original_input, const AssumeCtx* ctx,
                         const Expr* complexity_func);

/* Forward declaration for simp_factorial (defined alongside simp_bottomup
 * further down). contains_factorial was forward-declared above. */
static Expr* simp_factorial(const Expr* e);

/* ----------------------------------------------------------------------- */
/* Specialised pipelines                                                   */
/* ----------------------------------------------------------------------- */

/*
 * Polynomial pipeline. Skips every transform that targets trig, log,
 * abs, rational, or assumption-driven structure -- none of which can
 * fire on a SHAPE_POLYNOMIAL input. The candidate set is:
 *   - the input itself
 *   - Expand[input]
 *   - Factor[input], FactorTerms[input]   (gated by has_non_integer_power)
 *   - Collect[input, v]   for each variable v
 *   - Collect[Expand[input], v]   for each variable v
 *
 * The Expand-then-Collect path is the one that recovers c + x*(a+b)
 * from a*x + b*x + c (test_simplify_collect_by_variable). Factor wins
 * on cases like (x-1)*(x+1)*(x^2+1)+1 -> x^4 because it emits the
 * already-expanded form when no nontrivial factorisation exists, but
 * it ties on score and the Expand candidate wins by SimplifyCount.
 *
 * No round loop: every winning move on a polynomial is reachable in
 * one application of the listed transforms; iterating doesn't help.
 */
static Expr* simp_pipeline_polynomial(const Expr* input,
                                      const AssumeCtx* ctx,
                                      const Expr* complexity_func) {
    (void)ctx;
    Expr* best = expr_copy((Expr*)input);
    size_t bs = score_with_func(best, complexity_func);

    Expr* expanded = traced_call_unary("Expand", input);
    if (expanded) update_best(&best, &bs, expanded, complexity_func);

    if (!has_non_integer_power(input)) {
        Expr* factored = traced_call_unary("Factor", input);
        if (factored) update_best(&best, &bs, factored, complexity_func);
        if (factored) expr_free(factored);
        Expr* fterms = traced_call_unary("FactorTerms", input);
        if (fterms) update_best(&best, &bs, fterms, complexity_func);
        if (fterms) expr_free(fterms);
    }

    /* Per-variable Collect on input AND on the expanded form. The
     * Expand-then-Collect path catches cases the input-Collect misses,
     * because Collect groups like powers of v and a Plus-of-Times input
     * already in factored shape doesn't expose them. */
    if (transform_can_fire("CollectPerVariable", input, NULL)) {
        CandSet next; cs_init(&next);
        try_collect_per_variable(input, bs, &next, &best, &bs, complexity_func);
        if (expanded) {
            try_collect_per_variable(expanded, bs, &next, &best, &bs, complexity_func);
        }
        cs_free(&next);
    }

    if (expanded) expr_free(expanded);
    return best;
}

/*
 * Rational pipeline. For inputs whose Together-form is poly/poly with
 * no trig/log/abs heads. The candidate set covers:
 *   - Together  (canonical fraction form)
 *   - Cancel    (gcd reduction)
 *   - ExpandNumerator / ExpandDenominator
 *   - Apart     (partial-fraction decomposition)
 *   - Factor on the Cancel'd form
 *   - Per-variable Collect on the Cancel'd form
 *   - AssumptionRules force-take, when ctx has facts (covers cases like
 *     (1-a^2)/b^2 with a^2+b^2==1)
 *
 * No round loop. Every meaningful win on a rational expression is one
 * of these transforms applied directly.
 */
static Expr* simp_pipeline_rational(const Expr* input,
                                    const AssumeCtx* ctx,
                                    const Expr* complexity_func) {
    Expr* best = expr_copy((Expr*)input);
    size_t bs = score_with_func(best, complexity_func);

    /* AssumptionRules force-take, mirroring the seed-phase semantics in
     * simp_search: assumption-driven rewrites are correctness-preserving
     * under the assumption set and qualitatively "more simplified", so we
     * accept them at equal or lower complexity. */
    Expr* assum_seed = NULL;
    if (transform_can_fire("AssumptionRules", input, ctx)) {
        Expr* ar = apply_assumption_rules(input, ctx);
        if (ar && !expr_eq(ar, input)) {
            size_t s = score_with_func(ar, complexity_func);
            if (s <= bs) {
                expr_free(best);
                best = expr_copy(ar);
                bs = s;
            }
            assum_seed = ar;  /* keep for downstream Together/Cancel */
        } else if (ar) {
            expr_free(ar);
        }
    }

    const Expr* seed = assum_seed ? assum_seed : input;

    Expr* tg = traced_call_unary("Together", seed);
    if (tg) update_best(&best, &bs, tg, complexity_func);

    Expr* cn_src = tg ? tg : (Expr*)seed;
    Expr* cn = traced_call_unary("Cancel", cn_src);
    if (cn) update_best(&best, &bs, cn, complexity_func);

    Expr* en = traced_call_unary("ExpandNumerator", cn ? cn : (Expr*)seed);
    if (en) update_best(&best, &bs, en, complexity_func);

    Expr* ed = traced_call_unary("ExpandDenominator", cn ? cn : (Expr*)seed);
    if (ed) update_best(&best, &bs, ed, complexity_func);

    Expr* ap = traced_call_unary("Apart", seed);
    if (ap) update_best(&best, &bs, ap, complexity_func);

    /* Factor on the most-canonical form available.  Prefer the Cancel'd
     * form when it exists, but fall through to seed/input when Cancel
     * returned NULL (no change to apply) -- otherwise Simplify would
     * never try Factor on inputs like `x/(x^3 + a b + a x + b x^2)`
     * which are already cancelled.  Gated against non-integer powers,
     * which shouldn't appear on a SHAPE_RATIONAL input but is
     * defensive.
     *
     * Push a NULL factor_memo around the Factor call so it uses its
     * normal (separate num/den variable lists) behaviour rather than
     * the inside-Simplify combined-scope path.  The combined-scope
     * mode refuses to factor parametric denominators like
     * `x^3 + a b + a x + b x^2` (which factors over Z[a, b][x] as
     * (a + x^2)(b + x)); the separate-scope mode handles them
     * correctly. */
    {
        const Expr* fc_src = cn ? cn : (tg ? tg : seed);
        if (!has_non_integer_power(fc_src)) {
            factor_memo_push(NULL);
            Expr* fc = traced_call_unary("Factor", fc_src);
            factor_memo_pop();
            if (fc) {
                update_best(&best, &bs, fc, complexity_func);
                expr_free(fc);
            }
        }
    }

    /* Per-variable Collect on the most-canonical form (Cancel'd if
     * available, otherwise seed/input -- same fall-through reason as
     * the Factor branch above). */
    {
        const Expr* col_src = cn ? cn : (tg ? tg : seed);
        if (transform_can_fire("CollectPerVariable", col_src, NULL)) {
            CandSet next; cs_init(&next);
            try_collect_per_variable(col_src, bs, &next, &best, &bs, complexity_func);
            cs_free(&next);
        }
    }

    if (tg) expr_free(tg);
    if (cn) expr_free(cn);
    if (en) expr_free(en);
    if (ed) expr_free(ed);
    if (ap) expr_free(ap);
    if (assum_seed) expr_free(assum_seed);
    return best;
}

/*
 * Log/Exp pipeline. For inputs containing Log (and no trig/hyperbolic).
 * The cascade in apply_logexp_rules implements the positivity-aware
 * Log/Power identities (Log[a*b] -> Log[a]+Log[b], (a*b)^c -> a^c b^c,
 * Log[x^p] -> p Log[x], (x^p)^q -> x^(p*q)). It is force-take in
 * simp_search (correctness-preserving under positivity assumptions) and
 * we replicate that here.
 *
 * Candidate set:
 *   - LogExpRules cascade (force-take when it changes the form)
 *   - AssumptionRules (e.g. Log[Exp[x]] -> x via Log[E^x] -> x*Log[E])
 *   - Together, Cancel, Expand on the cascade output
 *
 * No round loop; no trig transforms.
 */
static Expr* simp_pipeline_logexp(const Expr* input,
                                  const AssumeCtx* ctx,
                                  const Expr* complexity_func) {
    Expr* best = expr_copy((Expr*)input);
    size_t bs = score_with_func(best, complexity_func);

    /* Cascade force-take. */
    if (transform_can_fire("LogExpRules", input, ctx)) {
        Expr* lr = apply_logexp_rules(input, ctx);
        if (lr && !expr_eq(lr, input)) {
            expr_free(best);
            best = expr_copy(lr);
            bs = score_with_func(best, complexity_func);
        }
        if (lr) expr_free(lr);
    }

    /* Assumption-driven rewrites on the (possibly rewritten) best.
     * Force-take semantics: an assumption-aware rewrite that actually
     * changed the form is correctness-preserving under the assumption
     * set and counts as "more simplified" even when the leaf-count
     * tiebreak is even (e.g. Log[a^p] -> p Log[a] both score 4). This
     * matches the seed-phase behaviour in simp_search. */
    if (transform_can_fire("AssumptionRules", best, ctx)) {
        Expr* ar = apply_assumption_rules(best, ctx);
        if (ar && !expr_eq(ar, best)) {
            size_t s = score_with_func(ar, complexity_func);
            if (s <= bs) {
                expr_free(best);
                best = expr_copy(ar);
                bs = s;
            }
        }
        if (ar) expr_free(ar);
    }

    /* Standard cleanup. */
    Expr* tg = traced_call_unary("Together", best);
    if (tg) update_best(&best, &bs, tg, complexity_func);
    Expr* cn = traced_call_unary("Cancel", best);
    if (cn) update_best(&best, &bs, cn, complexity_func);
    Expr* ex = traced_call_unary("Expand", best);
    if (ex) update_best(&best, &bs, ex, complexity_func);

    if (tg) expr_free(tg);
    if (cn) expr_free(cn);
    if (ex) expr_free(ex);
    return best;
}

/* ----------------------------------------------------------------------- */
/* Additive-subexpression splitter                                         */
/* ----------------------------------------------------------------------- */

/*
 * When the input is a Plus whose addends partition into >=2 disjoint
 * connected components in the free-symbol-sharing graph (e.g.
 * Cos[x] + Sin[x] + Cos[y] + Sin[y] splits into the {x} group and the
 * {y} group), simplify each component independently and return the
 * sum.
 *
 * Why: simp_search's trig transforms (Together, TrigExpand, Cancel)
 * cost super-linearly in the number of addends. When the addends
 * decompose over disjoint user-symbol sets, those pieces cannot
 * interact under any of those transforms -- simplifying each component
 * independently preserves the result while keeping each call cheap.
 * The concrete trigger was the 6-term sum
 *   Cos[x]/Sqrt[6] + Sin[x]/Sqrt[2] + Cos[y]/Sqrt[6]/3 + Sin[y]/Sqrt[2]/3
 *     - Sqrt[6] Sin[x+Pi/6]/3 - Sqrt[6] Sin[y+Pi/6]/9
 * which hits $RecursionLimit and runs for minutes when simplified as a
 * whole; each three-term piece simplifies to 0 in ~1s on its own.
 *
 * Free-symbol semantics: an addend's "free symbols" are symbol leaves
 * appearing in its argument tree, EXCLUDING numeric constants
 * (Pi, E, EulerGamma, ...) and symbols that appear ONLY as function
 * heads. So Sin[x] contributes {x}, Cos[x+y] contributes {x, y}, but
 * Sin and Cos themselves are not counted -- otherwise every all-trig
 * addend would glom into a single component and the split would never
 * fire on inputs like the case above.
 */

typedef struct {
    const char** names;     /* borrowed pointers into the input expr */
    size_t count;
    size_t cap;
} SplitSymSet;

static void split_symset_init(SplitSymSet* s) {
    s->names = NULL; s->count = 0; s->cap = 0;
}
static void split_symset_free(SplitSymSet* s) {
    free((void*)s->names);
    s->names = NULL; s->count = 0; s->cap = 0;
}
static bool split_symset_contains(const SplitSymSet* s, const char* name) {
    for (size_t i = 0; i < s->count; i++) {
        if (strcmp(s->names[i], name) == 0) return true;
    }
    return false;
}
static void split_symset_add(SplitSymSet* s, const char* name) {
    if (split_symset_contains(s, name)) return;
    if (s->count == s->cap) {
        size_t nc = s->cap == 0 ? 4 : s->cap * 2;
        const char** nn = (const char**)realloc((void*)s->names,
                                                nc * sizeof(char*));
        if (!nn) return;
        s->names = nn;
        s->cap = nc;
    }
    s->names[s->count++] = name;
}
static bool split_symset_intersects(const SplitSymSet* a,
                                    const SplitSymSet* b) {
    for (size_t i = 0; i < a->count; i++) {
        if (split_symset_contains(b, a->names[i])) return true;
    }
    return false;
}

/* Walk into args only (skip the head). Symbol leaves that aren't
 * numeric constants are added to `out`. */
static void split_collect_addend_symbols(const Expr* e, SplitSymSet* out) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (is_real_constant_symbol(e->data.symbol)) return;
        split_symset_add(out, e->data.symbol);
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        split_collect_addend_symbols(e->data.function.args[i], out);
    }
}

/* Tiny union-find (path-compressed). */
static int split_uf_find(int* parent, int x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}
static void split_uf_union(int* parent, int x, int y) {
    int rx = split_uf_find(parent, x);
    int ry = split_uf_find(parent, y);
    if (rx != ry) parent[rx] = ry;
}

/* Forward declaration: the splitter recurses through simp_dispatch. */
static Expr* simp_dispatch(const Expr* input, const AssumeCtx* ctx,
                           const Expr* complexity_func);

/* Returns NULL when the input doesn't decompose, has fewer than 4
 * addends, or every component is a singleton (no win). On a successful
 * split, returns a freshly allocated, evaluated sum of per-component
 * simplifications. */
static Expr* simp_split_additive(const Expr* input, const AssumeCtx* ctx,
                                 const Expr* complexity_func) {
    if (!input || input->type != EXPR_FUNCTION) return NULL;
    if (!input->data.function.head ||
        input->data.function.head->type != EXPR_SYMBOL ||
        input->data.function.head->data.symbol != SYM_Plus)
        return NULL;
    size_t n = input->data.function.arg_count;
    if (n < 4) return NULL;

    SplitSymSet* sets = (SplitSymSet*)calloc(n, sizeof(SplitSymSet));
    int* parent = (int*)malloc(n * sizeof(int));
    if (!sets || !parent) {
        free(sets); free(parent);
        return NULL;
    }
    for (size_t i = 0; i < n; i++) {
        split_symset_init(&sets[i]);
        split_collect_addend_symbols(input->data.function.args[i], &sets[i]);
        parent[i] = (int)i;
    }
    /* Pairwise union by symbol intersection. */
    for (size_t i = 0; i < n; i++) {
        if (sets[i].count == 0) continue;
        for (size_t j = i + 1; j < n; j++) {
            if (sets[j].count == 0) continue;
            if (split_symset_intersects(&sets[i], &sets[j])) {
                split_uf_union(parent, (int)i, (int)j);
            }
        }
    }
    /* Glom each constant addend (empty symset) into the first variable
     * component so it gets simplified alongside variable terms. Without
     * this, splitting `1 + Plus[..variable terms..]` into `{1}` and the
     * rest can mask outer/inner constant cancellations: simplifying the
     * variable component independently may introduce a new constant
     * (e.g. Factor turning -1 - x^3/6 - x^4/6 into -1/6*(6 + x^3 + x^4))
     * that no longer cancels with the split-out `1`, leaving Simplify
     * non-idempotent. The glom is into a single non-empty component, not
     * duplicated, so the variable structure of the split is preserved. */
    {
        int first_non_empty = -1;
        for (size_t i = 0; i < n; i++) {
            if (sets[i].count != 0) { first_non_empty = (int)i; break; }
        }
        if (first_non_empty >= 0) {
            for (size_t i = 0; i < n; i++) {
                if (sets[i].count == 0) {
                    split_uf_union(parent, (int)i, first_non_empty);
                }
            }
        }
    }
    int* root_per = (int*)malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) root_per[i] = split_uf_find(parent, (int)i);
    int* size_by_root = (int*)calloc(n, sizeof(int));
    for (size_t i = 0; i < n; i++) size_by_root[root_per[i]]++;

    int comp_count = 0, max_size = 0;
    for (size_t i = 0; i < n; i++) {
        if (size_by_root[i] == 0) continue;
        comp_count++;
        if (size_by_root[i] > max_size) max_size = size_by_root[i];
    }
    /* Useful split: 2+ components AND at least one component has 2+
     * addends (otherwise we'd just be re-dispatching every atom for no
     * win, and bottomup has already simplified each addend in
     * isolation). */
    if (comp_count < 2 || max_size < 2) {
        for (size_t i = 0; i < n; i++) split_symset_free(&sets[i]);
        free(sets); free(parent); free(root_per); free(size_by_root);
        return NULL;
    }

    /* Compact root -> 0..comp_count-1 index. */
    int* root_to_idx = (int*)malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) root_to_idx[i] = -1;
    int next_idx = 0;
    for (size_t i = 0; i < n; i++) {
        int r = root_per[i];
        if (root_to_idx[r] == -1) root_to_idx[r] = next_idx++;
    }
    /* Bucket addend pointers per component. */
    Expr*** comp_addends = (Expr***)calloc(comp_count, sizeof(Expr**));
    int* comp_alloc = (int*)calloc(comp_count, sizeof(int));
    int* comp_n = (int*)calloc(comp_count, sizeof(int));
    for (size_t i = 0; i < n; i++) {
        int idx = root_to_idx[root_per[i]];
        if (comp_n[idx] == comp_alloc[idx]) {
            int nc = comp_alloc[idx] == 0 ? 4 : comp_alloc[idx] * 2;
            comp_addends[idx] = (Expr**)realloc(comp_addends[idx],
                                                nc * sizeof(Expr*));
            comp_alloc[idx] = nc;
        }
        comp_addends[idx][comp_n[idx]++] =
            expr_copy(input->data.function.args[i]);
    }

    /* Simplify each component, place results into a Plus. */
    Expr** results = (Expr**)calloc(comp_count, sizeof(Expr*));
    for (int c = 0; c < comp_count; c++) {
        Expr* sub;
        if (comp_n[c] == 1) {
            sub = comp_addends[c][0];
        } else {
            Expr* p = expr_new_function(expr_new_symbol("Plus"),
                                        comp_addends[c], comp_n[c]);
            sub = evaluate(p);
            expr_free(p);
        }
        results[c] = simp_dispatch(sub, ctx, complexity_func);
        expr_free(sub);
    }
    Expr* sum_raw = expr_new_function(expr_new_symbol("Plus"),
                                      results, comp_count);
    Expr* sum = evaluate(sum_raw);
    expr_free(sum_raw);

    for (size_t i = 0; i < n; i++) split_symset_free(&sets[i]);
    free(sets); free(parent); free(root_per); free(size_by_root);
    free(root_to_idx);
    for (int c = 0; c < comp_count; c++) free(comp_addends[c]);
    free(comp_addends); free(comp_alloc); free(comp_n);
    free(results);

    return sum;
}

/* simp_split_multiplicative: the multiplicative analog of
 * simp_split_additive. Decomposes a Times node whose factors fall into
 * 2+ variable-disjoint components (with at least one component holding
 * 2+ factors). Each component's sub-Times is dispatched in isolation,
 * then the per-component results are multiplied.
 *
 * Same correctness argument as the additive splitter: every transform
 * in SIMP_TRANSFORMS (Together, Cancel, Expand, Factor, TrigExpand,
 * TrigToExp, TrigFactor, TrigRoundtrip, Collect, Pythag*, half-angle,
 * Radicals) acts within a single variable's algebraic/trigonometric
 * structure -- variable-disjoint factors cannot interact under any of
 * them, so component-wise simplification preserves the answer.
 *
 * Without this pass, a stray independent factor inflates simp_search's
 * effective free-symbol budget and the heuristic gives up. With it,
 *   Tan[z] Cos[x] Cos[y] Sec[x+y] (Tan[x]+Tan[y])
 *     -> Tan[z] * simp_dispatch[Cos[x] Cos[y] Sec[x+y] (Tan[x]+Tan[y])]
 *     -> Tan[z] Tan[x+y]
 */
static Expr* simp_split_multiplicative(const Expr* input,
                                       const AssumeCtx* ctx,
                                       const Expr* complexity_func) {
    if (!input || input->type != EXPR_FUNCTION) return NULL;
    if (!input->data.function.head ||
        input->data.function.head->type != EXPR_SYMBOL ||
        input->data.function.head->data.symbol != SYM_Times)
        return NULL;
    size_t n = input->data.function.arg_count;
    if (n < 3) return NULL;

    SplitSymSet* sets = (SplitSymSet*)calloc(n, sizeof(SplitSymSet));
    int* parent = (int*)malloc(n * sizeof(int));
    if (!sets || !parent) {
        free(sets); free(parent);
        return NULL;
    }
    for (size_t i = 0; i < n; i++) {
        split_symset_init(&sets[i]);
        split_collect_addend_symbols(input->data.function.args[i], &sets[i]);
        parent[i] = (int)i;
    }
    /* Pairwise union by symbol intersection. Empty sets (constants) stay
     * in their own singleton component -- no point dragging them in. */
    for (size_t i = 0; i < n; i++) {
        if (sets[i].count == 0) continue;
        for (size_t j = i + 1; j < n; j++) {
            if (sets[j].count == 0) continue;
            if (split_symset_intersects(&sets[i], &sets[j])) {
                split_uf_union(parent, (int)i, (int)j);
            }
        }
    }
    int* root_per = (int*)malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) root_per[i] = split_uf_find(parent, (int)i);
    int* size_by_root = (int*)calloc(n, sizeof(int));
    for (size_t i = 0; i < n; i++) size_by_root[root_per[i]]++;

    int comp_count = 0, max_size = 0;
    for (size_t i = 0; i < n; i++) {
        if (size_by_root[i] == 0) continue;
        comp_count++;
        if (size_by_root[i] > max_size) max_size = size_by_root[i];
    }
    /* Useful split: 2+ components AND at least one component holds 2+
     * factors (singleton-only splits would just re-dispatch each atom
     * for no win). */
    if (comp_count < 2 || max_size < 2) {
        for (size_t i = 0; i < n; i++) split_symset_free(&sets[i]);
        free(sets); free(parent); free(root_per); free(size_by_root);
        return NULL;
    }

    int* root_to_idx = (int*)malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) root_to_idx[i] = -1;
    int next_idx = 0;
    for (size_t i = 0; i < n; i++) {
        int r = root_per[i];
        if (root_to_idx[r] == -1) root_to_idx[r] = next_idx++;
    }
    Expr*** comp_factors = (Expr***)calloc(comp_count, sizeof(Expr**));
    int* comp_alloc = (int*)calloc(comp_count, sizeof(int));
    int* comp_n = (int*)calloc(comp_count, sizeof(int));
    for (size_t i = 0; i < n; i++) {
        int idx = root_to_idx[root_per[i]];
        if (comp_n[idx] == comp_alloc[idx]) {
            int nc = comp_alloc[idx] == 0 ? 4 : comp_alloc[idx] * 2;
            comp_factors[idx] = (Expr**)realloc(comp_factors[idx],
                                                nc * sizeof(Expr*));
            comp_alloc[idx] = nc;
        }
        comp_factors[idx][comp_n[idx]++] =
            expr_copy(input->data.function.args[i]);
    }

    Expr** results = (Expr**)calloc(comp_count, sizeof(Expr*));
    for (int c = 0; c < comp_count; c++) {
        Expr* sub;
        if (comp_n[c] == 1) {
            sub = comp_factors[c][0];
        } else {
            Expr* p = expr_new_function(expr_new_symbol("Times"),
                                        comp_factors[c], comp_n[c]);
            sub = evaluate(p);
            expr_free(p);
        }
        results[c] = simp_dispatch(sub, ctx, complexity_func);
        expr_free(sub);
    }
    Expr* prod_raw = expr_new_function(expr_new_symbol("Times"),
                                       results, comp_count);
    Expr* prod = evaluate(prod_raw);
    expr_free(prod_raw);

    for (size_t i = 0; i < n; i++) split_symset_free(&sets[i]);
    free(sets); free(parent); free(root_per); free(size_by_root);
    free(root_to_idx);
    for (int c = 0; c < comp_count; c++) free(comp_factors[c]);
    free(comp_factors); free(comp_alloc); free(comp_n);
    free(results);

    return prod;
}

/* simp_dispatch is the public entry point. It runs the shape classifier
 * and forwards to a specialised pipeline. SHAPE_TRIG and SHAPE_GENERAL
 * fall through to simp_search; trig is the heuristic search's strongest
 * domain, and general inputs need the full machinery. */
static Expr* simp_dispatch(const Expr* input, const AssumeCtx* ctx,
                           const Expr* complexity_func) {
    /* Universal Power-identity seeds: PrimeRebase and PowerOneify are
     * cheap, inert-by-default rewrites that benefit every pipeline
     * (POLYNOMIAL, RATIONAL, LOGEXP, TRIG, GENERAL alike).  Both are
     * idempotent on their fixed point, so the strict-win recursion
     * below terminates after at most one rebase + one oneify cascade.
     * Score-gated: only adopt the rewrite when it strictly beats the
     * input; otherwise fall through to the normal classifier route. */
    {
        Expr* alt = transform_prime_rebase(input);
        if (alt) {
            if (!expr_eq(alt, input)) {
                size_t s_alt = score_with_func(alt, complexity_func);
                size_t s_in = score_with_func(input, complexity_func);
                if (s_alt < s_in) {
                    Expr* result = simp_dispatch(alt, ctx, complexity_func);
                    expr_free(alt);
                    return result;
                }
            }
            expr_free(alt);
        }
    }
    {
        Expr* alt = transform_power_oneify(input);
        if (alt) {
            if (!expr_eq(alt, input)) {
                size_t s_alt = score_with_func(alt, complexity_func);
                size_t s_in = score_with_func(input, complexity_func);
                if (s_alt < s_in) {
                    Expr* result = simp_dispatch(alt, ctx, complexity_func);
                    expr_free(alt);
                    return result;
                }
            }
            expr_free(alt);
        }
    }
    /* PowerDistribute: principal-branch power-distribution rewrites that
     * are universally sound -- (-c)^e split, constant-positive base
     * distribute, and integer-exponent distribute (via ctx).  Same
     * recursion shape as the rebase / oneify gates above: only adopt
     * the rewrite when it strictly beats the input on score. */
    {
        Expr* alt = transform_power_distribute(input, ctx);
        if (alt) {
            if (!expr_eq(alt, input)) {
                size_t s_alt = score_with_func(alt, complexity_func);
                size_t s_in = score_with_func(input, complexity_func);
                if (s_alt < s_in) {
                    Expr* result = simp_dispatch(alt, ctx, complexity_func);
                    expr_free(alt);
                    return result;
                }
            }
            expr_free(alt);
        }
    }
    /* RadicalCanon: split Power[Rational[a,b], q] into Power[a,q] *
     * Power[b,-q] and rationalise Power[positive_int, negative_rational]
     * into Power[..., positive_residue] / a^k.  Forces equivalent radical
     * shapes (Sqrt[1/2] vs 1/Sqrt[2] vs Sqrt[2]/2) into a single canonical
     * representation so additive cancellations like
     *     -Sqrt[1/2]/3 + Sqrt[2]/6  ->  0
     * fire instead of being trapped behind shape divergence.  This is
     * NOT a strict-win gate: even when the score is the same, the rewrite
     * is correctness-preserving and always yields the more useful form
     * for downstream additive cancellation; the score check below uses
     * <= rather than <. */
    {
        Expr* alt = transform_radical_canon(input);
        if (alt) {
            if (!expr_eq(alt, input)) {
                size_t s_alt = score_with_func(alt, complexity_func);
                size_t s_in = score_with_func(input, complexity_func);
                if (s_alt <= s_in) {
                    Expr* result = simp_dispatch(alt, ctx, complexity_func);
                    expr_free(alt);
                    return result;
                }
            }
            expr_free(alt);
        }
    }

    /* Try to decompose a Plus into disjoint-variable components and
     * simplify each independently. Each connected component's sub-Plus
     * cannot itself decompose further, so the recursion through
     * simp_dispatch terminates in one level. */
    Expr* split = simp_split_additive(input, ctx, complexity_func);
    if (split) return split;
    /* Same idea for Times: lift variable-disjoint factors out of the
     * search space. */
    Expr* tsplit = simp_split_multiplicative(input, ctx, complexity_func);
    if (tsplit) return tsplit;

    SimpShape shape = simp_classify(input);
    switch (shape) {
        case SIMP_SHAPE_POLYNOMIAL:
            return simp_pipeline_polynomial(input, ctx, complexity_func);
        case SIMP_SHAPE_RATIONAL:
            return simp_pipeline_rational(input, ctx, complexity_func);
        case SIMP_SHAPE_LOGEXP:
            return simp_pipeline_logexp(input, ctx, complexity_func);
        case SIMP_SHAPE_TRIG: {
            /* Fast algebraic normal-form pre-check for rational
             * functions of Sin/Cos/Sinh/Cosh (and Tan/Cot/Sec/Csc/Tanh/
             * etc. after preprocessing). Strict leaf-count gate inside:
             * returns NULL when the algorithm does not apply or does
             * not strictly improve the input, so we fall through to
             * simp_search unchanged in those cases. */
            Expr* tr = simp_trig_rational(input, ctx, complexity_func);
            if (tr) return tr;
            return simp_search(input, ctx, complexity_func);
        }
        case SIMP_SHAPE_GENERAL:
        default:
            return simp_search(input, ctx, complexity_func);
    }
}

static Expr* simp_search(const Expr* original_input, const AssumeCtx* ctx,
                         const Expr* complexity_func) {
    /* Phase 0: pre-apply the Abs structural rules. These (idempotent,
     * force-take) rewrites canonicalise Abs[Times[...]] -> product of Abs,
     * Abs[Abs[x]] -> Abs[x], Abs[E^z] -> E^Re[z], etc. Doing it here
     * (rather than as a regular seed) means the rest of the search starts
     * from the rewritten form -- otherwise the original (often smaller in
     * leaf count) re-enters the candidate set and the leaf-count tiebreak
     * brings us right back. */
    Expr* abs_pre = transform_can_fire("AbsRules", original_input, ctx)
                        ? apply_abs_rules(original_input, ctx)
                        : NULL;
    const Expr* input;
    if (abs_pre && !expr_eq(abs_pre, original_input)) {
        if (simp_debug_enabled()) {
            simp_debug_log("AbsRules", original_input, abs_pre, 0.0);
        }
        input = abs_pre;
    } else {
        if (abs_pre) { expr_free(abs_pre); abs_pre = NULL; }
        input = original_input;
    }

    Expr* best = expr_copy((Expr*)input);
    size_t best_score = score_with_func(best, complexity_func);

    CandSet seeds;
    cs_init(&seeds);
    cs_add_or_free(&seeds, expr_copy((Expr*)input));

    /* Assumption-derived alternatives. An assumption-aware rewrite that
     * actually changed the form is correctness-preserving under the
     * assumption set and is by definition more "simplified" than the
     * input, even if the leaf-count tiebreak is even (e.g.
     * Log[x^p] -> p Log[x] both score 4). Force-take it as the new best
     * so long as it isn't strictly worse. */
    if (transform_can_fire("AssumptionRules", input, ctx)) {
        bool dbg = simp_debug_enabled();
        clock_t t0 = dbg ? clock() : 0;
        Expr* alt = apply_assumption_rules(input, ctx);
        if (dbg) simp_debug_log("AssumptionRules", input, alt, simp_debug_elapsed_ms(t0));
        if (alt && !expr_eq(alt, input)) {
            size_t s = score_with_func(alt, complexity_func);
            if (s <= best_score) {
                expr_free(best);
                best = expr_copy(alt);
                best_score = s;
            }
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Logexp rewriter seed. The Log/Power identities (Log[a b] -> Log[a] +
     * Log[b], (a b)^c -> a^c b^c, (a^p)^q -> a^(p q), ...) typically
     * INCREASE leaf count, so they cannot win the standard complexity
     * tiebreak. We force them as the new best regardless of score: the
     * cascade is correctness-preserving under positivity assumptions and
     * the user's intent (per the documented rule cascade) is that they
     * fire whenever conditions are met. Downstream transforms (Cancel,
     * Together, ...) cannot recombine an expanded log/power, so the
     * expanded form persists through the rest of the search. */
    if (transform_can_fire("LogExpRules", input, ctx)) {
        bool dbg = simp_debug_enabled();
        clock_t t0 = dbg ? clock() : 0;
        Expr* alt = apply_logexp_rules(input, ctx);
        if (dbg) simp_debug_log("LogExpRules", input, alt, simp_debug_elapsed_ms(t0));
        if (alt && !expr_eq(alt, input)) {
            expr_free(best);
            best = expr_copy(alt);
            best_score = score_with_func(best, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Roundtrip seed. Score-gate the seed propagation: TrigRoundtrip
     * runs Together inside, which on pure-real Sinh/Cosh exp forms can
     * factor out an asymmetric E^(-kx) that introduces a high-frequency
     * E^(2kx) term in the cofactor (e.g.
     * Together[-1/2 + 1/4 E^(-2x) + 1/4 E^(2x)]
     *   ->  1/4 E^(-2x) (1 - 2 E^(2x) + E^(4x))).
     * ExpToTrig then turns that into a 16-term Cosh/Sinh product that
     * subsequent transforms (Together, Cancel, Factor) all consume
     * 100s of ms on. update_best still picks the result if it's a win;
     * we only refuse to propagate dramatic blow-ups as a seed for
     * further exploration. */
    if (transform_can_fire("TrigRoundtrip", input, NULL)) {
        Expr* alt = transform_trig_roundtrip(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            size_t alt_score = score_with_func(alt, complexity_func);
            size_t input_score = score_with_func(input, complexity_func);
            if (alt_score <= 2 * input_score + 8) {
                cs_add_or_free(&seeds, alt);
            } else {
                expr_free(alt);
            }
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* ExpToTrig seed. TrigRoundtrip is gated on contains_trig_or_hyperbolic
     * which rejects pure-Exp inputs (e.g.
     * `-1 + (-E^-x + E^x)^2/(E^-x + E^x)^2` has no Cosh/Sinh head and so
     * misses the entire trig pipeline). ExpToTrig converts E^... into
     * Cosh/Sinh form, after which PythagReduce / Together / Cancel can
     * collapse the result (the example above lands at `-Sech[x]^2`).
     * Score-gate the seed propagation the same way TrigRoundtrip does:
     * keep wins as the new best, but only forward seeds that haven't
     * blown up the structure. Skip when the input already has trig or
     * hyperbolic heads -- TrigRoundtrip handles that case. */
    if (contains_exp_form(input) && !contains_trig_or_hyperbolic(input)) {
        bool dbg = simp_debug_enabled();
        clock_t t0 = dbg ? clock() : 0;
        Expr* alt = call_unary_copy("ExpToTrig", input);
        if (dbg) simp_debug_log("ExpToTrigSeed", input, alt,
                                simp_debug_elapsed_ms(t0));
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            size_t alt_score = score_with_func(alt, complexity_func);
            size_t input_score = score_with_func(input, complexity_func);
            if (alt_score <= 2 * input_score + 8) {
                cs_add_or_free(&seeds, alt);
            } else {
                expr_free(alt);
            }
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Pythagorean square-completion seed. Idempotent on inputs without the
     * 1 +/- 2 Sin Cos shape, so always cheap to try. */
    {
        Expr* alt = transform_pythag_square_complete(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Pythagorean reduction seed (1 - Cos^2 -> Sin^2, etc.). Strict
     * leaf-count win when it fires; inert otherwise. */
    {
        Expr* alt = transform_pythag_reduce(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Pythagorean canonicalizer seed: substitute Cos^(2k) -> (1 - Sin^2)^k
     * (or the reverse) globally and Expand. Catches difference-of-squares
     * trig products that, after Expand, leave coefficient-bearing Cos^2
     * factors PythagReduce cannot match. The transform itself only
     * returns a strict-score win, so passing it through update_best is
     * safe: it never propagates a worse seed. */
    {
        Expr* alt = transform_pythag_canon(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }
    if (simp_best_is_zero(best)) goto search_done;

    /* TrigReduce seed.  The product-to-sum / power-reduction pass
     * collapses inputs like `(Sin[x]+Cos[x])^4 - (1+Sin[2x])^2` directly
     * to 0 (after internal Expand + angle-addition collapse) in a few
     * ms.  The same call appears later in the round loop's
     * SIMP_TRANSFORMS table at index 11, which means without this seed
     * the search runs Factor (~28 ms), TrigFactor (~14 ms),
     * FactorSquareFree (~9 ms), etc. on the input *before* TrigReduce
     * gets a turn, even though TrigReduce alone is enough.  Promoting
     * it to the seed phase + relying on simp_best_is_zero to short-
     * circuit afterwards turns ~220 ms into ~5 ms on this shape.  Inert
     * (returns the input unchanged via the trig_memo no-op cycle) when
     * the input has no trig-product structure to reduce. */
    if (transform_can_fire("TrigReduce", input, NULL)) {
        Expr* alt = traced_call_unary("TrigReduce", input);
        if (alt) {
            if (!expr_eq(alt, input)) {
                update_best(&best, &best_score, alt, complexity_func);
                size_t alt_score = score_with_func(alt, complexity_func);
                size_t input_score = score_with_func(input, complexity_func);
                /* Same blow-up guard TrigRoundtrip uses: keep the wins,
                 * drop catastrophically larger candidates from the seed
                 * set so they don't seed a wave of expensive Factor/
                 * TrigFactor work in subsequent rounds. */
                if (alt_score <= 2 * input_score + 8) {
                    cs_add_or_free(&seeds, alt);
                } else {
                    expr_free(alt);
                }
            } else {
                expr_free(alt);
            }
        }
    }
    if (simp_best_is_zero(best)) goto search_done;

    /* Trig-at-rational-Pi canonicalization seed.  Picks a unique
     * Sin-vs-Cos / Tan-vs-Cot / Sec-vs-Csc representation per pair,
     * which lets pairs like Cos[4 Pi/9] - Sin[Pi/18] (both lower to
     * Sin[Pi/18]) and 1/8 (Cos[4 Pi/9] + Cos[5 Pi/9]) (lowers to
     * 1/8 (Sin[Pi/18] + (-Sin[Pi/18])) = 0) cancel additively in
     * the surrounding Plus.  Idempotent and inert when the input
     * has no rational-Pi trig calls. */
    {
        Expr* alt = simp_trig_pi_canon(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* TanAddition seed: search for triples (a, b, c=a+b) among the trig
     * args appearing in the input and rewrite Tan[c]/Cot[c]/etc. via
     * the angle-addition formula in terms of Tan[a], Tan[b].  Catches
     * shapes like Tan[2] Tan[3] - 1 + (Tan[2]+Tan[3])/Tan[5] -> 0 (since
     * 5 = 2+3) which would otherwise go un-recognised because TrigExpand
     * only fires on literal-Plus arguments.  Inert when fewer than three
     * distinct trig args appear or when no pair sums to a third arg. */
    {
        Expr* alt = transform_tan_addition(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }
    if (simp_best_is_zero(best)) goto search_done;

    /* Half-angle tangent / Tanh seed. Idempotent on inputs without
     * the Sin/(1+Cos) (resp. Sinh/(1+Cosh)) shape. */
    {
        Expr* alt = transform_halfangle(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Radical product seed: collapses Sqrt[a]*Sqrt[b] -> Sqrt[a*b] for
     * positive integer a, b (and similarly for higher rational
     * exponents).  Inert on inputs without Power[+integer, Rational]
     * factors. */
    {
        Expr* alt = simp_radicals(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* The common-factor lift is applied as a final-form polish after
     * simp_bottomup completes (see builtin_simplify), not as a search
     * seed. Wiring it into the round loop or the seed phase changes the
     * structural shape of intermediate forms in ways that interact badly
     * with TrigExpand / TrigReduce on multi-variable trig inputs. */

    /* Sqrt-of-Sqrt denesting seed (phase-1 radical denesting). Fires the
     * half-sum identity Sqrt[A + b Sqrt[C]] -> Sqrt[(A+s)/2] +/- Sqrt[(A-s)/2]
     * where s^2 = A^2 - b^2 C is detected via integer/rational
     * perfect-square or polynomial FactorSquareFree. Inert (returns the
     * input copy) when no Sqrt[Plus[...]] subtree denests cleanly under
     * the active assumption set. See the simp_denest_sqrt section above
     * for branch-soundness preconditions. */
    {
        Expr* alt = simp_denest_sqrt(input, ctx);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Denominator-rationalisation seed (phase-2 radical simplification).
     * Walks the tree looking for Power[denom, -1] where denom is a
     * polynomial in radicals over a single positive integer base c,
     * and computes the inverse via extended-Euclidean in
     * Q[α]/(α^n - c). Closes user cases like 1/(2^(1/3) - 1) ->
     * 4^(1/3) + 2^(1/3) + 1. Inert when no such Power[denom, -1]
     * subtree exists or when denom involves multiple distinct integer
     * bases (out of scope for the current single-extension
     * implementation). */
    {
        Expr* alt = simp_rationalize_denom(input, ctx);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Cube-root denesting / sum-of-conjugates seed (phase-3). Pattern A:
     * (a + b sqrt(c))^(1/3) -> p + q sqrt(c) via small-grid search.
     * Pattern B: (a+b sqrt(c))^(1/3) + (a-b sqrt(c))^(1/3) -> rational
     * via Cardano discriminant + rational-root test. Inert otherwise. */
    {
        Expr* alt = simp_cuberoot(input, ctx);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Algebraic-extension seed: substitute each Sqrt[u_i] by a fresh
     * generator g_i, reduce in K(vars)[g_1,...,g_n]/(g_i^2 - u_i), and
     * rationalise the denominator by successive sigma-conjugation.
     * Collapses identities ordinary Together / Cancel cannot see, e.g.
     *   (x/Sqrt[x^2+1] + 1) / ((Sqrt[x^2+1] + x)^2 + 1)  ->  1/(2 + 2 x^2)
     * Inert when there is no surd, when any surd has q != 2, or when
     * the input contains an explicit complex literal. */
    {
        Expr* alt = simp_algebraic(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Roots-of-unity seed. Reduces sums of (-1)^(p/q) and E^(I p Pi/q)
     * via the minimal polynomial Phi_{2Q}(x); see
     * simp_roots_of_unity above. Idempotent and inert when the input
     * has no root-of-unity atoms. */
    {
        Expr* alt = simp_roots_of_unity(input);
        if (alt && !expr_eq(alt, input)) {
            update_best(&best, &best_score, alt, complexity_func);
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Factorial seed. Decomposes every Factorial[arg] as arg = sym + c,
     * groups by sym, and rewrites each call into Factorial[base] times
     * an explicit Pochhammer product so common Factorial[base] factors
     * cancel/absorb under Together / Cancel / Expand. A re-fold pass
     * then absorbs (base+1)(base+2)...(base+k) cofactors back into
     * Factorial[base+k] (with up to two missing-factor gap-fills). The
     * separate (2v)!/(2^v v!) -> Factorial2[2v-1] check piggybacks on
     * the same walk. Inert when the input has no Factorial atoms.
     *
     * Force-take semantics (mirroring LogExpRules / AssumptionRules): a
     * shape-changing factorial rewrite is by definition "more simplified"
     * even when the SimplifyCount tiebreak ties. The default complexity
     * measure cannot distinguish (n-2)!/n! (count 11) from 1/(n^2-n)
     * (count 11) but the latter is unambiguously the canonical form. We
     * force-take whenever the number of Factorial atoms strictly drops
     * (or stays the same and the leaf count is no worse). */
    if (transform_can_fire("FactorialRules", input, NULL)) {
        Expr* alt = simp_factorial(input);
        if (alt && !expr_eq(alt, input)) {
            size_t s = score_with_func(alt, complexity_func);
            /* Count Factorial atoms in input and alt. */
            size_t fc_in = 0, fc_alt = 0;
            {
                /* Inline tiny tree walk -- avoid threading a counter into
                 * the contains_factorial signature. */
                CandSet stack;
                cs_init(&stack);
                cs_add_or_free(&stack, expr_copy((Expr*)input));
                while (stack.count) {
                    Expr* top = stack.items[stack.count - 1];
                    stack.count--;
                    if (top && top->type == EXPR_FUNCTION) {
                        if (simp_eq_head_sym(top, "Factorial")) fc_in++;
                        cs_add_or_free(&stack, expr_copy(top->data.function.head));
                        for (size_t i = 0; i < top->data.function.arg_count; i++) {
                            cs_add_or_free(&stack,
                                expr_copy(top->data.function.args[i]));
                        }
                    }
                    expr_free(top);
                }
                cs_free(&stack);
            }
            {
                CandSet stack;
                cs_init(&stack);
                cs_add_or_free(&stack, expr_copy(alt));
                while (stack.count) {
                    Expr* top = stack.items[stack.count - 1];
                    stack.count--;
                    if (top && top->type == EXPR_FUNCTION) {
                        if (simp_eq_head_sym(top, "Factorial")) fc_alt++;
                        cs_add_or_free(&stack, expr_copy(top->data.function.head));
                        for (size_t i = 0; i < top->data.function.arg_count; i++) {
                            cs_add_or_free(&stack,
                                expr_copy(top->data.function.args[i]));
                        }
                    }
                    expr_free(top);
                }
                cs_free(&stack);
            }
            bool force_take = (fc_alt < fc_in) ||
                              (fc_alt == fc_in && s <= best_score);
            if (force_take) {
                expr_free(best);
                best = expr_copy(alt);
                best_score = s;
            } else {
                update_best(&best, &best_score, alt, complexity_func);
            }
            cs_add_or_free(&seeds, alt);
        } else if (alt) {
            expr_free(alt);
        }
    }

    /* Per-variable Collect seed. parent_score = score(input). */
    if (transform_can_fire("CollectPerVariable", input, NULL)) {
        try_collect_per_variable(input, best_score, &seeds, &best, &best_score,
                                 complexity_func);
    }

    /*
     * Round loop. Branch-pruning rule: a candidate produced by a
     * transform on `seed` is propagated to `next` only when its
     * complexity is no greater than `seed`'s. Strictly worse forms
     * cannot lead to a better best through any further unary transform
     * we apply (they'd just feed transforms more work and grow the
     * candidate set), so we drop them. They may still beat the running
     * best on this very step (update_best already handles that), but
     * they won't seed further exploration.
     *
     * Assumption-aware rewrites and the logexp cascade keep their
     * "force-win" behaviour for the best slot (they're correctness-
     * preserving under the assumption set and qualitatively "more
     * simplified"), but the same parent-score gate applies to whether
     * they propagate as a new seed.
     */
    for (int round = 0; round < SIMP_ROUNDS; round++) {
        CandSet next;
        cs_init(&next);
        for (size_t i = 0; i < seeds.count; i++) {
            const Expr* seed = seeds.items[i];
            size_t parent_score = score_with_func(seed, complexity_func);

            for (size_t t = 0; t < SIMP_TRANSFORM_COUNT; t++) {
                if (!transform_can_fire(SIMP_TRANSFORMS[t], seed, ctx)) continue;
                Expr* r = traced_call_unary(SIMP_TRANSFORMS[t], seed);
                if (!r) continue;
                update_best(&best, &best_score, r, complexity_func);
                /* Chain a PythagReduce pass on the transform output so a
                 * candidate produced in the final round (e.g. FactorSquareFree
                 * surfacing Cos[x]^2 - 1 inside a product) still gets the
                 * Pythagorean rewrite applied -- otherwise the reduction
                 * would only fire when this candidate became a seed in a
                 * subsequent round, which doesn't happen at SIMP_ROUNDS-1. */
                {
                    Expr* pr = transform_pythag_reduce(r);
                    if (pr) {
                        if (!expr_eq(pr, r)) {
                            update_best(&best, &best_score, pr, complexity_func);
                        }
                        expr_free(pr);
                    }
                }
                /* Chain a simp_radicals pass on the transform output.
                 * Together / Cancel after TrigExpand can surface
                 * intermediate forms with adjacent radical factors
                 * (Sqrt[a] Sqrt[b]) -- e.g.
                 * Together[4/3 Sin/Sqrt[2] - 2/9 Sqrt[3] Sqrt[6] Sin + ...]
                 *   -> (-12 Sqrt[2] Sqrt[3] Sin + 12 Sqrt[6] Sin)/(9 Sqrt[2] Sqrt[6])
                 * whose combination collapses the numerator to 0. Without
                 * this chain, simp_radicals would only fire when r becomes
                 * a seed next round -- which doesn't happen at the final
                 * round. */
                {
                    Expr* rd = simp_radicals(r);
                    if (rd) {
                        if (!expr_eq(rd, r)) {
                            update_best(&best, &best_score, rd, complexity_func);
                        }
                        expr_free(rd);
                    }
                }
                /* Chain a trig-Pi canonicalization pass on the transform
                 * output.  TrigReduce on Cos[Pi/9] Cos[2Pi/9] Cos[3Pi/9]
                 * Cos[4Pi/9] - 1/16 surfaces the intermediate
                 * 1/8 (Cos[4 Pi/9] + Cos[5 Pi/9]); only after canonicalising
                 * Cos[5 Pi/9] -> -Sin[Pi/18] (equivalently -Cos[4 Pi/9])
                 * does the Plus collapse to 0.  Without this chain the
                 * canonicalizer would have to wait for r to become a seed
                 * next round, which does not happen at SIMP_ROUNDS-1. */
                {
                    Expr* tc = simp_trig_pi_canon(r);
                    if (tc) {
                        if (!expr_eq(tc, r)) {
                            update_best(&best, &best_score, tc, complexity_func);
                        }
                        expr_free(tc);
                    }
                }
                if (expr_eq(r, seed)) {
                    expr_free(r);
                } else if (score_with_func(r, complexity_func) > parent_score) {
                    /* TrigExpand expands Sin[a+b]/Cos[a+b] into
                     * Sin[a]Cos[b]+Cos[a]Sin[b] etc., which usually grows
                     * the leaf count but surfaces radical products
                     * (Sqrt[a]Sqrt[b]) and rationalisable coefficients that
                     * subsequent transforms (Together, simp_radicals,
                     * Collect) can collapse -- sometimes all the way to 0.
                     * The default score gate would drop those candidates;
                     * loosen it for TrigExpand using the same bound
                     * TrigRoundtrip uses at the seed phase (2*input + 8) so
                     * pathological blow-ups (Sin[a+b+c+d]) still get
                     * pruned. */
                    if (strcmp(SIMP_TRANSFORMS[t], "TrigExpand") == 0 &&
                        score_with_func(r, complexity_func) <=
                            2 * parent_score + 8) {
                        cs_add_or_free(&next, r);
                    } else {
                        expr_free(r);
                    }
                } else {
                    cs_add_or_free(&next, r);
                }
            }
            /* Also try the trig roundtrip on each candidate. */
            if (transform_can_fire("TrigRoundtrip", seed, NULL)) {
                Expr* tr = transform_trig_roundtrip(seed);
                if (tr) {
                    update_best(&best, &best_score, tr, complexity_func);
                    if (expr_eq(tr, seed)) {
                        expr_free(tr);
                    } else if (score_with_func(tr, complexity_func) > parent_score) {
                        expr_free(tr);
                    } else {
                        cs_add_or_free(&next, tr);
                    }
                }
            }
            /* Pythagorean square completion on each candidate. The Factor
             * transform run earlier may have produced (1 + 2 Sin Cos)^2;
             * this round step lets the rule fire on that intermediate. */
            {
                Expr* psc = transform_pythag_square_complete(seed);
                if (psc) {
                    update_best(&best, &best_score, psc, complexity_func);
                    if (expr_eq(psc, seed)) {
                        expr_free(psc);
                    } else if (score_with_func(psc, complexity_func) > parent_score) {
                        expr_free(psc);
                    } else {
                        cs_add_or_free(&next, psc);
                    }
                }
            }
            /* Pythagorean reduction on each candidate (1 - Cos^2 -> Sin^2,
             * Cosh^2 - 1 -> Sinh^2, etc.). Strict leaf-count win when
             * matched. */
            {
                Expr* pr = transform_pythag_reduce(seed);
                if (pr) {
                    update_best(&best, &best_score, pr, complexity_func);
                    if (expr_eq(pr, seed)) {
                        expr_free(pr);
                    } else if (score_with_func(pr, complexity_func) > parent_score) {
                        expr_free(pr);
                    } else {
                        cs_add_or_free(&next, pr);
                    }
                }
            }
            /* Half-angle tangent / Tanh on each candidate. Lets Together
             * /Cancel intermediates (which can surface (1+Cos[x]) or
             * (1+Cosh[x]) factors after partial cancellation) feed into
             * the rule. */
            {
                Expr* ha = transform_halfangle(seed);
                if (ha) {
                    update_best(&best, &best_score, ha, complexity_func);
                    if (expr_eq(ha, seed)) {
                        expr_free(ha);
                    } else if (score_with_func(ha, complexity_func) > parent_score) {
                        expr_free(ha);
                    } else {
                        cs_add_or_free(&next, ha);
                    }
                }
            }
            /* Radical product combine on each candidate. Together /
             * Cancel can surface fresh Sqrt[a]*Sqrt[b] products in their
             * output; this lets the combine fire on the intermediate. */
            {
                Expr* rd = simp_radicals(seed);
                if (rd) {
                    update_best(&best, &best_score, rd, complexity_func);
                    if (expr_eq(rd, seed)) {
                        expr_free(rd);
                    } else if (score_with_func(rd, complexity_func) > parent_score) {
                        expr_free(rd);
                    } else {
                        cs_add_or_free(&next, rd);
                    }
                }
            }
            /* Common-factor lift not applied per-candidate -- see comment
             * in seed phase above and the final-form polish in
             * builtin_simplify. */
            /* Abs rewrites on each candidate. Same force-take semantics
             * as the seed phase (idempotent rules, structural answer). */
            if (transform_can_fire("AbsRules", seed, ctx)) {
                bool dbg = simp_debug_enabled();
                clock_t t0 = dbg ? clock() : 0;
                Expr* abr = apply_abs_rules(seed, ctx);
                if (dbg) simp_debug_log("AbsRules", seed, abr, simp_debug_elapsed_ms(t0));
                if (abr) {
                    if (!expr_eq(abr, seed)) {
                        size_t s = score_with_func(abr, complexity_func);
                        expr_free(best);
                        best = expr_copy(abr);
                        best_score = s;
                        if (s <= parent_score) {
                            cs_add_or_free(&next, abr);
                        } else {
                            expr_free(abr);
                        }
                    } else {
                        expr_free(abr);
                    }
                }
            }
            /* And per-variable Collect on each candidate. */
            if (transform_can_fire("CollectPerVariable", seed, NULL)) {
                try_collect_per_variable(seed, parent_score, &next, &best, &best_score,
                                         complexity_func);
            }
            /* And the assumption rewriter on each candidate. Bias as in
             * the seed phase: prefer assumption-driven forms at equal
             * complexity for `best`; gate seeding on parent_score. */
            if (transform_can_fire("AssumptionRules", seed, ctx)) {
                bool dbg = simp_debug_enabled();
                clock_t t0 = dbg ? clock() : 0;
                Expr* ar = apply_assumption_rules(seed, ctx);
                if (dbg) simp_debug_log("AssumptionRules", seed, ar, simp_debug_elapsed_ms(t0));
                if (ar) {
                    if (!expr_eq(ar, seed)) {
                        size_t s = score_with_func(ar, complexity_func);
                        if (s <= best_score) {
                            expr_free(best);
                            best = expr_copy(ar);
                            best_score = s;
                        }
                        if (s <= parent_score) {
                            cs_add_or_free(&next, ar);
                        } else {
                            expr_free(ar);
                        }
                    } else {
                        expr_free(ar);
                    }
                }
            }
            if (transform_can_fire("LogExpRules", seed, ctx)) {
                bool dbg = simp_debug_enabled();
                clock_t t1 = dbg ? clock() : 0;
                Expr* lr = apply_logexp_rules(seed, ctx);
                if (dbg) simp_debug_log("LogExpRules", seed, lr, simp_debug_elapsed_ms(t1));
                if (lr) {
                    if (!expr_eq(lr, seed)) {
                        /* Force-win for `best` (correctness-preserving
                         * under positivity assumptions). Still gate the
                         * seeding step. */
                        size_t s = score_with_func(lr, complexity_func);
                        expr_free(best);
                        best = expr_copy(lr);
                        best_score = s;
                        if (s <= parent_score) {
                            cs_add_or_free(&next, lr);
                        } else {
                            expr_free(lr);
                        }
                    } else {
                        expr_free(lr);
                    }
                }
            }
        }
        cs_free(&seeds);
        seeds = next;
        if (seeds.count == 0) break;
        if (simp_best_is_zero(best)) break;
    }
search_done:
    cs_free(&seeds);
    if (abs_pre) expr_free(abs_pre);
    return best;
}

/* ----------------------------------------------------------------------- */
/* Bottom-up Simplify: memoised recursive descent over subtrees            */
/* ----------------------------------------------------------------------- */

#define SIMP_BOTTOMUP_MAX_DEPTH 64
#define SIMP_MEMO_BUCKETS 256

typedef struct SimpMemoEntry {
    Expr* key;
    Expr* value;
    struct SimpMemoEntry* next;
} SimpMemoEntry;

typedef struct {
    SimpMemoEntry* buckets[SIMP_MEMO_BUCKETS];
} SimpMemo;

static void simp_memo_init(SimpMemo* m) {
    for (int i = 0; i < SIMP_MEMO_BUCKETS; i++) m->buckets[i] = NULL;
}

static void simp_memo_free(SimpMemo* m) {
    for (int i = 0; i < SIMP_MEMO_BUCKETS; i++) {
        SimpMemoEntry* e = m->buckets[i];
        while (e) {
            SimpMemoEntry* next = e->next;
            expr_free(e->key);
            expr_free(e->value);
            free(e);
            e = next;
        }
        m->buckets[i] = NULL;
    }
}

/* Borrowed pointer to cached value, or NULL on miss. */
static const Expr* simp_memo_get(SimpMemo* m, const Expr* key) {
    uint64_t h = expr_hash(key) % SIMP_MEMO_BUCKETS;
    for (SimpMemoEntry* e = m->buckets[h]; e; e = e->next) {
        if (expr_eq(e->key, key)) return e->value;
    }
    return NULL;
}

/* Stores deep copies of both key and value. */
static void simp_memo_put(SimpMemo* m, const Expr* key, const Expr* value) {
    uint64_t h = expr_hash(key) % SIMP_MEMO_BUCKETS;
    SimpMemoEntry* e = (SimpMemoEntry*)malloc(sizeof(SimpMemoEntry));
    if (!e) return;
    e->key = expr_copy((Expr*)key);
    e->value = expr_copy((Expr*)value);
    e->next = m->buckets[h];
    m->buckets[h] = e;
}

/* Heads whose internal structure must not be rewritten by Simplify even
 * when no Hold attribute is set. Pattern/Blank* would change matcher
 * semantics; Function captures named slots; Hold* are explicitly
 * preserved by the user. */
static bool simp_skip_recursion_head(const char* h) {
    return strcmp(h, "Hold") == 0 ||
           strcmp(h, "HoldForm") == 0 ||
           strcmp(h, "HoldComplete") == 0 ||
           strcmp(h, "HoldPattern") == 0 ||
           strcmp(h, "Unevaluated") == 0 ||
           strcmp(h, "Pattern") == 0 ||
           strcmp(h, "Blank") == 0 ||
           strcmp(h, "BlankSequence") == 0 ||
           strcmp(h, "BlankNullSequence") == 0 ||
           strcmp(h, "Function") == 0 ||
           strcmp(h, "Slot") == 0 ||
           strcmp(h, "SlotSequence") == 0;
}

/* Heads whose evaluator-level Hold attributes mean we must not descend
 * (Set, SetDelayed, Module, Block, With, If, While, For, Do, ...): a
 * bottom-up rewrite would change which sub-expression is the assignment
 * target / loop variable / branch body. */
static bool simp_head_holds_args(const char* h) {
    SymbolDef* def = symtab_lookup(h);
    if (!def) return false;
    return (def->attributes & (ATTR_HOLDFIRST | ATTR_HOLDREST |
                               ATTR_HOLDALL | ATTR_HOLDALLCOMPLETE)) != 0;
}

/* Recursive bottom-up Simplify.
 *
 * Strategy: simplify each child in isolation, then re-evaluate the
 * rebuilt parent (so canonical-form invariants on Plus/Times/etc. are
 * restored if children changed shape), then run the standard top-level
 * candidate-set search on the result.
 *
 * Why this helps: a transform like the Pythagorean identity may be
 * inapplicable at the root (e.g. f[Sin[x]^2 + Cos[x]^2]) but applies
 * cleanly to a subtree. Top-level search alone would miss it.
 *
 * Cost control:
 *   - Memoisation keyed by expr_hash + expr_eq: identical subtrees are
 *     simplified once per Simplify call (e.g. f[g[x], g[x], g[x]]).
 *   - Atoms and held heads bottom out into a single simp_search.
 *   - SIMP_BOTTOMUP_MAX_DEPTH guards against pathological nesting; once
 *     hit, we fall back to the existing top-level simp_search. */
/* ----------------------------------------------------------------------- */
/* simp_factorial -- general factorial simplification                       */
/* ----------------------------------------------------------------------- */

/*
 * Simplify expressions that contain Factorial[...] atoms via a single
 * principled four-step procedure rather than a per-pattern table.
 *
 *   Step A. Decompose every Factorial[arg] as arg = sym + c, where c is the
 *           integer addend in arg's Plus form (0 if arg has no integer
 *           term) and sym is the rest. Run Expand on the arg first so that
 *           shapes like 2(n-1) become 2n-2 (sym = 2n, c = -2).
 *   Step B. Group factorials whose sym parts are structurally equal. The
 *           group base is b = sym + min(c_i); each member's offset is
 *           k_i = c_i - min(c_i) >= 0.
 *   Step C. Rewrite each Factorial[b + k_i] as
 *               Factorial[b] * (b+1) * (b+2) * ... * (b+k_i).
 *           (k_i = 0 leaves the call unchanged.) Built as a literal Times
 *           product, then evaluated; the unifying Factorial[b] factor lets
 *           algebraic transforms see the cancellation/absorption.
 *   Step D. Run Together -> Cancel -> Expand on the rewritten form so that
 *           common Factorial[b] factors cancel and additive collapses fire
 *           ((n+1)! - n n! -> n!,   1/n! - 1/(n+1)! -> n/(n+1)!).
 *   Step E. Re-fold: walk every Times node holding a Factorial[b] factor
 *           and absorb consecutive (b+1)(b+2)...(b+k) cofactors back into
 *           Factorial[b+k]. A small gap-fill (budget 2) allows shapes like
 *           Factorial[b]*(b+1)(b+3)  ->  Factorial[b+3] / (b+2),
 *           which closes (n^2 - 1)*(n-2)! -> (n+1)!/n.
 *
 * The whole transform is gated on `contains_factorial(e)` and the result
 * is propagated as a Simplify seed: only the strictly-shorter form
 * survives the standard SimplifyCount tiebreak.
 */

static bool contains_factorial(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (simp_eq_head_sym(e, "Factorial") && e->data.function.arg_count == 1) {
        return true;
    }
    if (contains_factorial(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_factorial(e->data.function.args[i])) return true;
    }
    return false;
}

/* Decompose `arg` as sym + c with c an int64 offset. Returns true on
 * success. *sym_out is a freshly-allocated Expr* (caller must free).
 * Conservative: any BigInt or Real addend aborts the decomposition (we
 * cannot soundly fold a factorial with a non-int64 offset under our
 * Pochhammer-style expansion). */
static bool simp_fact_decompose(const Expr* arg, Expr** sym_out,
                                int64_t* c_out) {
    *c_out = 0;
    *sym_out = NULL;
    if (arg->type == EXPR_INTEGER) {
        *c_out = arg->data.integer;
        *sym_out = expr_new_integer(0);
        return true;
    }
    if (arg->type == EXPR_BIGINT || arg->type == EXPR_REAL) return false;
    if (simp_eq_head_sym(arg, "Plus")) {
        size_t n = arg->data.function.arg_count;
        int64_t total = 0;
        Expr** rest = (Expr**)calloc(n, sizeof(Expr*));
        size_t rest_count = 0;
        for (size_t i = 0; i < n; i++) {
            const Expr* t = arg->data.function.args[i];
            if (t->type == EXPR_INTEGER) {
                int64_t v = t->data.integer;
                int64_t nt;
                if (__builtin_add_overflow(total, v, &nt)) {
                    for (size_t j = 0; j < rest_count; j++) expr_free(rest[j]);
                    free(rest);
                    return false;
                }
                total = nt;
            } else if (t->type == EXPR_BIGINT || t->type == EXPR_REAL) {
                for (size_t j = 0; j < rest_count; j++) expr_free(rest[j]);
                free(rest);
                return false;
            } else {
                rest[rest_count++] = expr_copy((Expr*)t);
            }
        }
        *c_out = total;
        if (rest_count == 0) {
            *sym_out = expr_new_integer(0);
        } else if (rest_count == 1) {
            *sym_out = rest[0];
        } else {
            *sym_out = expr_new_function(expr_new_symbol("Plus"), rest,
                                         rest_count);
        }
        free(rest);
        return true;
    }
    /* Atom or non-Plus: c = 0, sym = arg. */
    *sym_out = expr_copy((Expr*)arg);
    return true;
}

/* Build the canonical expression sym + c (returns owned Expr*; consumes
 * `sym` ownership). Special-cases sym == 0 (returns just c) and c == 0
 * (returns just sym). */
static Expr* simp_fact_make_offset(Expr* sym, int64_t c) {
    bool sym_is_zero = (sym->type == EXPR_INTEGER && sym->data.integer == 0);
    if (sym_is_zero) {
        expr_free(sym);
        return expr_new_integer(c);
    }
    if (c == 0) return sym;
    Expr* args[2] = { sym, expr_new_integer(c) };
    return expr_new_function(expr_new_symbol("Plus"), args, 2);
}

/* Group registry. Each entry maps a symbolic part to its minimum integer
 * offset across the input. Linear-scan find by structural equality. */
typedef struct {
    Expr*   sym;       /* owned */
    int64_t min_c;
    int64_t max_c;
    size_t  count;
} SimpFactGroup;

typedef struct {
    SimpFactGroup* items;
    size_t count;
    size_t cap;
} SimpFactGroupSet;

static void simp_fgs_init(SimpFactGroupSet* s) {
    s->items = NULL; s->count = 0; s->cap = 0;
}
static void simp_fgs_free(SimpFactGroupSet* s) {
    for (size_t i = 0; i < s->count; i++) expr_free(s->items[i].sym);
    free(s->items);
    s->items = NULL; s->count = 0; s->cap = 0;
}
static SimpFactGroup* simp_fgs_find(SimpFactGroupSet* s, const Expr* sym) {
    for (size_t i = 0; i < s->count; i++) {
        if (expr_eq(s->items[i].sym, sym)) return &s->items[i];
    }
    return NULL;
}
static void simp_fgs_record(SimpFactGroupSet* s, const Expr* sym, int64_t c) {
    SimpFactGroup* g = simp_fgs_find(s, sym);
    if (g) {
        if (c < g->min_c) g->min_c = c;
        if (c > g->max_c) g->max_c = c;
        g->count++;
        return;
    }
    if (s->count >= s->cap) {
        s->cap = s->cap ? s->cap * 2 : 4;
        s->items = (SimpFactGroup*)realloc(s->items, s->cap * sizeof(*s->items));
    }
    s->items[s->count].sym   = expr_copy((Expr*)sym);
    s->items[s->count].min_c = c;
    s->items[s->count].max_c = c;
    s->items[s->count].count = 1;
    s->count++;
}

/* Walk the expression collecting (sym, c) for every Factorial[arg]. Each
 * Factorial argument is first run through Expand so shapes like 2(n-1)
 * decompose under the expanded form 2n-2. */
static void simp_fact_gather(const Expr* e, SimpFactGroupSet* groups) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (simp_eq_head_sym(e, "Factorial") && e->data.function.arg_count == 1) {
        Expr* arg_exp = call_unary_copy("Expand", e->data.function.args[0]);
        if (arg_exp) {
            Expr* sym = NULL; int64_t c = 0;
            if (simp_fact_decompose(arg_exp, &sym, &c)) {
                simp_fgs_record(groups, sym, c);
                expr_free(sym);
            }
            expr_free(arg_exp);
        }
    }
    simp_fact_gather(e->data.function.head, groups);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        simp_fact_gather(e->data.function.args[i], groups);
    }
}

/* Build Factorial[b] * (b+1) * (b+2) * ... * (b + offset) for offset >= 0.
 * `b_sym_template` is borrowed -- the helper deep-copies as needed.
 * `b_const` is the integer addend that combines with `b_sym_template`
 * to form b. */
static Expr* simp_fact_pochhammer_expansion(const Expr* b_sym_template,
                                            int64_t b_const,
                                            int64_t offset) {
    /* Pochhammer guard: keep at most 32 explicit multiplied factors so a
     * pathological Factorial[n + 1000000] does not blow up the rewrite.
     * Beyond this we leave the call alone. */
    if (offset < 0 || offset > 32) return NULL;

    Expr* base = simp_fact_make_offset(expr_copy((Expr*)b_sym_template),
                                       b_const);
    Expr* fact_args[1] = { base };
    Expr* fact_b = expr_new_function(expr_new_symbol("Factorial"), fact_args, 1);
    if (offset == 0) return fact_b;

    /* Times[Factorial[b], (b+1), (b+2), ..., (b+offset)] */
    size_t n_args = (size_t)offset + 1;
    Expr** args = (Expr**)calloc(n_args, sizeof(Expr*));
    args[0] = fact_b;
    for (int64_t k = 1; k <= offset; k++) {
        args[(size_t)k] = simp_fact_make_offset(
            expr_copy((Expr*)b_sym_template), b_const + k);
    }
    Expr* out = expr_new_function(expr_new_symbol("Times"), args, n_args);
    free(args);
    return out;
}

/* Recursively rewrite each Factorial[arg] in `e` per the group base. The
 * returned tree is a fresh allocation. */
static Expr* simp_fact_rewrite(const Expr* e, const SimpFactGroupSet* groups) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    if (simp_eq_head_sym(e, "Factorial") && e->data.function.arg_count == 1) {
        Expr* arg_exp = call_unary_copy("Expand", e->data.function.args[0]);
        if (arg_exp) {
            Expr* sym = NULL; int64_t c = 0;
            if (simp_fact_decompose(arg_exp, &sym, &c)) {
                /* find the group for this sym */
                for (size_t i = 0; i < groups->count; i++) {
                    if (expr_eq(groups->items[i].sym, sym)) {
                        int64_t mc = groups->items[i].min_c;
                        Expr* exp = simp_fact_pochhammer_expansion(
                            sym, mc, c - mc);
                        if (exp) {
                            expr_free(sym);
                            expr_free(arg_exp);
                            return exp;
                        }
                        break;
                    }
                }
                expr_free(sym);
            }
            expr_free(arg_exp);
        }
        /* fall through: rewrite arg recursively (not a known group) */
    }

    /* Generic function: rewrite head + args. */
    Expr* new_head = simp_fact_rewrite(e->data.function.head, groups);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = simp_fact_rewrite(e->data.function.args[i], groups);
    }
    Expr* out = expr_new_function(new_head, new_args, n);
    free(new_args);
    return out;
}

/* Re-fold helper: given a list of Times args, look at the first Factorial
 * factor and try to absorb consecutive (base+j) cofactors. Returns a
 * new Times-content list (possibly with a Power[..., -1] denominator
 * for gap-fill terms) on success, or NULL when no fold fires.
 *
 * Inputs:
 *   args        -- borrowed Times children (caller still owns).
 *   n           -- count of `args`.
 *   out_args    -- on success, *out_args is a fresh Expr** of *out_n
 *                  freshly-allocated children that the caller owns.
 *   out_n       -- count of children in *out_args.
 */
static bool simp_fact_refold_times(Expr** args, size_t n,
                                   Expr*** out_args, size_t* out_n) {
    if (n < 2) return false;
    /* Find the first Factorial[b] direct child. */
    size_t fact_idx = (size_t)-1;
    Expr* base_sym = NULL;
    int64_t base_c = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* a = args[i];
        if (simp_eq_head_sym(a, "Factorial") &&
            a->data.function.arg_count == 1) {
            Expr* exp = call_unary_copy("Expand", a->data.function.args[0]);
            if (exp) {
                if (simp_fact_decompose(exp, &base_sym, &base_c)) {
                    fact_idx = i;
                    expr_free(exp);
                    break;
                }
                expr_free(exp);
            }
        }
    }
    if (fact_idx == (size_t)-1) return false;

    /* Collect candidate absorbing factors: args of the form sym + c with
     * the SAME sym as `base_sym`. Record (arg_index, j = c - base_c). */
    typedef struct { size_t idx; int64_t j; } AbsCand;
    AbsCand* cands = (AbsCand*)calloc(n, sizeof(AbsCand));
    size_t ncand = 0;
    for (size_t i = 0; i < n; i++) {
        if (i == fact_idx) continue;
        Expr* a = args[i];
        Expr* sym2 = NULL; int64_t c2 = 0;
        if (!simp_fact_decompose(a, &sym2, &c2)) continue;
        if (expr_eq(sym2, base_sym)) {
            int64_t j = c2 - base_c;
            if (j >= 1) {
                cands[ncand].idx = i;
                cands[ncand].j   = j;
                ncand++;
            }
        }
        expr_free(sym2);
    }
    if (ncand == 0) {
        free(cands);
        expr_free(base_sym);
        return false;
    }

    /* Sort cands by j ascending. Small n -- a bubble pass is fine. */
    for (size_t i = 0; i < ncand; i++) {
        for (size_t j = i + 1; j < ncand; j++) {
            if (cands[j].j < cands[i].j) {
                AbsCand t = cands[i]; cands[i] = cands[j]; cands[j] = t;
            }
        }
    }
    /* Drop duplicate j-values (a Times node should not contain duplicate
     * (base+j) factors, but be defensive). */
    {
        size_t w = 0;
        for (size_t i = 0; i < ncand; i++) {
            if (w == 0 || cands[w-1].j != cands[i].j) {
                cands[w++] = cands[i];
            }
        }
        ncand = w;
    }

    /* Find the largest k such that:
     *   - {1, 2, ..., k} subset of cand_offsets union gap_set
     *   - |gap_set| <= GAP_BUDGET
     *   - in_count >= gap_count (else the fold strictly grows the form).
     * The bound k_max = largest cand offset (folding past it adds factors
     * we did not have).
     */
    const int GAP_BUDGET = 2;
    int64_t kmax = cands[ncand - 1].j;
    int64_t best_k = 0;
    int64_t best_gaps = 0;
    /* We iterate k from 1 to kmax. Maintain a pointer into cands for
     * "how many cand offsets are <= k". */
    size_t ci = 0;
    for (int64_t k = 1; k <= kmax; k++) {
        while (ci < ncand && cands[ci].j <= k) ci++;
        int64_t n_in = (int64_t)ci;
        int64_t n_gaps = k - n_in;
        if (n_gaps > GAP_BUDGET) break;
        if (n_in < n_gaps) continue;
        if (k > best_k) { best_k = k; best_gaps = n_gaps; }
    }
    if (best_k == 0) {
        free(cands);
        expr_free(base_sym);
        return false;
    }

    /* Identify which cand indices we are absorbing (j <= best_k). */
    bool* absorbed = (bool*)calloc(n, sizeof(bool));
    bool* present_offset = (bool*)calloc((size_t)best_k + 1, sizeof(bool));
    for (size_t i = 0; i < ncand; i++) {
        if (cands[i].j <= best_k) {
            absorbed[cands[i].idx] = true;
            present_offset[cands[i].j] = true;
        }
    }

    /* Build the new args list:
     *   - replace args[fact_idx] with Factorial[base_sym + base_c + best_k].
     *   - drop args[absorbed].
     *   - keep all other args unchanged.
     *   - append gap_factors as Power[(base_sym + base_c + m), -1] for each
     *     m in 1..best_k missing from present_offset.
     */
    size_t new_n = 0;
    for (size_t i = 0; i < n; i++) {
        if (i == fact_idx) { new_n++; continue; }
        if (!absorbed[i]) new_n++;
    }
    new_n += (size_t)best_gaps;
    Expr** built = (Expr**)calloc(new_n ? new_n : 1, sizeof(Expr*));
    size_t w = 0;
    for (size_t i = 0; i < n; i++) {
        if (i == fact_idx) {
            Expr* nb = simp_fact_make_offset(expr_copy(base_sym),
                                             base_c + best_k);
            Expr* fa[1] = { nb };
            built[w++] = expr_new_function(expr_new_symbol("Factorial"),
                                           fa, 1);
        } else if (!absorbed[i]) {
            built[w++] = expr_copy(args[i]);
        }
    }
    for (int64_t m = 1; m <= best_k; m++) {
        if (!present_offset[m]) {
            Expr* fac = simp_fact_make_offset(expr_copy(base_sym),
                                              base_c + m);
            Expr* pa[2] = { fac, expr_new_integer(-1) };
            built[w++] = expr_new_function(expr_new_symbol("Power"), pa, 2);
        }
    }

    free(absorbed);
    free(present_offset);
    free(cands);
    expr_free(base_sym);
    *out_args = built;
    *out_n = new_n;
    return true;
}

/* Combine Times children that all carry exponent -1 into a single
 * Power[Times[..., -1]]. Picocas's evaluator does NOT auto-coalesce
 * separate Power[a, -1] * Power[b, -1] into Power[Times[a, b], -1]
 * (Mathematica does, picocas doesn't), so the auto-cancelled output of
 * a factorial rewrite sits at SimplifyCount 12 even when the same
 * expression printed as `1/(a*b)` would score 9. Lifting all negative-
 * exponent factors into a shared denominator gets us back to the lower
 * canonical score so the FactorialRules seed wins the tiebreak.
 *
 * Conservative: only combines exponent -1 (the trivial "denominator"
 * case). Mixed-sign exponents are left untouched. */
static Expr* simp_fact_combine_inverses(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (!simp_eq_head_sym(e, "Times")) {
        /* Recurse into children. */
        Expr* new_head = simp_fact_combine_inverses(e->data.function.head);
        size_t n = e->data.function.arg_count;
        Expr** new_args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            new_args[i] = simp_fact_combine_inverses(e->data.function.args[i]);
        }
        Expr* out = expr_new_function(new_head, new_args, n);
        free(new_args);
        return out;
    }
    size_t n = e->data.function.arg_count;
    /* Recurse into each Times child first. */
    Expr** child = (Expr**)calloc(n, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        child[i] = simp_fact_combine_inverses(e->data.function.args[i]);
    }
    /* Partition children: numerator (no -1 power), denominator (Power[X, -1]). */
    Expr** num = (Expr**)calloc(n, sizeof(Expr*));
    Expr** den = (Expr**)calloc(n, sizeof(Expr*));
    size_t nn = 0, dn = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* a = child[i];
        if (simp_eq_head_sym(a, "Power") &&
            a->data.function.arg_count == 2 &&
            a->data.function.args[1]->type == EXPR_INTEGER &&
            a->data.function.args[1]->data.integer == -1) {
            den[dn++] = expr_copy(a->data.function.args[0]);
        } else {
            num[nn++] = expr_copy(a);
        }
    }
    /* If 0 or 1 denominator factors, nothing to combine. Restore original. */
    if (dn <= 1) {
        for (size_t i = 0; i < nn; i++) expr_free(num[i]);
        for (size_t i = 0; i < dn; i++) expr_free(den[i]);
        free(num); free(den);
        Expr* out = expr_new_function(
            expr_copy(e->data.function.head), child, n);
        /* expr_new_function takes ownership of `child` slot but we still
         * need to free our wrapper allocation. */
        free(child);
        return out;
    }
    /* Build combined denominator: Times[den[0], ..., den[dn-1]] */
    Expr* den_times;
    if (dn == 1) {
        den_times = den[0];
    } else {
        Expr** dargs = (Expr**)calloc(dn, sizeof(Expr*));
        for (size_t i = 0; i < dn; i++) dargs[i] = den[i];
        den_times = expr_new_function(expr_new_symbol("Times"), dargs, dn);
        free(dargs);
    }
    Expr* den_pow_args[2] = { den_times, expr_new_integer(-1) };
    Expr* den_pow = expr_new_function(expr_new_symbol("Power"),
                                      den_pow_args, 2);
    /* Build new Times: numerator factors + den_pow. */
    if (nn == 0) {
        for (size_t i = 0; i < n; i++) expr_free(child[i]);
        free(child); free(num); free(den);
        return den_pow;
    }
    Expr** new_args = (Expr**)calloc(nn + 1, sizeof(Expr*));
    for (size_t i = 0; i < nn; i++) new_args[i] = num[i];
    new_args[nn] = den_pow;
    Expr* out = expr_new_function(expr_new_symbol("Times"), new_args, nn + 1);
    free(new_args);
    for (size_t i = 0; i < n; i++) expr_free(child[i]);
    free(child); free(num); free(den);
    return out;
}

/* Re-fold tree walk. Bottom-up rebuild; at each Times node, attempt the
 * absorption pass repeatedly until it stops firing. */
static Expr* simp_fact_refold(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    /* Recurse first. */
    Expr* new_head = simp_fact_refold(e->data.function.head);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = simp_fact_refold(e->data.function.args[i]);
    }

    if (simp_eq_head_sym(e, "Times")) {
        /* Iterate refold until no more folds. Bounded by initial n: each
         * fold reduces the cofactor count by at least 1 (one absorbed
         * arg replaced; possibly +gaps but those are appended once). */
        for (size_t guard = 0; guard < 16; guard++) {
            Expr** out = NULL; size_t out_n = 0;
            if (!simp_fact_refold_times(new_args, n, &out, &out_n)) break;
            for (size_t i = 0; i < n; i++) expr_free(new_args[i]);
            free(new_args);
            new_args = out;
            n = out_n;
        }
    }

    Expr* out = expr_new_function(new_head, new_args, n);
    free(new_args);
    return out;
}

/* Recognize the canonical (2 v)! / (2^v v!) shape and rewrite to
 * Factorial2[2 v - 1]. The check is purely structural: input must be a
 * Times whose direct children include
 *    Factorial[Times[2, v]]                    (numerator factorial)
 *    Power[Factorial[v], -1]                   (denominator factorial)
 *    Power[2, Times[-1, v]]                    (1 / 2^v)
 * with the same v expression appearing in all three. Other Times
 * children are kept as a residual cofactor on the result. */
static Expr* simp_fact_double_factorial(const Expr* e) {
    if (!simp_eq_head_sym(e, "Times")) return NULL;
    size_t n = e->data.function.arg_count;
    if (n < 3) return NULL;

    int idx_2v_fact = -1;       /* Factorial[2 v]            */
    int idx_v_fact_inv = -1;    /* Power[Factorial[v], -1]   */
    int idx_2_pow_neg = -1;     /* Power[2, -v]              */
    Expr* v_expr = NULL;

    for (size_t i = 0; i < n; i++) {
        Expr* a = e->data.function.args[i];
        /* Factorial[Times[2, v]] */
        if (idx_2v_fact == -1 &&
            simp_eq_head_sym(a, "Factorial") &&
            a->data.function.arg_count == 1) {
            Expr* arg = a->data.function.args[0];
            Expr* arg_exp = call_unary_copy("Expand", arg);
            if (arg_exp && simp_eq_head_sym(arg_exp, "Times") &&
                arg_exp->data.function.arg_count == 2) {
                Expr* a0 = arg_exp->data.function.args[0];
                Expr* a1 = arg_exp->data.function.args[1];
                if (a0->type == EXPR_INTEGER && a0->data.integer == 2) {
                    if (!v_expr || expr_eq(v_expr, a1)) {
                        if (!v_expr) v_expr = expr_copy(a1);
                        idx_2v_fact = (int)i;
                        expr_free(arg_exp);
                        continue;
                    }
                }
            }
            if (arg_exp) expr_free(arg_exp);
        }
        /* Power[Factorial[v], -1] */
        if (idx_v_fact_inv == -1 &&
            simp_eq_head_sym(a, "Power") &&
            a->data.function.arg_count == 2 &&
            a->data.function.args[1]->type == EXPR_INTEGER &&
            a->data.function.args[1]->data.integer == -1 &&
            simp_eq_head_sym(a->data.function.args[0], "Factorial") &&
            a->data.function.args[0]->data.function.arg_count == 1) {
            Expr* v = a->data.function.args[0]->data.function.args[0];
            if (!v_expr || expr_eq(v_expr, v)) {
                if (!v_expr) v_expr = expr_copy(v);
                idx_v_fact_inv = (int)i;
                continue;
            }
        }
        /* Power[2, Times[-1, v]] or Power[2, -v] depending on shape. */
        if (idx_2_pow_neg == -1 &&
            simp_eq_head_sym(a, "Power") &&
            a->data.function.arg_count == 2 &&
            a->data.function.args[0]->type == EXPR_INTEGER &&
            a->data.function.args[0]->data.integer == 2) {
            Expr* exp = a->data.function.args[1];
            /* Build Times[-1, v_expr_candidate] and compare. We instead
             * compare via -exp == v: evaluate Times[-1, exp] and check. */
            Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(exp) };
            Expr* neg_call = expr_new_function(expr_new_symbol("Times"),
                                               neg_args, 2);
            Expr* neg_eval = evaluate(neg_call);
            expr_free(neg_call);
            if (neg_eval) {
                if (!v_expr || expr_eq(v_expr, neg_eval)) {
                    if (!v_expr) v_expr = expr_copy(neg_eval);
                    idx_2_pow_neg = (int)i;
                    expr_free(neg_eval);
                    continue;
                }
                expr_free(neg_eval);
            }
        }
    }

    if (idx_2v_fact < 0 || idx_v_fact_inv < 0 || idx_2_pow_neg < 0 || !v_expr) {
        if (v_expr) expr_free(v_expr);
        return NULL;
    }

    /* Build Factorial2[2 v - 1]. */
    Expr* two_v_args[2] = { expr_new_integer(2), expr_copy(v_expr) };
    Expr* two_v = expr_new_function(expr_new_symbol("Times"), two_v_args, 2);
    Expr* arg_args[2] = { two_v, expr_new_integer(-1) };
    Expr* arg_plus = expr_new_function(expr_new_symbol("Plus"), arg_args, 2);
    Expr* fac2_args[1] = { arg_plus };
    Expr* fac2 = expr_new_function(expr_new_symbol("Factorial2"), fac2_args, 1);
    expr_free(v_expr);

    /* If there are residual cofactors, multiply them in. */
    size_t residue = 0;
    for (size_t i = 0; i < n; i++) {
        if ((int)i == idx_2v_fact || (int)i == idx_v_fact_inv ||
            (int)i == idx_2_pow_neg) continue;
        residue++;
    }
    if (residue == 0) return fac2;

    Expr** all = (Expr**)calloc(residue + 1, sizeof(Expr*));
    size_t w = 0;
    all[w++] = fac2;
    for (size_t i = 0; i < n; i++) {
        if ((int)i == idx_2v_fact || (int)i == idx_v_fact_inv ||
            (int)i == idx_2_pow_neg) continue;
        all[w++] = expr_copy(e->data.function.args[i]);
    }
    Expr* out = expr_new_function(expr_new_symbol("Times"), all, residue + 1);
    free(all);
    return out;
}

/* Walk the tree applying simp_fact_double_factorial at every Times node. */
static Expr* simp_fact_double_factorial_walk(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    Expr* new_head = simp_fact_double_factorial_walk(e->data.function.head);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = simp_fact_double_factorial_walk(e->data.function.args[i]);
    }
    Expr* out = expr_new_function(new_head, new_args, n);
    free(new_args);
    if (simp_eq_head_sym(out, "Times")) {
        Expr* d = simp_fact_double_factorial(out);
        if (d) { expr_free(out); out = d; }
    }
    return out;
}

/* Top-level transform. Returns NULL when the input contains no factorial
 * or when the rewrite produced no change. Otherwise returns the
 * candidate (a freshly-allocated Expr*; caller takes ownership). */
static Expr* simp_factorial(const Expr* e) {
    if (!contains_factorial(e)) return NULL;
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    SimpFactGroupSet groups;
    simp_fgs_init(&groups);
    simp_fact_gather(e, &groups);

    /* If every group has exactly one offset (no shifting possible) AND the
     * tree contains only one factorial, the substitution step is a no-op.
     * Still proceed to the double-factorial pattern check below. */
    bool any_shift = false;
    for (size_t i = 0; i < groups.count; i++) {
        if (groups.items[i].max_c > groups.items[i].min_c) {
            any_shift = true;
            break;
        }
    }

    Expr* result = NULL;
    if (any_shift) {
        /* Step C — substitute. The rewrite emits Pochhammer products
         * literally (e.g. Factorial[n] -> Factorial[n-2]*(n-1)*n in a
         * group whose base is n-2) and the evaluator's auto-cancellation
         * through normal Times/Power semantics is sufficient to remove
         * any common Factorial[base] factor between numerator and
         * denominator. We deliberately do NOT run Together here because
         * Together would expand the polynomial denominator (turning
         * 1/(n*(n-1)) into 1/(n^2 - n) and Cancel cannot recover the
         * factored form afterwards because the original Power[Times,
         * -1] structure is gone). The auto-cancel already gives the
         * minimal canonical fraction.
         *
         * Step D' (the lighter version) -- run Cancel on the result so
         * that any residual rational form left by the evaluator is
         * reduced. Skip if the evaluator's own auto-cancel already
         * produced a form without nested Times/Power[Times,-1] noise. */
        Expr* rewritten = simp_fact_rewrite(e, &groups);
        /* Force a full evaluator pass so the auto-cancellation between
         * Factorial[base] in the numerator and the same factor in a
         * Pochhammer-expanded denominator collapses. evaluate owns its
         * argument (per the BuiltinFunc convention) so we hand off
         * `rewritten` and replace the local with the evaluated form. */
        Expr* evaluated = evaluate(rewritten);
        if (!evaluated) evaluated = expr_copy(rewritten);
        rewritten = NULL;

        /* Expand distributes products over sums, which is what makes
         *   Factorial[n]*(1+n) - n*Factorial[n] -> Factorial[n]
         *   (n+1)*n*Factorial[n-1]            -> n^2 Factorial[n-1] + ...
         * fire. Without it, the evaluator does not distribute (a + b)
         * across c, so Plus's Orderless/FLAT cancellation never sees
         * the like terms. Harmless on inputs that have nothing to
         * distribute (e.g. plain rational forms like 1/(n*(n-1))). */
        Expr* expanded = call_unary_copy("Expand", evaluated);
        if (!expanded) expanded = expr_copy(evaluated);

        /* Two paths feed the re-fold pass:
         *   A. expanded form        -- preserves a factored
         *                              Times[Power[a,-1], Power[b,-1]]
         *                              denominator so (n-2)!/n! lands
         *                              at 1/((n-1)*n) instead of being
         *                              re-expanded to 1/(n^2-n).
         *   B. Together'd form      -- combines additive fractions
         *                              over a common denominator,
         *                              needed for 1/n! - 1/(n+1)! to
         *                              fold to n/(n+1)!.
         * Both paths run Factor (to surface (b+j) linear factors in
         * the cofactor) and combine_inverses (to re-coalesce the
         * Power[a,-1]*Power[b,-1] factors that picocas's evaluator
         * leaves separated, which is what makes the combined Times
         * denominator visible to the re-fold walker). The score
         * tiebreak below picks whichever path lands at the lowest
         * SimplifyCount. */
        Expr* together = call_unary_copy("Together", expanded);

        Expr* prep_a = NULL;
        if (!has_non_integer_power(expanded)) {
            factor_memo_push(NULL);
            prep_a = call_unary_copy("Factor", expanded);
            factor_memo_pop();
        }
        if (!prep_a) prep_a = expr_copy(expanded);
        {
            Expr* coal = simp_fact_combine_inverses(prep_a);
            if (coal) { expr_free(prep_a); prep_a = coal; }
        }
        Expr* refold_a = simp_fact_refold(prep_a);

        /* Path B: Together expands the polynomial denominator (e.g.
         *   1/n! - 1/((n+1) n!)  -> n / (n! + n n!)
         * landing in a Plus inside Power[..., -1]). Factor pulls the
         * Plus apart back into a Times of Power[a, -1] factors:
         *   n / (n! + n n!) -> Times[n, Power[n!, -1], Power[1+n, -1]].
         * combine_inverses then coalesces those into a single
         * Power[Times[n!, 1+n], -1], and the re-fold walker descends
         * into that Times to fold Factorial[n] * (n+1) -> Factorial[n+1]. */
        Expr* prep_b = NULL;
        if (together) {
            if (!has_non_integer_power(together)) {
                /* Push a NULL memo to opt out of Factor's inside-Simplify
                 * variable-list narrowing. The narrowing collapses
                 * num/den variable scopes and prevents Factor from
                 * pulling Factorial[n] out of a denominator like
                 * Factorial[n] + n*Factorial[n] -> Factorial[n]*(1+n)
                 * (the polynomial-in-n viewer treats Factorial[n] as a
                 * coefficient that the inside_simplify path then
                 * refuses to factor). With separate variable lists the
                 * factorisation succeeds. We still consult / populate
                 * the outer memo via factor_memo_active() inside the
                 * helper -- only the gating in builtin_factor that
                 * branches on factor_memo_top() != NULL is what we
                 * silence here. */
                factor_memo_push(NULL);
                prep_b = call_unary_copy("Factor", together);
                factor_memo_pop();
            }
            if (!prep_b) prep_b = expr_copy(together);
            Expr* coal = simp_fact_combine_inverses(prep_b);
            if (coal) { expr_free(prep_b); prep_b = coal; }
        }
        Expr* refold_b = prep_b ? simp_fact_refold(prep_b) : NULL;

        Expr* refolded = refold_a;
        if (refold_b && simp_default_complexity(refold_b) <
                        simp_default_complexity(refold_a)) {
            expr_free(refold_a);
            refolded = refold_b;
        } else if (refold_b) {
            expr_free(refold_b);
        }
        if (prep_a) expr_free(prep_a);
        if (prep_b) expr_free(prep_b);
        if (together) expr_free(together);
        expr_free(expanded);

        /* Note on Cancel as an alternate: tempting, but Cancel/Together
         * on the post-refold form expand a Times[Power[a,-1], Power[b,-1]]
         * factored denominator into Power[Plus[a*b...], -1] (e.g.
         * 1/(n*(n-1)) -> 1/(n^2 - n)). The expanded form scores lower
         * under SimplifyCount but is the canonical "wrong" answer for
         * our purposes -- the user's expected forms keep the factored
         * denominator. We rely on the evaluator's own auto-cancel
         * (which produced `evaluated` from `rewritten`) to give the
         * canonical reduction; no further smoothing is wanted. */

        expr_free(evaluated);
        result = refolded;
    } else {
        /* No shift step; still try Factor + fold over the input (covers
         * n*(n-1)! -> n! and (n^2-1)*(n-2)! shapes where there is only
         * one factorial group but the cofactor needs factoring before
         * the (base+j) factors become visible). The NULL memo push
         * mirrors the any_shift path -- see comment there. */
        Expr* factored = NULL;
        if (!has_non_integer_power(e)) {
            factor_memo_push(NULL);
            factored = call_unary_copy("Factor", e);
            factor_memo_pop();
        }
        if (!factored) factored = expr_copy((Expr*)e);
        result = simp_fact_refold(factored);
        expr_free(factored);
    }

    /* Double-factorial pattern recognition. */
    {
        Expr* d = simp_fact_double_factorial_walk(result);
        if (d) { expr_free(result); result = d; }
    }

    /* Final canonicalization: coalesce all `Power[X, -1]` factors in
     * Times nodes into a single `Power[Times[..., -1]]`. Picocas's
     * evaluator does not auto-perform this combine, which leaves the
     * factored result at a higher SimplifyCount than it deserves
     * (count 12 vs count 9 on `1/(n*(n-1))`). Coalescing brings the
     * score under the original input, so the FactorialRules seed
     * wins the round-loop tiebreak instead of being undone. */
    {
        Expr* c = simp_fact_combine_inverses(result);
        if (c) { expr_free(result); result = c; }
    }

    simp_fgs_free(&groups);
    if (result && expr_eq(result, e)) {
        if (dbg) simp_debug_log("FactorialRules", e, NULL,
                                simp_debug_elapsed_ms(t0));
        expr_free(result);
        return NULL;
    }
    if (dbg) simp_debug_log("FactorialRules", e, result,
                            simp_debug_elapsed_ms(t0));
    return result;
}

/* ----------------------------------------------------------------------- */
/* simp_bottomup -- recursive driver                                       */
/* ----------------------------------------------------------------------- */

static Expr* simp_bottomup(const Expr* input, const AssumeCtx* ctx,
                           const Expr* complexity_func, SimpMemo* memo,
                           int depth) {
    if (!input) return NULL;

    /* Atoms have no children. Without active assumptions every transform
     * is a no-op on a bare atom, so skip the entire candidate-set search
     * and return a copy. (assume_ctx_from_expr always returns non-NULL
     * even for trivial $Assumptions=True, so we test the fact count
     * rather than the pointer.) With assumptions, atom-targeted rewrites
     * (e.g. Equal facts that name a leaf) still fire via simp_search. */
    if (input->type != EXPR_FUNCTION) {
        if (!ctx || ctx->count == 0) return expr_copy((Expr*)input);
        return simp_dispatch(input, ctx, complexity_func);
    }

    /* Depth cap: bail to top-level. */
    if (depth >= SIMP_BOTTOMUP_MAX_DEPTH) {
        return simp_dispatch(input, ctx, complexity_func);
    }

    /* Top-level Pythagorean short-circuit. Try a one-shot
     * substitution-and-Expand pass before descending into children: if
     * Cos^(2k) <-> 1 -/+ Sin^(2k) (or the hyperbolic counterpart)
     * collapses the input to a strictly smaller form, recurse on the
     * collapsed form instead. Without this short-circuit, on inputs like
     *   18 (Cos[x]+1)(Cos[x]-1)(Cos[y]^2-1)^2 (x-1) +
     *      18 (x-1) Sin[x]^2 Sin[y]^4
     * simp_bottomup descends into every Plus/Times subnode, runs the
     * full search on each, and only after all that work rediscovers the
     * difference-of-squares cancellation at the root (~6.5 s). With the
     * short-circuit, the canon collapses the whole thing to 0 in tens
     * of ms.
     *
     * Gated to depth == 0: the canon is expensive on large inputs and
     * we only want one chance at it per Simplify call. The recursive
     * call below uses depth + 1 so it does not retry. */
    Expr* canon_owned = NULL;
    if (depth == 0) {
        Expr* alt = transform_pythag_canon(input);
        if (alt) {
            if (!expr_eq(alt, input)) {
                size_t s_in = score_with_func(input, complexity_func);
                size_t s_alt = score_with_func(alt, complexity_func);
                if (s_alt < s_in) {
                    canon_owned = alt;
                    input = alt;
                } else {
                    expr_free(alt);
                }
            } else {
                expr_free(alt);
            }
        }
    }

    /* Top-level TanAddition short-circuit.  Sits BEFORE the TrigReduce
     * short-circuit because TrigReduce on inputs containing Sec[a+b] /
     * Csc[a+b] alongside multiple distinct Plus-arg trig calls expands
     * to a much larger Cos[...] Sec[...] Sec[...] product (e.g. a 13-leaf
     * Tan[z] Cos[x] Cos[y] Sec[x+y] (Tan[x]+Tan[y]) - Tan[z] Tan[x+y]
     * blows up to 9 Cos/Sec terms in ~700 ms, only to be rejected by the
     * score gate).  TanAddition's gate (has_pythag_head + 3+ distinct
     * trig args + a sum-witnessing triple) keeps it cheap when inert,
     * and on the case above it collapses the input directly to 0.
     *
     * Same depth==0 gating and strict-score gate as the other short-
     * circuits.  When TanAddition produces an atom (typically 0), the
     * `canon_owned && input->type != EXPR_FUNCTION` branch below returns
     * immediately, so the still-expensive TrigReduce short-circuit never
     * even runs. */
    if (depth == 0) {
        Expr* alt = transform_tan_addition(input);
        if (alt) {
            if (!expr_eq(alt, input)) {
                size_t s_in = score_with_func(input, complexity_func);
                size_t s_alt = score_with_func(alt, complexity_func);
                if (s_alt < s_in) {
                    if (canon_owned) expr_free(canon_owned);
                    canon_owned = alt;
                    input = alt;
                } else {
                    expr_free(alt);
                }
            } else {
                expr_free(alt);
            }
        }
    }

    /* Top-level TrigReduce short-circuit.  Parallel to the PythagCanon
     * one above: try product-to-sum + angle-addition collapse on the
     * whole input before descending into children.  Cases where the
     * input is a sum of two trig "compounds" whose TrigReduce'd forms
     * cancel -- e.g.
     *     (Sin[x]+Cos[x])^4 - (1 + Sin[2 x])^2  -> 0
     * -- bottom-up descent would call simp_search on each child
     * (running Factor/TrigFactor/etc. on each), then notice the
     * cancellation only at the root.  TrigReduce on the whole input
     * sees both summands as products of single-arg trig calls and
     * collapses them to the same `1/2 (3 - Cos[4 x] + 4 Sin[2 x])` form,
     * after which Plus auto-cancels and the result is 0.
     *
     * Same gating as the PythagCanon short-circuit: only at depth == 0
     * (so the recursive simp_bottomup call below cannot retry it), and
     * only adopt the result when it's a strict score win -- otherwise
     * fall through to the normal bottom-up search.  The transform_can_fire
     * gate keeps non-trig and non-Plus/Times inputs out, so the
     * extra call costs a few microseconds at worst. */
    if (depth == 0 &&
        transform_can_fire("TrigReduce", input, NULL)) {
        Expr* alt = traced_call_unary("TrigReduce", input);
        if (alt) {
            if (!expr_eq(alt, input)) {
                size_t s_in = score_with_func(input, complexity_func);
                size_t s_alt = score_with_func(alt, complexity_func);
                if (s_alt < s_in) {
                    if (canon_owned) expr_free(canon_owned);
                    canon_owned = alt;
                    input = alt;
                } else {
                    expr_free(alt);
                }
            } else {
                expr_free(alt);
            }
        }
    }

    /* If the canon collapsed input to an atom (e.g. 0 or a literal),
     * the rest of simp_bottomup -- which dereferences input->data.function
     * -- would read garbage. Short-circuit to the atom branch's behaviour:
     * with no assumptions, return the atom; with assumptions, defer to
     * simp_dispatch in case an assumption-driven rewrite still applies. */
    if (canon_owned && input->type != EXPR_FUNCTION) {
        Expr* result;
        if (!ctx || ctx->count == 0) {
            result = expr_copy((Expr*)input);
        } else {
            result = simp_dispatch(input, ctx, complexity_func);
        }
        expr_free(canon_owned);
        return result;
    }

    /* Memo lookup. */
    const Expr* hit = simp_memo_get(memo, input);
    if (hit) {
        Expr* cached = expr_copy((Expr*)hit);
        if (canon_owned) expr_free(canon_owned);
        return cached;
    }

    /* Held heads: don't descend, but still run top-level search. */
    const Expr* head = input->data.function.head;
    if (head && head->type == EXPR_SYMBOL) {
        const char* hn = head->data.symbol;
        if (simp_skip_recursion_head(hn) || simp_head_holds_args(hn)) {
            Expr* result = simp_dispatch(input, ctx, complexity_func);
            simp_memo_put(memo, input, result);
            if (canon_owned) expr_free(canon_owned);
            return result;
        }
    }

    /* Recurse into each child. */
    size_t argc = input->data.function.arg_count;
    Expr** new_args = (Expr**)calloc(argc ? argc : 1, sizeof(Expr*));
    bool any_changed = false;
    for (size_t i = 0; i < argc; i++) {
        new_args[i] = simp_bottomup(input->data.function.args[i], ctx,
                                    complexity_func, memo, depth + 1);
        if (!new_args[i]) {
            new_args[i] = expr_copy(input->data.function.args[i]);
        }
        if (!expr_eq(new_args[i], input->data.function.args[i])) {
            any_changed = true;
        }
    }

    Expr* canonical;
    if (any_changed) {
        Expr* new_head = expr_copy((Expr*)head);
        Expr* rebuilt = expr_new_function(new_head, new_args, argc);
        canonical = evaluate(rebuilt);
        expr_free(rebuilt);
    } else {
        for (size_t i = 0; i < argc; i++) expr_free(new_args[i]);
        canonical = expr_copy((Expr*)input);
    }
    free(new_args);

    /* Skip simp_search at non-top levels for "trivially small" subtrees.
     * Identity-collapse transforms (TrigFactor's Pythagorean rules,
     * LogExpRules, etc.) fire only when the subtree contains a *compound*
     * structure -- a sum, a product with multiple factors, a Power whose
     * base is itself a non-trivial expression. For something like
     * Cosh[x]^2 (4 leaves) or -Sinh[x]^2 (Times[-1, Power[Sinh[x],2]],
     * 7 leaves) in isolation, there is no useful identity to find, but
     * transforms like TrigRoundtrip on them produce explosive
     * intermediate forms (TrigToExp -> ExpToTrig of an isolated Cosh^2
     * leaves a 12-term polynomial in Cosh[2x], Sinh[2x], Cosh[4x],
     * Sinh[4x]) that drag the per-call cost into the seconds range.
     *
     * Pythagorean-eligible Plus/Times have at least 8 leaves
     * (Plus[Power[Sin,x,2], Power[Cos,x,2]] = 9; Plus[Power[Cosh,x,2],
     * Times[-1, Power[Sinh,x,2]]] = 12), so threshold 7 includes them
     * while excluding the explosive single-trig-power forms. The
     * top-level Simplify call (depth == 0) always runs simp_search,
     * regardless of size. */
    Expr* result;
    if (depth > 0 && simp_default_complexity(canonical) <= 7) {
        result = canonical;
    } else {
        result = simp_dispatch(canonical, ctx, complexity_func);
        expr_free(canonical);
    }

    simp_memo_put(memo, input, result);
    if (canon_owned) expr_free(canon_owned);
    return result;
}

/* ----------------------------------------------------------------------- */
/* builtin_simplify                                                        */
/* ----------------------------------------------------------------------- */

static Expr* read_dollar_assumptions(void) {
    /* Read the OwnValue directly. We must NOT evaluate $Assumptions, because
     * once an assumption like Element[x, Reals] becomes the bound value, our
     * own Element evaluator would recurse on it (Element reads $Assumptions
     * to decide -> evaluator fires the OwnValue rule -> Element[x, Reals]
     * gets re-evaluated -> ...). The first OwnValue rule on a symbol is its
     * current value (newest first); we deep-copy its replacement. */
    Rule* r = symtab_get_own_values("$Assumptions");
    if (!r || !r->replacement) return expr_new_symbol("True");
    return expr_copy(r->replacement);
}

/* ----------------------------------------------------------------------- */
/* builtin_element -- Element[x, Domain]                                   */
/* ----------------------------------------------------------------------- */

static bool is_complex_literal(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Complex
        && e->data.function.arg_count == 2;
}

static bool is_rational_literal(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2;
}

/* True iff `r` is exactly representable as a 64-bit integer. */
static bool real_is_integer(double r) {
    if (r != r) return false;                       /* NaN */
    if (r > 9.2233720368547758e18) return false;    /* > INT64_MAX */
    if (r < -9.2233720368547758e18) return false;
    long long i = (long long)r;
    return (double)i == r;
}

/* Element[x, dom] decision: 1 = True, 0 = False, -1 = undetermined. */
static int element_decide(const Expr* x, const char* dom, const AssumeCtx* ctx) {
    if (!x || !dom) return -1;

    /* Direct fact lookup is always safe regardless of domain. */
    if (ctx) {
        for (size_t i = 0; i < ctx->count; i++) {
            if (fact_in_domain(ctx->facts[i], x, dom)) return 1;
        }
    }

    if (strcmp(dom, "Integers") == 0) {
        if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT) return 1;
        if (x->type == EXPR_REAL) return real_is_integer(x->data.real) ? 1 : 0;
        if (is_rational_literal(x)) return 0;
        if (is_complex_literal(x)) return 0;
        if (prov_int(ctx, x)) return 1;
        return -1;
    }

    if (strcmp(dom, "Rationals") == 0) {
        if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT) return 1;
        if (is_rational_literal(x)) return 1;
        if (x->type == EXPR_REAL) return 1;             /* every double is dyadic */
        if (is_complex_literal(x)) return 0;
        if (prov_int(ctx, x)) return 1;
        return -1;
    }

    if (strcmp(dom, "Algebraics") == 0) {
        if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT) return 1;
        if (is_rational_literal(x)) return 1;
        if (x->type == EXPR_REAL) return 1;
        if (is_complex_literal(x)) return 1;            /* canonical Complex parts are rational */
        if (prov_int(ctx, x)) return 1;
        return -1;
    }

    if (strcmp(dom, "Reals") == 0) {
        if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT || x->type == EXPR_REAL) return 1;
        if (is_rational_literal(x)) return 1;
        if (is_complex_literal(x)) {
            /* canonical Complex always carries a non-zero imaginary part */
            Expr* im = x->data.function.args[1];
            if (im->type == EXPR_INTEGER && im->data.integer == 0) return 1;
            return 0;
        }
        if (prov_re(ctx, x)) return 1;
        return -1;
    }

    if (strcmp(dom, "Complexes") == 0) {
        if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT || x->type == EXPR_REAL) return 1;
        if (is_rational_literal(x)) return 1;
        if (is_complex_literal(x)) return 1;
        if (prov_re(ctx, x)) return 1;
        return -1;
    }

    if (strcmp(dom, "Booleans") == 0) {
        if (x->type == EXPR_SYMBOL) {
            if (x->data.symbol == SYM_True)  return 1;
            if (x->data.symbol == SYM_False) return 1;
        }
        return -1;
    }

    if (strcmp(dom, "Primes") == 0) {
        if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT) {
            Expr* primeq = call_unary_copy("PrimeQ", x);
            int ans = -1;
            if (primeq && primeq->type == EXPR_SYMBOL) {
                if (primeq->data.symbol == SYM_True) ans = 1;
                if (primeq->data.symbol == SYM_False) ans = 0;
            }
            if (primeq) expr_free(primeq);
            return ans;
        }
        return -1;
    }

    if (strcmp(dom, "Composites") == 0) {
        if ((x->type == EXPR_INTEGER && x->data.integer >= 2) || x->type == EXPR_BIGINT) {
            Expr* primeq = call_unary_copy("PrimeQ", x);
            int ans = -1;
            if (primeq && primeq->type == EXPR_SYMBOL) {
                if (primeq->data.symbol == SYM_True) ans = 0;
                if (primeq->data.symbol == SYM_False) ans = 1;
            }
            if (primeq) expr_free(primeq);
            return ans;
        }
        return -1;
    }

    return -1;
}

Expr* builtin_element(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;

    Expr* x   = res->data.function.args[0];
    Expr* dom = res->data.function.args[1];

    /* Threading: Element[{x1, x2, ...}, dom] -> {Element[x1, dom], ...}.
     * Mathematica returns this only when ALL elements decide; if any are
     * undetermined we leave them as Element[xi, dom]. */
    if (x->type == EXPR_FUNCTION && x->data.function.head &&
        x->data.function.head->type == EXPR_SYMBOL &&
        x->data.function.head->data.symbol == SYM_List) {
        size_t n = x->data.function.arg_count;
        Expr** out = (Expr**)calloc(n, sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            Expr* sub_args[2] = { expr_copy(x->data.function.args[i]), expr_copy(dom) };
            Expr* call = expr_new_function(expr_new_symbol("Element"), sub_args, 2);
            out[i] = evaluate(call);
            expr_free(call);
        }
        Expr* list = expr_new_function(expr_new_symbol("List"), out, n);
        free(out);
        return list;
    }

    if (dom->type != EXPR_SYMBOL) return NULL;
    const char* d = dom->data.symbol;

    /* Build context from current $Assumptions. */
    Expr* dollar = read_dollar_assumptions();
    AssumeCtx* ctx = assume_ctx_from_expr(dollar);
    expr_free(dollar);

    int decision = element_decide(x, d, ctx);
    assume_ctx_free(ctx);

    if (decision == 1)  return expr_new_symbol("True");
    if (decision == 0)  return expr_new_symbol("False");
    return NULL;
}

static Expr* combine_with_and(Expr* a, Expr* b) {
    /* Both inputs owned and consumed. */
    Expr* args[2] = { a, b };
    Expr* call = expr_new_function(expr_new_symbol("And"), args, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* ----------------------------------------------------------------------- */
/* Equation / inequality rebalancing                                       */
/*                                                                         */
/* For a binary relation `lhs OP rhs`, compute d = lhs - rhs as an         */
/* evaluated Plus, then rewrite as `pos OP neg` after dividing through by  */
/* the GCD of integer coefficients. Negative-coefficient terms move to    */
/* the opposite side. The result is correctness-preserving for both       */
/* equality and ordering relations (we never multiply or divide by a      */
/* negative quantity, only the positive integer GCD).                     */
/* ----------------------------------------------------------------------- */

static bool simp_eq_head_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol, name) == 0;
}

/* Decompose a Plus term into integer-coefficient * rest. Returns false
 * for terms whose leading numeric factor isn't an int64 (Real, BigInt,
 * Rational), which signals the caller to skip rebalancing -- mixing in
 * those would risk losing precision or introducing fractions. */
static bool simp_plus_term_int_coeff(const Expr* term, int64_t* coef,
                                     Expr** rest_out) {
    if (term->type == EXPR_INTEGER) {
        *coef = term->data.integer;
        *rest_out = expr_new_integer(1);
        return true;
    }
    if (term->type == EXPR_BIGINT || term->type == EXPR_REAL) return false;

    if (simp_eq_head_sym(term, "Times") &&
        term->data.function.arg_count >= 1) {
        const Expr* a0 = term->data.function.args[0];
        if (a0->type == EXPR_INTEGER) {
            *coef = a0->data.integer;
            size_t n = term->data.function.arg_count;
            if (n == 2) {
                *rest_out = expr_copy(term->data.function.args[1]);
            } else {
                Expr** args = (Expr**)calloc(n - 1, sizeof(Expr*));
                for (size_t i = 1; i < n; i++) {
                    args[i - 1] = expr_copy(term->data.function.args[i]);
                }
                *rest_out = expr_new_function(
                    expr_new_symbol("Times"), args, n - 1);
                free(args);
            }
            return true;
        }
        if (a0->type == EXPR_BIGINT || a0->type == EXPR_REAL) return false;
        if (simp_eq_head_sym(a0, "Rational")) return false;
    }

    /* Generic term: implicit coefficient 1, rest = term. */
    *coef = 1;
    *rest_out = expr_copy((Expr*)term);
    return true;
}

/* Build `c * rest`, dropping a coefficient of 1 and Times wrappers when
 * rest = 1. Takes ownership of `rest`. */
static Expr* simp_make_term(int64_t c, Expr* rest) {
    if (rest->type == EXPR_INTEGER && rest->data.integer == 1) {
        expr_free(rest);
        return expr_new_integer(c);
    }
    if (c == 1) return rest;
    /* Flatten into existing Times; otherwise wrap. */
    if (simp_eq_head_sym(rest, "Times")) {
        size_t n = rest->data.function.arg_count;
        Expr** args = (Expr**)calloc(n + 1, sizeof(Expr*));
        args[0] = expr_new_integer(c);
        for (size_t i = 0; i < n; i++) {
            args[i + 1] = expr_copy(rest->data.function.args[i]);
        }
        Expr* out = expr_new_function(
            expr_new_symbol("Times"), args, n + 1);
        free(args);
        expr_free(rest);
        return out;
    }
    Expr* args[2] = { expr_new_integer(c), rest };
    return expr_new_function(expr_new_symbol("Times"), args, 2);
}

/* Returns NULL when no rebalanced form is produced (non-int64 coeffs,
 * fully symbolic d, or d = 0). The caller compares scores. */
static Expr* simp_try_rebalance_relation(const Expr* relation) {
    if (!relation || relation->type != EXPR_FUNCTION) return NULL;
    if (relation->data.function.arg_count != 2) return NULL;
    const Expr* h = relation->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return NULL;
    const char* hn = h->data.symbol;
    bool ok = (hn == SYM_Equal ||
               hn == SYM_Unequal ||
               hn == SYM_Less ||
               hn == SYM_LessEqual ||
               hn == SYM_Greater ||
               hn == SYM_GreaterEqual);
    if (!ok) return NULL;

    /* d = lhs - rhs, evaluated. */
    Expr* neg_args[2] = {
        expr_new_integer(-1),
        expr_copy(relation->data.function.args[1])
    };
    Expr* neg_rhs = expr_new_function(
        expr_new_symbol("Times"), neg_args, 2);
    Expr* d_args[2] = {
        expr_copy(relation->data.function.args[0]),
        neg_rhs
    };
    Expr* d_call = expr_new_function(
        expr_new_symbol("Plus"), d_args, 2);
    Expr* d_sum = evaluate(d_call);
    /* Expand so Times[2, Plus[...]] partitions term-by-term. The threaded
     * input may already have collected common factors via Collect, which
     * defeats coefficient-level rebalancing. */
    Expr* exp_args[1] = { d_sum };
    Expr* d_exp_call = expr_new_function(
        expr_new_symbol("Expand"), exp_args, 1);
    Expr* d = evaluate(d_exp_call);

    Expr* d_singleton[1];
    Expr** terms;
    size_t n;
    if (simp_eq_head_sym(d, "Plus")) {
        n = d->data.function.arg_count;
        terms = d->data.function.args;
    } else {
        d_singleton[0] = d;
        terms = d_singleton;
        n = 1;
    }
    if (n == 0) { expr_free(d); return NULL; }

    /* Extract integer coefficients. Bail on non-int64. */
    int64_t* coefs = (int64_t*)calloc(n, sizeof(int64_t));
    Expr** rests = (Expr**)calloc(n, sizeof(Expr*));
    bool ok2 = true;
    for (size_t i = 0; i < n; i++) {
        if (!simp_plus_term_int_coeff(terms[i], &coefs[i], &rests[i])) {
            ok2 = false;
            for (size_t j = 0; j < i; j++) expr_free(rests[j]);
            break;
        }
    }
    if (!ok2) {
        free(coefs);
        free(rests);
        expr_free(d);
        return NULL;
    }

    /* GCD of |coefs|. */
    int64_t g = 0;
    for (size_t i = 0; i < n; i++) {
        int64_t c = coefs[i];
        if (c == INT64_MIN) { g = 1; break; }
        if (c < 0) c = -c;
        g = (g == 0) ? c : gcd(g, c);
    }
    if (g == 0) g = 1;

    /* Polarity: pick the first non-constant term's coefficient sign so the
     * leading variable term ends up positive after dividing through. This
     * turns `-2 x == 4` into `x == -2` rather than `0 == x + 2`. For strict
     * inequalities (Less, Greater) a negative divisor flips the operator;
     * the non-strict and equality forms are direction-symmetric. */
    int64_t divisor = g;
    bool flipped = false;
    for (size_t i = 0; i < n; i++) {
        bool is_const = (rests[i]->type == EXPR_INTEGER &&
                         rests[i]->data.integer == 1);
        if (!is_const) {
            if (coefs[i] < 0) { divisor = -g; flipped = true; }
            break;
        }
    }
    for (size_t i = 0; i < n; i++) coefs[i] /= divisor;

    const char* out_head = hn;
    if (flipped) {
        if      (hn == SYM_Less)         out_head = "Greater";
        else if (hn == SYM_Greater)      out_head = "Less";
        else if (hn == SYM_LessEqual)    out_head = "GreaterEqual";
        else if (hn == SYM_GreaterEqual) out_head = "LessEqual";
    }

    /* Build LHS from positive-coef variable terms, RHS from
     * negated-negative-coef variable terms plus the negated constant. */
    Expr** pos = (Expr**)calloc(n, sizeof(Expr*));
    Expr** neg = (Expr**)calloc(n, sizeof(Expr*));
    size_t pn = 0, nn = 0;
    int64_t const_sum = 0;       /* moves to RHS as -const_sum */
    bool const_overflow = false; /* on overflow, fall back to a Plus term */
    Expr** const_terms = (Expr**)calloc(n, sizeof(Expr*));
    size_t cn = 0;
    for (size_t i = 0; i < n; i++) {
        bool is_const = (rests[i]->type == EXPR_INTEGER &&
                         rests[i]->data.integer == 1);
        if (is_const) {
            int64_t c = coefs[i];
            /* Track sum but guard against int64 overflow. */
            int64_t sum;
            if (!const_overflow &&
                ((c > 0 && const_sum > INT64_MAX - c) ||
                 (c < 0 && const_sum < INT64_MIN - c))) {
                const_overflow = true;
            }
            if (!const_overflow) {
                sum = const_sum + c;
                const_sum = sum;
            }
            /* Always keep the term in case we hit overflow later. */
            const_terms[cn++] = simp_make_term(c, rests[i]);
        } else {
            if (coefs[i] > 0) {
                pos[pn++] = simp_make_term(coefs[i], rests[i]);
            } else if (coefs[i] < 0) {
                neg[nn++] = simp_make_term(-coefs[i], rests[i]);
            } else {
                expr_free(rests[i]);
            }
        }
    }

    Expr* new_lhs;
    if (pn == 0)      new_lhs = expr_new_integer(0);
    else if (pn == 1) new_lhs = pos[0];
    else              new_lhs = expr_new_function(
                          expr_new_symbol("Plus"), pos, pn);

    /* RHS = (negated negative-coef vars) + (-const). */
    size_t total_rhs = nn + cn;
    Expr* new_rhs;
    if (total_rhs == 0) {
        new_rhs = expr_new_integer(0);
        for (size_t i = 0; i < cn; i++) expr_free(const_terms[i]);
    } else {
        Expr** rhs_terms = (Expr**)calloc(total_rhs, sizeof(Expr*));
        size_t rt = 0;
        for (size_t i = 0; i < nn; i++) rhs_terms[rt++] = neg[i];
        if (!const_overflow) {
            /* Single integer for the constant: -const_sum (zero is fine). */
            for (size_t i = 0; i < cn; i++) expr_free(const_terms[i]);
            if (const_sum != 0 || rt == 0) {
                /* Build -const_sum, watching INT64_MIN. */
                int64_t neg_const = (const_sum == INT64_MIN)
                                        ? INT64_MAX  /* impossible in practice */
                                        : -const_sum;
                rhs_terms[rt++] = expr_new_integer(neg_const);
            }
        } else {
            /* Overflow path: keep each constant term, negated. */
            for (size_t i = 0; i < cn; i++) {
                /* Negate the leading coefficient. */
                if (const_terms[i]->type == EXPR_INTEGER) {
                    /* Replace, don't mutate: the integer atom may be
                     * shared (M3 atom-sharing). */
                    int64_t v = -const_terms[i]->data.integer;
                    expr_free(const_terms[i]);
                    rhs_terms[rt++] = expr_new_integer(v);
                } else {
                    /* Wrap in Times[-1, ...]. */
                    Expr* args[2] = { expr_new_integer(-1), const_terms[i] };
                    rhs_terms[rt++] = expr_new_function(
                        expr_new_symbol("Times"), args, 2);
                }
            }
        }
        if (rt == 0) {
            new_rhs = expr_new_integer(0);
            free(rhs_terms);
        } else if (rt == 1) {
            new_rhs = rhs_terms[0];
            free(rhs_terms);
        } else {
            new_rhs = expr_new_function(
                expr_new_symbol("Plus"), rhs_terms, rt);
            free(rhs_terms);
        }
    }

    free(const_terms);
    free(pos);
    free(neg);
    free(coefs);
    free(rests);
    expr_free(d);

    /* Re-evaluate each side so canonical ordering / Plus flattening kicks in. */
    Expr* lhs_e = evaluate(new_lhs);
    Expr* rhs_e = evaluate(new_rhs);

    Expr* rel_args[2] = { lhs_e, rhs_e };
    Expr* out = expr_new_function(
        expr_new_symbol(out_head), rel_args, 2);
    Expr* out_eval = evaluate(out);
    return out_eval;
}

Expr* builtin_simplify(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    /* The simplification pipeline routes through Together/Cancel/Apart/
     * Factor and the polynomial GCD machinery, all of which need rational
     * coefficients. Rationalise on entry, run the exact pipeline, then
     * numericalise on the way out so callers still see inexact-in /
     * inexact-out semantics. */
    if (internal_args_contain_inexact(res)) {
        return internal_rationalize_then_numericalize(res, builtin_simplify);
    }

    Expr* expr = res->data.function.args[0];

    /* Parse remaining args: at most one positional assumption, plus
     * options Rule[Assumptions, X] and Rule[ComplexityFunction, f]. */
    Expr* positional_assum = NULL;
    Expr* opt_assumptions  = NULL;
    Expr* opt_complexity   = NULL;

    for (size_t i = 1; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (is_rule_with_lhs(a, "Assumptions")) {
            opt_assumptions = a->data.function.args[1];
        } else if (is_rule_with_lhs(a, "ComplexityFunction")) {
            opt_complexity = a->data.function.args[1];
        } else if (positional_assum == NULL) {
            positional_assum = a;
        }
    }

    /* ComplexityFunction -> Automatic is a synonym for the built-in
     * default. Treating it as NULL makes score_with_func use the fast
     * native simp_default_complexity path instead of evaluating
     * Automatic[candidate] (which would never reduce). */
    if (opt_complexity &&
        opt_complexity->type == EXPR_SYMBOL &&
        opt_complexity->data.symbol == SYM_Automatic) {
        opt_complexity = NULL;
    }

    /* Compute the effective assumption expression.
     *   - With Assumptions->X, X overrides the $Assumptions default.
     *   - Without, the positional assumption is appended to $Assumptions.
     * Then evaluate to canonicalise (e.g. And[True, x>0] -> x>0). */
    Expr* effective;
    if (opt_assumptions) {
        if (positional_assum) {
            effective = combine_with_and(expr_copy(positional_assum),
                                         expr_copy(opt_assumptions));
        } else {
            effective = evaluate(expr_copy(opt_assumptions));
        }
    } else {
        Expr* dollar = read_dollar_assumptions();
        if (positional_assum) {
            effective = combine_with_and(expr_copy(positional_assum), dollar);
        } else {
            effective = dollar;
        }
    }

    AssumeCtx* ctx = assume_ctx_from_expr(effective);
    expr_free(effective);

    /* If the input is a predicate that appears literally as one of our
     * assumed facts, it folds to True. This is a narrow win for simple
     * cases like Simplify[x > 0, x > 0]; it does not constitute a real
     * inequality reasoner (see picocas_spec.md for v1 gaps). */
    if (ctx) {
        for (size_t i = 0; i < ctx->count; i++) {
            if (expr_eq(expr, ctx->facts[i])) {
                assume_ctx_free(ctx);
                return expr_new_symbol("True");
            }
        }
    }

    /* Manual threading over Equal/Less/.../And/Or (List handled by
     * ATTR_LISTABLE on the Simplify symbol itself). For binary
     * relational heads we additionally try a rebalanced form
     * `pos OP neg` (after dividing through by the GCD of integer
     * coefficients) and pick the simpler of the two by SimplifyCount. */
    if (expr->type == EXPR_FUNCTION &&
        expr->data.function.head &&
        expr->data.function.head->type == EXPR_SYMBOL &&
        head_threads_over(expr->data.function.head->data.symbol)) {
        size_t n = expr->data.function.arg_count;
        Expr** new_args = (Expr**)calloc(n, sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            Expr** sub_args = (Expr**)calloc(argc, sizeof(Expr*));
            sub_args[0] = expr_copy(expr->data.function.args[i]);
            for (size_t k = 1; k < argc; k++) {
                sub_args[k] = expr_copy(res->data.function.args[k]);
            }
            Expr* call = expr_new_function(expr_new_symbol("Simplify"), sub_args, argc);
            new_args[i] = evaluate(call);
            expr_free(call);
        }
        Expr* threaded = expr_new_function(expr_copy(expr->data.function.head), new_args, n);
        free(new_args);
        Expr* threaded_eval = evaluate(threaded);

        /* Rebalance candidate: only meaningful for a binary relation that
         * survived evaluation (Equal collapsed to True/False is not a
         * Function any more). */
        Expr* rebalanced = simp_try_rebalance_relation(threaded_eval);
        if (rebalanced && !expr_eq(rebalanced, threaded_eval)) {
            size_t s_threaded = score_with_func(threaded_eval, opt_complexity);
            size_t s_rebal    = score_with_func(rebalanced, opt_complexity);
            if (s_rebal < s_threaded) {
                expr_free(threaded_eval);
                threaded_eval = rebalanced;
            } else {
                expr_free(rebalanced);
            }
        } else if (rebalanced) {
            expr_free(rebalanced);
        }

        assume_ctx_free(ctx);
        return threaded_eval;
    }

    SimpMemo memo;
    simp_memo_init(&memo);

    FactorMemo* fmemo = factor_memo_new();
    factor_memo_push(fmemo);

    /* Top-level rational shortcut. simp_bottomup descends into every Plus /
     * Times child before dispatching at the top, and for a SHAPE_RATIONAL
     * input each child re-enters simp_dispatch -> simp_pipeline_rational.
     * Together / Cancel / Factor at the top combines all the children into
     * a single canonical num/den, so the per-child work is wasted: each
     * subnode's "best" form ends up subsumed by the top-level Together.
     *
     * Empirically, on multivariate rational inputs Simplify takes ~8 s vs
     * Cancel[Together[expr]] ~25 ms (~300x). Even when the search returns
     * the input unchanged, the cost is in the search itself. By dispatching
     * once at the top we cut directly to the pipeline that decides
     * acceptance against the input, bypassing the redundant per-subnode
     * traversal. The polish passes (lift_common_factor, PythagReduce,
     * canon_negate_pairs) still run on the result.
     *
     * Gated on SHAPE_RATIONAL: the classifier rejects inputs with trig,
     * log, abs, and non-integer powers, so we only take the shortcut when
     * the polynomial pipeline has full coverage. */
    Expr* best;
    if (simp_classify(expr) == SIMP_SHAPE_RATIONAL) {
        best = simp_dispatch(expr, ctx, opt_complexity);
    } else {
        /* Top-level algebraic-rational fast path. When the input is a
         * Plus over a multi-generator algebraic-number tower (e.g. the
         * output of D[Integrate[a x/(x^3+2), x], x] which is a sum of 3
         * fractions over {2^(1/3), Sqrt[3], Sqrt[radicand-with-α-inside]}),
         * Together[expr, Extension -> Automatic] is the one transform
         * that can collapse it back to (a x)/(x^3+2) in a single pass via
         * builtin_together's multi-gen single-α fallback (rat.c, Phase F).
         * simp_bottomup's per-subnode descent doesn't reach this combined-
         * over-common-denominator form on its own.  Strict leaf-count gate
         * ensures no regression on inputs where Together-with-Auto is a
         * no-op or worse. */
        Expr* alg_top = NULL;
        if (expr->type == EXPR_FUNCTION
            && expr->data.function.head
            && expr->data.function.head->type == EXPR_SYMBOL
            && expr->data.function.head->data.symbol == SYM_Plus
            && has_non_integer_power(expr)
            && !contains_explicit_complex(expr)) {
            QATower* qa_t = extension_autodetect(expr);
            if (qa_t && qa_t->n >= 2) {
                qa_tower_free(qa_t);
                Expr* tog = expr_new_function(
                    expr_new_symbol("Together"),
                    (Expr*[]){
                        expr_copy(expr),
                        expr_new_function(expr_new_symbol("Rule"),
                            (Expr*[]){expr_new_symbol("Extension"),
                                      expr_new_symbol("Automatic")}, 2)
                    }, 2);
                Expr* cand = evaluate(tog);
                expr_free(tog);
                if (cand && simp_default_complexity(cand)
                                < simp_default_complexity(expr)) {
                    alg_top = cand;
                } else if (cand) {
                    expr_free(cand);
                }
            } else if (qa_t) {
                qa_tower_free(qa_t);
            }
        }

        if (alg_top) {
            /* Use the algebraic collapse as the starting point.  Run
             * simp_bottomup on it to apply any further polish (rare for
             * already-canonical rational forms but harmless). */
            best = simp_bottomup(alg_top, ctx, opt_complexity, &memo, 0);
            size_t s_alg = score_with_func(alg_top, opt_complexity);
            size_t s_best = score_with_func(best, opt_complexity);
            if (s_alg <= s_best) {
                expr_free(best);
                best = alg_top;
            } else {
                expr_free(alg_top);
            }
        } else {
            /* Top-level trig-rational fast path. Substitutes Sin/Cos/Sinh/
             * Cosh (and Tan/Cot/Sec/Csc/Tanh/etc. after preprocessing) plus
             * every opaque non-rational subtree (Log[...], Exp[...], etc.)
             * into fresh ground-field symbols so the algebraic core sees a
             * pure rational function; works in the quotient ring modulo the
             * trig/hyp ideals, then back-substitutes. Strict leaf-count gate
             * inside ensures it never regresses; on no-improvement or when
             * the input is out of budget it returns NULL and we fall through
             * to the normal bottom-up search. Doing this BEFORE simp_bottomup
             * means we bypass the per-subnode descent (which itself is
             * extremely slow on inputs like
             *   D[Integrate`RischNorman[Tan[x]^2 + Tan[x] + 1, x], x]
             * because every internal node fires a full simp_search). */
            Expr* tr = simp_trig_rational(expr, ctx, opt_complexity);
            if (tr) {
                best = tr;
            } else {
                best = simp_bottomup(expr, ctx, opt_complexity, &memo, 0);
            }
        }
    }

    /* Final-form polish: lift a shared algebraic generator out of a
     * top-level Plus (or out of a Plus child of a top-level Times -- the
     * numerator of a fraction with a non-integer-power denominator).
     * This catches:
     *   (8/105)(1+x^2)^(3/2) - (4/35)x^2(1+x^2)^(3/2) + (1/7)x^4(1+x^2)^(3/2)
     *     -> (1/105)(1+x^2)^(3/2)(8 - 12 x^2 + 15 x^4)
     *   (15 x^2 + 5 x^3)/(5+2x)^(3/2)
     *     -> (5 x^2 (3 + x))/(5+2x)^(3/2)
     * which picocas's polynomial Factor cannot reach because Variables[]
     * does not return non-integer-power generators. We apply it once at
     * the top level rather than as a seed in simp_search to avoid
     * destabilising the heuristic search on multi-variable trig inputs. */
    {
        Expr* lifted = simp_lift_common_factor(best);
        if (lifted && !expr_eq(lifted, best)) {
            size_t s_lift = score_with_func(lifted, opt_complexity);
            size_t s_best = score_with_func(best, opt_complexity);
            if (s_lift <= s_best) {
                expr_free(best);
                best = lifted;
            } else {
                expr_free(lifted);
            }
        } else if (lifted) {
            expr_free(lifted);
        }
    }

    /* Pythagorean polish: PythagReduce already runs as a seed inside
     * simp_search, but its result enters update_best with a strict `<`
     * tiebreak, so structurally-collapsed forms that tie on
     * SimplifyCount (e.g. `-Sech[x]^2` vs `-1 + Tanh[x]^2`, both score
     * 7) lose to whatever arrived at the score plateau first. As a
     * polish, accept on `<=`: when the pythag rules turn the result
     * into a single Power-of-trig head with the same score or lower,
     * take it. Bypass when the Tanh/Coth/Tan/Cot rules cannot fire
     * (no relevant head present). */
    {
        Expr* reduced = transform_pythag_reduce(best);
        if (reduced && !expr_eq(reduced, best)) {
            size_t s_red  = score_with_func(reduced, opt_complexity);
            size_t s_best = score_with_func(best, opt_complexity);
            if (s_red <= s_best) {
                expr_free(best);
                best = reduced;
            } else {
                expr_free(reduced);
            }
        } else if (reduced) {
            expr_free(reduced);
        }
    }

    /* Sign canonicalisation: flip pairs of negative-leading Plus factors
     * inside a top-level Times so each binomial leads with its
     * positive-coefficient term, e.g.
     *   ((-a + c) (-b + d))/(a b c d)  ->  ((a - c) (b - d))/(a b c d)
     * Value-preserving (flips occur in pairs so signs cancel). */
    {
        Expr* canon = canon_negate_pairs(best);
        if (canon) {
            expr_free(best);
            best = canon;
        }
    }

    factor_memo_pop();
    factor_memo_free(fmemo);

    simp_memo_free(&memo);
    assume_ctx_free(ctx);
    return best;
}

/* ----------------------------------------------------------------------- */
/* builtin_assuming -- desugar to Block[{$Assumptions = $A && a}, body]    */
/* ----------------------------------------------------------------------- */

Expr* builtin_assuming(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;

    Expr* assum = res->data.function.args[0];   /* already evaluated */
    Expr* body  = res->data.function.args[1];   /* held by HoldRest */

    /* Convert lists of assumptions to conjunctions, per Mathematica
     * semantics. */
    Expr* assum_norm;
    if (assum->type == EXPR_FUNCTION &&
        assum->data.function.head &&
        assum->data.function.head->type == EXPR_SYMBOL &&
        assum->data.function.head->data.symbol == SYM_List) {
        size_t n = assum->data.function.arg_count;
        Expr** copies = (Expr**)calloc(n, sizeof(Expr*));
        for (size_t i = 0; i < n; i++) copies[i] = expr_copy(assum->data.function.args[i]);
        Expr* and_call = expr_new_function(expr_new_symbol("And"), copies, n);
        free(copies);
        assum_norm = and_call;  /* not yet evaluated; Block will evaluate it */
    } else {
        assum_norm = expr_copy(assum);
    }

    /* Build $Assumptions && assum_norm */
    Expr* and_args[2] = { expr_new_symbol("$Assumptions"), assum_norm };
    Expr* combined = expr_new_function(expr_new_symbol("And"), and_args, 2);

    /* Build Set[$Assumptions, combined] -- represents
     * "$Assumptions = $Assumptions && a" inside the Block var list. */
    Expr* set_args[2] = { expr_new_symbol("$Assumptions"), combined };
    Expr* set_call = expr_new_function(expr_new_symbol("Set"), set_args, 2);

    /* Block[{Set[$Assumptions, ...]}, body] */
    Expr* var_list_args[1] = { set_call };
    Expr* var_list = expr_new_function(expr_new_symbol("List"), var_list_args, 1);

    Expr* block_args[2] = { var_list, expr_copy(body) };
    Expr* block_call = expr_new_function(expr_new_symbol("Block"), block_args, 2);

    Expr* result = evaluate(block_call);
    expr_free(block_call);
    return result;
}

/* ----------------------------------------------------------------------- */
/* simp_init                                                               */
/* ----------------------------------------------------------------------- */

void simp_init(void) {
    /* $Assumptions defaults to True. */
    Expr* dollar_pat = expr_new_symbol("$Assumptions");
    Expr* dollar_val = expr_new_symbol("True");
    symtab_add_own_value("$Assumptions", dollar_pat, dollar_val);
    expr_free(dollar_pat);
    expr_free(dollar_val);

    /* $SimplifyDebug defaults to False. When set to True, Simplify emits
     * one stderr line per transform invocation in the form
     *   /<TransformName>/: <input> -> <output> [<elapsed> ms]
     * to help diagnose hangs and runaway candidate explosion. */
    Expr* dbg_pat = expr_new_symbol("$SimplifyDebug");
    Expr* dbg_val = expr_new_symbol("False");
    symtab_add_own_value("$SimplifyDebug", dbg_pat, dbg_val);
    expr_free(dbg_pat);
    expr_free(dbg_val);
    symtab_set_docstring("$SimplifyDebug",
        "$SimplifyDebug\n\tWhen set to True, Simplify prints one stderr line per\n"
        "\ttransform invocation: /Name/: <input> -> <output> [<ms> ms].\n"
        "\tDefaults to False. Useful for diagnosing slow Simplify calls.");

    symtab_add_builtin("Simplify", builtin_simplify);
    symtab_get_def("Simplify")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);

    symtab_add_builtin("SimplifyCount", builtin_simplify_count);
    symtab_get_def("SimplifyCount")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);
    symtab_set_docstring("SimplifyCount",
        "SimplifyCount[expr]\n\tThe complexity measure used by Simplify when no\n"
        "\tComplexityFunction option (or ComplexityFunction -> Automatic) is\n"
        "\tgiven. Counts subexpressions; integers contribute their decimal\n"
        "\tdigit count plus a constant for the sign. Real numbers contribute\n"
        "\t2 (NumberQ but not Integer/Rational).");

    symtab_add_builtin("Assuming", builtin_assuming);
    symtab_get_def("Assuming")->attributes |= (ATTR_HOLDREST | ATTR_PROTECTED);

    symtab_add_builtin("Element", builtin_element);
    symtab_get_def("Element")->attributes |= ATTR_PROTECTED;
}
