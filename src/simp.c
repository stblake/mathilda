#include "simp.h"
#include "arithmetic.h"
#include "attr.h"
#include "eval.h"
#include "facpoly.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "expr.h"
#include "rationalize.h"
#include "sym_names.h"
#include "sym_intern.h"

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
 * the leaf level so we don't need a separate fact for them. v1 does not
 * propagate evenness through Times/Plus -- those would need a stronger
 * "two of these are even" reasoner. */
static bool prov_even(const AssumeCtx* ctx, const Expr* x) {
    if (!x) return false;
    if (x->type == EXPR_INTEGER) return (x->data.integer % 2) == 0;
    if (x->type == EXPR_BIGINT) return mpz_even_p(x->data.bigint) != 0;
    if (!ctx) return false;
    for (size_t i = 0; i < ctx->count; i++) {
        if (fact_implies_even(ctx->facts[i], x)) return true;
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
        if (fact_in_domain(ctx->facts[i], x, "Integers")) return true;
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

static size_t score_with_func(const Expr* e, const Expr* complexity_func) {
    if (!complexity_func) return simp_default_complexity(e);
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
    }

    for (size_t i = 0; i < nneg; i++) {
        const char* x = negatives[i];
        /* Abs[x] -> -x  for x < 0 */
        SEP(); EMIT("Abs[%s] :> -%s", x, x);
        /* Sqrt[x^2] -> -x  for x < 0  (Power[Power[x,2], 1/2]) */
        SEP(); EMIT("Power[Power[%s, 2], Rational[1, 2]] :> -%s", x, x);
        /* Sqrt[x^2 * rest] -> -x * Sqrt[rest] for x < 0. */
        SEP(); EMIT("Power[Times[Power[%s, 2], rest___], Rational[1, 2]] :> -%s Power[Times[rest], Rational[1, 2]]", x, x);
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
    if (is_rational((Expr*)e, &n, &d)) {
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

/* Cheap structural check: does `e` contain any of Cos/Sin/Cosh/Sinh as
 * a function head?  PythagReduce's rules cannot match anything else,
 * so when the answer is no we can skip the ReplaceRepeated walk
 * entirely.  Walks the tree once; cheaper by orders of magnitude than
 * the pattern-matching pass it gates. */
static bool has_pythag_head(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL) {
        const char* h = head->data.symbol;
        if (h == SYM_Cos || h == SYM_Sin ||
            h == SYM_Cosh || h == SYM_Sinh ||
            h == SYM_Tan || h == SYM_Cot ||
            h == SYM_Tanh || h == SYM_Coth) {
            return true;
        }
    }
    if (has_pythag_head(head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_pythag_head(e->data.function.args[i])) return true;
    }
    return false;
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
            "  -1 - Cot[x_]^2 + r___  :> -Csc[x]^2 + r }");
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

    /* Power[Times[u1,...,un], a]  -> Times[ui^a]  when every ui positive.
     * Power[Power[x, p], q]       -> Power[x, p*q] when x positive, p real. */
    if (h == SYM_Power && n == 2) {
        Expr* base = a[0];
        Expr* exp_  = a[1];

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
    if (!ctx) return NULL;
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
    /* LogExp identity cascade: nothing fires unless Log or Power is present
     * AND the assumption ctx has at least one fact (the cascade reads
     * positivity/reality of operands from ctx). */
    if (strcmp(name, "LogExpRules") == 0) {
        if (!ctx_has_facts(ctx)) return false;
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
 * external) if any step fails. */
static bool simp_is_rational(const Expr* e) {
    Expr* tg = call_unary_copy("Together", e);
    if (!tg) return false;
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

    /* Factor on the Cancel'd form (gated against non-integer powers, which
     * shouldn't appear on a SHAPE_RATIONAL input but is defensive). */
    if (cn && !has_non_integer_power(cn)) {
        Expr* fc = traced_call_unary("Factor", cn);
        if (fc) {
            update_best(&best, &bs, fc, complexity_func);
            expr_free(fc);
        }
    }

    /* Per-variable Collect on the Cancel'd form. */
    if (cn && transform_can_fire("CollectPerVariable", cn, NULL)) {
        CandSet next; cs_init(&next);
        try_collect_per_variable(cn, bs, &next, &best, &bs, complexity_func);
        cs_free(&next);
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
    /* Pairwise union by symbol intersection. Empty sets (constants)
     * stay in their own singleton component; we don't glom them onto
     * an arbitrary first match. */
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
        case SIMP_SHAPE_TRIG:
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
    }
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

    /* Memo lookup. */
    const Expr* hit = simp_memo_get(memo, input);
    if (hit) return expr_copy((Expr*)hit);

    /* Held heads: don't descend, but still run top-level search. */
    const Expr* head = input->data.function.head;
    if (head && head->type == EXPR_SYMBOL) {
        const char* hn = head->data.symbol;
        if (simp_skip_recursion_head(hn) || simp_head_holds_args(hn)) {
            Expr* result = simp_dispatch(input, ctx, complexity_func);
            simp_memo_put(memo, input, result);
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
        best = simp_bottomup(expr, ctx, opt_complexity, &memo, 0);
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
