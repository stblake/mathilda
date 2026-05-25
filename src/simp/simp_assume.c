#include "simp.h"
#include "simp_internal.h"
#include "arithmetic.h"
#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expand.h"
#include "facpoly.h"
#include "numeric.h"
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
#include "simp_log.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif


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
        /* Element[{x1, ..., xN}, dom] and Element[x1|x2|...|xN, dom] are
         * both shorthand for the conjunction of Element[xi, dom]. Split
         * into per-variable facts so prov_re / fact_in_domain see each
         * variable individually. */
        if (h == SYM_Element &&
            a->data.function.arg_count == 2 &&
            a->data.function.args[0]->type == EXPR_FUNCTION &&
            a->data.function.args[0]->data.function.head &&
            a->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
            (a->data.function.args[0]->data.function.head->data.symbol == SYM_List ||
             a->data.function.args[0]->data.function.head->data.symbol == SYM_Alternatives)) {
            const Expr* xs  = a->data.function.args[0];
            const Expr* dom = a->data.function.args[1];
            for (size_t i = 0; i < xs->data.function.arg_count; i++) {
                Expr* sub_args[2] = {
                    expr_copy(xs->data.function.args[i]),
                    expr_copy((Expr*)dom)
                };
                Expr* fact = expr_new_function(expr_new_symbol("Element"),
                                               sub_args, 2);
                ctx_push(ctx, fact);
                expr_free(fact);
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
int numeric_sign(const Expr* e) {
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

bool fact_is_function(const Expr* f, const char* head, size_t arity) {
    return head_is(f, intern_symbol(head))
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
bool fact_in_domain(const Expr* f, const Expr* x, const char* dom) {
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
bool is_positive_constant_symbol(const char* s) {
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
bool is_real_constant_symbol(const char* s) {
    if (is_positive_constant_symbol(s)) return true;
    return strcmp(s, "MachineEpsilon") == 0;
}

/* Forward declarations for mutual recursion across the predicate family. */
bool prov_pos(const AssumeCtx* ctx, const Expr* x);
bool prov_nn (const AssumeCtx* ctx, const Expr* x);
bool prov_neg(const AssumeCtx* ctx, const Expr* x);
bool prov_np (const AssumeCtx* ctx, const Expr* x);
bool prov_int(const AssumeCtx* ctx, const Expr* x);
bool prov_re (const AssumeCtx* ctx, const Expr* x);

/* Local helper: true iff `e` contains a free (unknown) symbol leaf as
 * an argument position — i.e. a symbol other than the named numeric
 * constants Pi / E / EulerGamma / GoldenRatio / Catalan / Degree /
 * Glaisher / Khinchin appearing somewhere a numeric value would be
 * required. Function heads (Cos, Sin, Log, ...) do NOT count even
 * though they are non-constant symbol leaves: `Cos[2]` is fully
 * numeric. This is the cheap gate before paying for `numericalize`. */
static bool nsf_has_free_symbol(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        return !is_real_constant_symbol(e->data.symbol);
    }
    if (e->type != EXPR_FUNCTION) return false;
    /* Skip the head — Cos/Sin/Log/etc. are not "free variables" even
     * though they aren't recognised numeric constants. */
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (nsf_has_free_symbol(e->data.function.args[i])) return true;
    }
    return false;
}

/* Numerical sign fallback used by prov_pos / prov_neg when structural
 * analysis can't decide. Returns +1 / -1 if N[x] is a finite real number
 * whose magnitude is far enough from zero to rule out accumulated MPFR
 * roundoff (we use ~1e-20 with ~100 bits of working precision, leaving
 * a ~5-decade safety margin). Returns 0 ("undecided") when:
 *   - x contains a free symbol (no numeric value to compute),
 *   - numericalize couldn't produce a numeric result,
 *   - the magnitude is within the safety margin (the expression may
 *     algebraically be exactly zero — caller must NOT assume a sign),
 *   - the numeric result isn't a real number (e.g. Complex[2.0, 3.0]).
 *
 * This is critical for Simplify[Sqrt[e^2] - e] when e contains numeric
 * sub-expressions like Cos[2] whose negativity defeats the structural
 * Plus-of-nonneg rule but whose contribution to e is small enough that
 * e is still numerically positive. The reason this is sound: if e is
 * numerically e.g. ~1.24 with 100-bit precision, no rearrangement of
 * the algebraic identity can make it secretly zero. */
static int numeric_sign_fallback(const Expr* x) {
    if (!x) return 0;
    /* Integer / bigint / real / mpfr literals already handled by the
     * inline numeric_sign() in prov_pos / prov_neg; skip to avoid
     * pointless allocation. */
    if (x->type == EXPR_INTEGER || x->type == EXPR_BIGINT ||
        x->type == EXPR_REAL) return 0;
    /* Performance gate: skip if x has any free symbol — numericalize
     * would just rebuild the same tree without producing a numeric
     * result, but the recursive walk is expensive for deep inputs. */
    if (nsf_has_free_symbol(x)) return 0;
#ifdef USE_MPFR
    if (x->type == EXPR_MPFR) return 0;
    NumericSpec spec;
    spec.mode = NUMERIC_MODE_MPFR;
    spec.bits = 100;  /* ~30 decimal digits */
#else
    NumericSpec spec = numeric_machine_spec();
#endif
    Expr* num = numericalize(x, spec);
    if (!num) return 0;
    int sign = 0;
#ifdef USE_MPFR
    if (num->type == EXPR_MPFR) {
        if (mpfr_number_p(num->data.mpfr)) {
            int s = mpfr_sgn(num->data.mpfr);
            if (s != 0) {
                /* |num| > 1e-20 ? */
                mpfr_t thresh;
                mpfr_init2(thresh, 64);
                mpfr_set_str(thresh, "1e-20", 10, MPFR_RNDN);
                if (mpfr_cmpabs(num->data.mpfr, thresh) > 0) {
                    sign = s;
                }
                mpfr_clear(thresh);
            }
        }
    } else
#endif
    if (num->type == EXPR_REAL) {
        double v = num->data.real;
        if (isfinite(v) && fabs(v) > 1e-10) {
            sign = (v > 0) ? 1 : -1;
        }
    } else if (num->type == EXPR_INTEGER) {
        if (num->data.integer > 0) sign = 1;
        else if (num->data.integer < 0) sign = -1;
    } else if (num->type == EXPR_BIGINT) {
        sign = mpz_sgn(num->data.bigint);
    }
    /* Any other shape (symbol, function, Complex, ...) leaves sign at 0. */
    expr_free(num);
    return sign;
}

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
bool prov_even(const AssumeCtx* ctx, const Expr* x) {
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

bool prov_pos(const AssumeCtx* ctx, const Expr* x) {
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
        /* Times: positive iff every factor positive. We DON'T early-return
         * `false` on a structural mismatch — falling through gives the
         * numerical-sign fallback below a chance to catch expressions that
         * are positive in spite of a locally-unprovable factor (e.g.
         * a product whose net sign is decidable numerically). */
        if (h == SYM_Times && n > 0) {
            bool all_pos = true;
            for (size_t i = 0; i < n; i++) {
                if (!prov_pos(ctx, a[i])) { all_pos = false; break; }
            }
            if (all_pos) return true;
        }
        /* Plus: at least one strictly positive, all others non-negative.
         * Same fall-through-on-failure as Times so the numeric fallback
         * gets a shot — needed for sums like 1 + Cos[2] + ... where one
         * term is locally negative but the total is unambiguously positive. */
        if (h == SYM_Plus && n > 0) {
            bool plus_ok = true;
            bool any = false;
            for (size_t i = 0; i < n; i++) {
                if (prov_pos(ctx, a[i])) { any = true; continue; }
                if (prov_nn(ctx, a[i])) continue;
                plus_ok = false; break;
            }
            if (plus_ok && any) return true;
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
    /* Last resort: numerical sign check. Handles the case where the
     * structural Plus-of-nonneg rule fails because some term (e.g. Cos[2])
     * is locally negative, but the overall expression is still positive
     * — see numeric_sign_fallback for the safety analysis. */
    if (numeric_sign_fallback(x) == 1) return true;
    return false;
}

bool prov_nn(const AssumeCtx* ctx, const Expr* x) {
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

bool prov_neg(const AssumeCtx* ctx, const Expr* x) {
    if (!x) return false;
    if (numeric_sign(x) == -1) return true;
    if (fact_directly_negative(ctx, x)) return true;
    /* Times: even number of negatives among factors, with the rest positive,
     * gives positive (not negative). For "negative" we need an odd number of
     * negative factors and the rest positive. v1 keeps this simple. */
    /* Last resort: numerical sign check — mirrors the fallback in
     * prov_pos. Catches expressions like 1 + Cos[2] + Cos[3] - 3 whose
     * negativity isn't visible structurally. */
    if (numeric_sign_fallback(x) == -1) return true;
    return false;
}

bool prov_np(const AssumeCtx* ctx, const Expr* x) {
    if (!x) return false;
    int s = numeric_sign(x);
    if (s == -1 || s == 0) return true;
    if (prov_neg(ctx, x)) return true;
    if (fact_directly_nonpos(ctx, x)) return true;
    return false;
}

bool prov_int(const AssumeCtx* ctx, const Expr* x) {
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

bool prov_re(const AssumeCtx* ctx, const Expr* x) {
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

