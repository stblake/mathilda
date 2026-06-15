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
Expr* apply_assumption_rules(const Expr* input, const AssumeCtx* ctx) {
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
         * Mathilda's matcher does NOT perform orderless-Plus subset matching
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
                              expr_new_function(expr_new_symbol(SYM_Times),
                                  (Expr*[]){ expr_new_integer(-1), expr_copy(rhs) }, 2) };
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), sub_args, 2);
        Expr* diff = eval_and_free(sum);
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
        all[fill++] = expr_new_function(expr_new_symbol(SYM_Rule), direct, 2);

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
                    other_args[oi++] = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1),
                                   expr_copy(diff->data.function.args[k]) }, 2);
                }
                Expr* iso_rhs;
                if (n - 1 == 1) {
                    iso_rhs = other_args[0];
                    free(other_args);
                } else {
                    iso_rhs = expr_new_function(expr_new_symbol(SYM_Plus), other_args, n - 1);
                    free(other_args);
                }
                Expr* iso[2] = { expr_copy(term), iso_rhs };
                all[fill++] = expr_new_function(expr_new_symbol(SYM_Rule), iso, 2);
            }
        }
    }
    for (size_t i = 0; i < ctx->count; i++) if (eq_diffs[i]) expr_free(eq_diffs[i]);
    free(eq_diffs);
    if (string_rules) expr_free(string_rules);

    Expr* rules_list = expr_new_function(expr_new_symbol(SYM_List), all, fill);
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

