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
    if (e->data.function.head->data.symbol.name != SYM_Power) return false;
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
            Expr* mul = expr_new_function(expr_new_symbol(SYM_Times),
                                          mul_args, 2);
            prod_base = evaluate(mul);
            expr_free(mul);
            consumed[j] = true;
            group++;
        }
        consumed[i] = true;

        if (group >= 2) {
            Expr* pow_args[2] = { prod_base, expr_copy((Expr*)exp_i) };
            Expr* pow = expr_new_function(expr_new_symbol(SYM_Power),
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
    Expr* rebuilt = expr_new_function(expr_new_symbol(SYM_Times),
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
        && target->data.function.head->data.symbol.name == SYM_Times) {
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

Expr* simp_radicals(const Expr* e) {
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
Expr* transform_radical_canon(const Expr* e);

/* Forward declaration -- defined just below. */
bool is_sqrt(const Expr* e);

/* Rewrite every Power[Rational[m, n], 1/2] subtree to Times[Sqrt[m*n],
 * Rational[1, n]] (i.e. Sqrt[m*n]/n) for positive m, n. Returns a fresh
 * tree (caller owns); the input is not mutated.
 *
 * `transform_radical_canon` already attempts this rewrite, but the
 * Mathilda evaluator re-merges Times[Power[m, 1/2], Power[n, -1/2]]
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
                        expr_new_symbol(SYM_Power), sqrt_args, 2);
                    Expr* sqrt_e = evaluate(sqrt_call);
                    expr_free(sqrt_call);
                    Expr* prod_args[2] = {
                        sqrt_e,
                        make_rational(1, nn)
                    };
                    Expr* prod = expr_new_function(
                        expr_new_symbol(SYM_Times), prod_args, 2);
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
bool is_sqrt(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_Power) return false;
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
        && e->data.function.head->data.symbol.name == SYM_Power
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
                Expr* pc = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
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
        e->data.function.head->data.symbol.name != SYM_Times) return false;

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
            top->node->data.function.head->data.symbol.name != SYM_Times) {
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
        Expr* tn = expr_new_function(expr_new_symbol(SYM_Times),
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
        Expr* tn = expr_new_function(expr_new_symbol(SYM_Times), other, oi);
        b = evaluate(tn);
        expr_free(tn);
    }

    *out_b = b;
    *out_C = C;
    return true;
}

/* ----------------------------------------------------------------------- */
/* Phase-2 multi-extension denesting: budget, memo, structural helpers.    */
/*                                                                         */
/* Phase 1 (the original denester) handles Sqrt[A + b*Sqrt[C]] with A in Q */
/* and a single sqrt-bearing term. Phase 2 generalises both directions:    */
/*   (a) the splitter accepts radicands with multiple sqrt-bearing terms   */
/*       and iterates over candidates for the "outer" sqrt;                */
/*   (b) sqrt_if_clean_square recurses into the denester when its input    */
/*       discriminant D is itself a Q(Sqrt[k]) element of the same shape.  */
/* A small file-static depth budget bounds the recursion, and a memo keyed */
/* on the radicand suppresses repeated work across the bottom-up walker.   */
/* ----------------------------------------------------------------------- */

#define DENEST_MAX_DEPTH 4

typedef struct DenestMemoEntry {
    Expr* key;   /* radicand (the input to Sqrt[...]); owned */
    Expr* val;   /* denested Sqrt[key]; NULL = negative cache; owned */
} DenestMemoEntry;

typedef struct DenestMemo {
    DenestMemoEntry* entries;
    size_t count;
    size_t capacity;
} DenestMemo;

/* Single-threaded REPL: file-static state is safe. The top-level
 * simp_denest_sqrt entry seeds these and tears them down; re-entry from
 * an evaluator-driven sub-Simplify just shares the existing slot so we
 * never bypass the depth cap. */
static int g_denest_budget = 0;
static DenestMemo* g_denest_memo = NULL;

static bool denest_memo_lookup(const Expr* radicand,
                               Expr** out_val,
                               bool* out_negative) {
    *out_val = NULL;
    *out_negative = false;
    if (!g_denest_memo) return false;
    for (size_t i = 0; i < g_denest_memo->count; i++) {
        if (expr_eq(g_denest_memo->entries[i].key, (Expr*)radicand)) {
            if (g_denest_memo->entries[i].val) {
                *out_val = expr_copy(g_denest_memo->entries[i].val);
            } else {
                *out_negative = true;
            }
            return true;
        }
    }
    return false;
}

static void denest_memo_insert(const Expr* radicand, const Expr* val) {
    if (!g_denest_memo) return;
    if (g_denest_memo->count == g_denest_memo->capacity) {
        size_t newcap = g_denest_memo->capacity ? g_denest_memo->capacity * 2 : 8;
        DenestMemoEntry* na = (DenestMemoEntry*)realloc(
            g_denest_memo->entries, sizeof(DenestMemoEntry) * newcap);
        if (!na) return;  /* drop the insert on OOM */
        g_denest_memo->entries = na;
        g_denest_memo->capacity = newcap;
    }
    g_denest_memo->entries[g_denest_memo->count].key =
        expr_copy((Expr*)radicand);
    g_denest_memo->entries[g_denest_memo->count].val =
        val ? expr_copy((Expr*)val) : NULL;
    g_denest_memo->count++;
}

static void denest_memo_clear(DenestMemo* memo) {
    if (!memo) return;
    for (size_t i = 0; i < memo->count; i++) {
        expr_free(memo->entries[i].key);
        if (memo->entries[i].val) expr_free(memo->entries[i].val);
    }
    free(memo->entries);
    memo->entries = NULL;
    memo->count = 0;
    memo->capacity = 0;
}

/* Maximum number of stacked Sqrt heads along any root-to-leaf path.
 * Primary key in the outer-sqrt ranking heuristic: phase 2 wants the
 * discriminant to live one extension up rather than at the same level. */
static int sqrt_nesting_depth(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return 0;
    int best = 0;
    if (is_sqrt(e)) {
        int inner = sqrt_nesting_depth(e->data.function.args[0]);
        if (inner + 1 > best) best = inner + 1;
    }
    if (e->data.function.head) {
        int d = sqrt_nesting_depth(e->data.function.head);
        if (d > best) best = d;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int d = sqrt_nesting_depth(e->data.function.args[i]);
        if (d > best) best = d;
    }
    return best;
}

/* True when `plus_node` is a Plus[...] containing at least one
 * sqrt-bearing term (in the extract_sqrt_term sense). Gates the
 * phase-2 recursion in sqrt_if_clean_square. */
static bool plus_has_sqrt_bearing_term(const Expr* plus_node,
                                       const AssumeCtx* ctx) {
    if (!plus_node || plus_node->type != EXPR_FUNCTION) return false;
    if (!plus_node->data.function.head ||
        plus_node->data.function.head->type != EXPR_SYMBOL ||
        plus_node->data.function.head->data.symbol.name != SYM_Plus) return false;
    for (size_t i = 0; i < plus_node->data.function.arg_count; i++) {
        Expr* b = NULL; Expr* C = NULL;
        if (extract_sqrt_term(plus_node->data.function.args[i], ctx, &b, &C)) {
            expr_free(b); expr_free(C);
            return true;
        }
    }
    return false;
}

/* Structural check: does `e` contain any Sqrt[Plus[...]] subtree? The
 * phase-2 recursion accepts a candidate result only when this returns
 * false -- i.e., every nested sqrt was unrolled, leaving at worst
 * Sqrt[atom] leaves. */
static bool contains_sqrt_of_plus(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (is_sqrt(e)) {
        Expr* inner = e->data.function.args[0];
        if (inner && inner->type == EXPR_FUNCTION &&
            inner->data.function.head &&
            inner->data.function.head->type == EXPR_SYMBOL &&
            inner->data.function.head->data.symbol.name == SYM_Plus) {
            return true;
        }
    }
    if (e->data.function.head && contains_sqrt_of_plus(e->data.function.head))
        return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_sqrt_of_plus(e->data.function.args[i])) return true;
    }
    return false;
}

/* Forward declaration: the walker is defined further below near the
 * top-level entry; the phase-2 recursive helper needs it. */
static Expr* simp_denest_sqrt_walk(const Expr* e, const AssumeCtx* ctx);

/* Recursive denester entry used by sqrt_if_clean_square during phase-2.
 * Decrements the depth budget around a single walker call. Returns NULL
 * when the budget is exhausted or the walker found nothing to denest --
 * the caller specifically wants to know whether a rewrite fired, so
 * this does NOT fall back to expr_copy(). */
static Expr* simp_denest_sqrt_recursive(const Expr* e, const AssumeCtx* ctx) {
    if (g_denest_budget <= 1) return NULL;
    g_denest_budget--;
    Expr* r = simp_denest_sqrt_walk(e, ctx);
    g_denest_budget++;
    return r;
}

/* ----------------------------------------------------------------------- */

/* Given a clean square D, return a closed-form expression s with
 * s^2 = D (sign-agnostic; the caller validates the downstream
 * nonnegativity of (A+s)/2 and (A-s)/2 -- if both are provably nonneg,
 * the sign of s itself doesn't matter for the identity).
 *
 * Returns NULL when no closed form is detected. The cases handled are:
 *   - integer/bigint nonneg perfect square                 -> integer sqrt
 *   - rational with perfect-square num and den             -> rational sqrt
 *   - Power[u, 2k]                                          -> Power[u, k]
 *   - polynomial Plus whose Expand+FactorSquareFree is a
 *     pure even-power Power[u, 2k]                          -> Power[u, k]
 *   - phase-2: Plus carrying a sqrt-bearing term, where the
 *     recursive denester succeeds on Sqrt[D] with no surviving
 *     Sqrt[Plus] residue (memoised by D).
 *
 * The polynomial path covers the symbolic discriminants that case 6
 * ((x-y)^2) and case 7 (y^2) produce after expansion.
 */
static Expr* sqrt_if_clean_square(const Expr* D, const AssumeCtx* ctx) {

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
        Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), inv_args, 2);
        Expr* inv_e = evaluate(inv);
        expr_free(inv);
        Expr* prod_args[2] = { sn, inv_e };
        Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), prod_args, 2);
        Expr* out = evaluate(prod);
        expr_free(prod);
        return out;
    }

    /* Power[u, 2k] direct. */
    if (D->type == EXPR_FUNCTION
        && D->data.function.head
        && D->data.function.head->type == EXPR_SYMBOL
        && D->data.function.head->data.symbol.name == SYM_Power
        && D->data.function.arg_count == 2) {
        Expr* exp = D->data.function.args[1];
        if (exp->type == EXPR_INTEGER &&
            exp->data.integer >= 2 &&
            (exp->data.integer % 2) == 0) {
            int64_t k = exp->data.integer / 2;
            Expr* u = expr_copy(D->data.function.args[0]);
            if (k == 1) return u;
            Expr* pa[2] = { u, expr_new_integer(k) };
            Expr* pc = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
            Expr* out = evaluate(pc);
            expr_free(pc);
            return out;
        }
    }

    /* Polynomial Plus: try Expand + FactorSquareFree. The factoriser
     * stalls on inputs containing non-integer Power exponents (e.g. an
     * inner Sqrt that survived expansion), so gate on that. When the
     * gate trips, fall through to the phase-2 branch below rather than
     * returning -- a sqrt-bearing D is the phase-2 case by definition. */
    if (D->type == EXPR_FUNCTION
        && D->data.function.head
        && D->data.function.head->type == EXPR_SYMBOL
        && D->data.function.head->data.symbol.name == SYM_Plus
        && !has_non_integer_power(D)) {
        Expr* expanded = call_unary_copy("Expand", D);
        if (expanded) {
            if (has_non_integer_power(expanded)) {
                expr_free(expanded);
                /* fall through to phase 2 */
            } else {
                Expr* fsf = call_unary_owned("FactorSquareFree", expanded);
                if (fsf) {
                    if (fsf->type == EXPR_FUNCTION
                        && fsf->data.function.head
                        && fsf->data.function.head->type == EXPR_SYMBOL
                        && fsf->data.function.head->data.symbol.name == SYM_Power
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
                            Expr* pc = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
                            Expr* out = evaluate(pc);
                            expr_free(pc);
                            return out;
                        }
                    }
                    expr_free(fsf);
                }
            }
        }
    }

    /* Phase-2 multi-extension fallback: when D is itself a Plus carrying a
     * sqrt-bearing term (a Q(Sqrt[k]) element), recursively denest Sqrt[D].
     * Accept the result only if every nested Sqrt[Plus] was unrolled.
     * Memoised by D so the bottom-up walker and the discriminant path
     * don't re-denest the same radicand. */
    if (D->type == EXPR_FUNCTION &&
        D->data.function.head &&
        D->data.function.head->type == EXPR_SYMBOL &&
        D->data.function.head->data.symbol.name == SYM_Plus &&
        plus_has_sqrt_bearing_term(D, ctx)) {
        Expr* cached_val = NULL;
        bool cached_neg = false;
        if (denest_memo_lookup(D, &cached_val, &cached_neg)) {
            return cached_val;  /* NULL on negative cache, else owned copy */
        }

        Expr* sqrtD_args[2] = { expr_copy((Expr*)D), make_rational(1, 2) };
        Expr* sqrtD = expr_new_function(expr_new_symbol(SYM_Power),
                                         sqrtD_args, 2);
        Expr* sqrtD_e = evaluate(sqrtD);
        expr_free(sqrtD);

        Expr* denested = simp_denest_sqrt_recursive(sqrtD_e, ctx);
        expr_free(sqrtD_e);

        if (!denested || contains_sqrt_of_plus(denested)) {
            if (denested) expr_free(denested);
            denest_memo_insert(D, NULL);
            return NULL;
        }
        denest_memo_insert(D, denested);
        return denested;
    }

    return NULL;
}

/* Number of sqrt-bearing terms in a Plus expression. Phase-1 inputs
 * have exactly 1 (single-extension); phase-2 inputs have 2 or more.
 * Returns 0 for non-Plus inputs. */
static size_t count_outer_sqrt_candidates(const Expr* plus_node,
                                          const AssumeCtx* ctx) {
    if (!plus_node || plus_node->type != EXPR_FUNCTION) return 0;
    if (!plus_node->data.function.head ||
        plus_node->data.function.head->type != EXPR_SYMBOL ||
        plus_node->data.function.head->data.symbol.name != SYM_Plus) return 0;
    size_t k = 0;
    for (size_t i = 0; i < plus_node->data.function.arg_count; i++) {
        Expr* b = NULL; Expr* C = NULL;
        if (extract_sqrt_term(plus_node->data.function.args[i], ctx, &b, &C)) {
            expr_free(b); expr_free(C);
            k++;
        }
    }
    return k;
}

/* Given a Plus expression, partition its arguments into (A, b, C) such
 * that the Plus equals A + b * Sqrt[C]. `candidate_idx` selects which
 * sqrt-bearing term to treat as the "outer" sqrt: phase-1 single-sqrt
 * radicands have exactly one candidate (idx == 0), and phase-2
 * multi-sqrt radicands expose all of them so the caller can iterate.
 *
 * Candidates are ranked outer-first by sqrt nesting depth (deeper wins)
 * with leaf-complexity as the tie-break (smaller wins). For phase 2 this
 * lifts the discriminant A^2 - b^2*C one extension up, which is exactly
 * what the recursive sqrt_if_clean_square needs to fire. On single-sqrt
 * inputs the ranking is a no-op (one candidate at idx 0) and behaviour
 * is byte-identical to phase 1.
 *
 * Returns false when `candidate_idx` is out of range, the input isn't a
 * Plus, or no sqrt-bearing terms exist. *out_A, *out_b, *out_C are
 * caller-owned allocations on success. */
static bool split_plus_into_a_plus_b_sqrt_c(const Expr* plus_node,
                                            const AssumeCtx* ctx,
                                            size_t candidate_idx,
                                            Expr** out_A,
                                            Expr** out_b,
                                            Expr** out_C) {
    *out_A = NULL; *out_b = NULL; *out_C = NULL;
    if (!plus_node || plus_node->type != EXPR_FUNCTION) return false;
    if (!plus_node->data.function.head ||
        plus_node->data.function.head->type != EXPR_SYMBOL ||
        plus_node->data.function.head->data.symbol.name != SYM_Plus) return false;

    size_t n = plus_node->data.function.arg_count;
    if (n < 2) return false;

    /* Collect every sqrt-bearing term's index, nesting depth, and
     * complexity. Bound the candidate set; pathological multi-sqrt
     * inputs (>16 sqrt terms) degrade to phase-1 single-candidate
     * behaviour rather than blowing up the search. */
    enum { MAX_OUTER_CANDIDATES = 16 };
    struct OuterCand { size_t idx; int depth; size_t complexity; };
    struct OuterCand cands[MAX_OUTER_CANDIDATES];
    size_t k = 0;
    for (size_t i = 0; i < n && k < MAX_OUTER_CANDIDATES; i++) {
        Expr* bi = NULL; Expr* Ci = NULL;
        if (extract_sqrt_term(plus_node->data.function.args[i], ctx,
                              &bi, &Ci)) {
            cands[k].idx = i;
            cands[k].depth = sqrt_nesting_depth(Ci);
            cands[k].complexity = simp_default_complexity(Ci);
            expr_free(bi); expr_free(Ci);
            k++;
        }
    }
    if (k == 0) return false;
    if (candidate_idx >= k) return false;

    /* Selection sort: deeper first, then smaller complexity. k is
     * bounded by MAX_OUTER_CANDIDATES so this is O(k^2) <= O(256). */
    for (size_t i = 0; i + 1 < k; i++) {
        size_t best_j = i;
        for (size_t j = i + 1; j < k; j++) {
            if (cands[j].depth > cands[best_j].depth ||
                (cands[j].depth == cands[best_j].depth &&
                 cands[j].complexity < cands[best_j].complexity)) {
                best_j = j;
            }
        }
        if (best_j != i) {
            struct OuterCand tmp = cands[i];
            cands[i] = cands[best_j];
            cands[best_j] = tmp;
        }
    }

    size_t chosen_idx = cands[candidate_idx].idx;

    Expr* b = NULL; Expr* C = NULL;
    if (!extract_sqrt_term(plus_node->data.function.args[chosen_idx],
                           ctx, &b, &C)) {
        return false;
    }

    /* Build A as Plus of remaining args, evaluating to canonicalise. */
    if (n == 2) {
        *out_A = expr_copy(plus_node->data.function.args[1 - chosen_idx]);
    } else {
        Expr** a_args = (Expr**)malloc(sizeof(Expr*) * (n - 1));
        if (!a_args) { expr_free(b); expr_free(C); return false; }
        size_t ai = 0;
        for (size_t i = 0; i < n; i++) {
            if (i == chosen_idx) continue;
            a_args[ai++] = expr_copy(plus_node->data.function.args[i]);
        }
        Expr* pn = expr_new_function(expr_new_symbol(SYM_Plus), a_args, ai);
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
        const char* h = f->data.function.head->data.symbol.name;
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
        const char* h = f->data.function.head->data.symbol.name;
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
        x->data.function.head->data.symbol.name == SYM_Plus) {
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
                    arg->data.function.head->data.symbol.name == SYM_Times &&
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
                    const char* h = f->data.function.head->data.symbol.name;
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
     * forms Times[1/2, Plus[x, -y]] that Mathilda's evaluator produces
     * before Expand distributes them. */
    if (x && x->type == EXPR_FUNCTION && x->data.function.head &&
        x->data.function.head->type == EXPR_SYMBOL &&
        x->data.function.head->data.symbol.name == SYM_Times) {
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

/* True if `e` is a numeric rational literal (integer, bigint, or
 * Rational[n, d]). */
static bool is_numeric_rational(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    return is_rational_literal(e);
}

/* Sign of a numeric rational: -1, 0, or +1. Returns 0 for non-rationals
 * (caller is expected to gate via is_numeric_rational). */
static int rational_sign(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) {
        return (e->data.integer > 0) - (e->data.integer < 0);
    }
    if (e->type == EXPR_BIGINT) {
        return mpz_sgn(e->data.bigint);
    }
    if (is_rational_literal(e)) {
        const Expr* num = e->data.function.args[0];
        if (num->type == EXPR_INTEGER) {
            return (num->data.integer > 0) - (num->data.integer < 0);
        }
        if (num->type == EXPR_BIGINT) {
            return mpz_sgn(num->data.bigint);
        }
    }
    return 0;
}

/* Classify one Plus term as a Q(Sqrt[gamma]) contribution.
 *   - pure rational               => *alpha_inc_out (owned)
 *   - rational * Sqrt[rational]   => *beta_inc_out (owned), *gamma_out (owned)
 *   - anything else               => return false
 * gamma must be a nonneg rational; the caller treats *gamma_out == NULL
 * as "this term contributes only to alpha". */
static bool classify_q_sqrt_term(const Expr* term,
                                  Expr** alpha_inc_out,
                                  Expr** beta_inc_out,
                                  Expr** gamma_out) {
    *alpha_inc_out = NULL; *beta_inc_out = NULL; *gamma_out = NULL;
    if (!term) return false;
    if (is_numeric_rational(term)) {
        *alpha_inc_out = expr_copy((Expr*)term);
        return true;
    }
    if (is_sqrt(term)) {
        Expr* g = term->data.function.args[0];
        if (!is_numeric_rational(g) || rational_sign(g) < 0) return false;
        *beta_inc_out = expr_new_integer(1);
        *gamma_out = expr_copy(g);
        return true;
    }
    if (term->type == EXPR_FUNCTION && term->data.function.head &&
        term->data.function.head->type == EXPR_SYMBOL &&
        term->data.function.head->data.symbol.name == SYM_Times) {
        size_t n = term->data.function.arg_count;
        if (n == 0) return false;
        Expr** coef_args = (Expr**)malloc(sizeof(Expr*) * n);
        if (!coef_args) return false;
        size_t ci = 0;
        const Expr* sqrt_factor = NULL;
        for (size_t i = 0; i < n; i++) {
            const Expr* a = term->data.function.args[i];
            if (is_numeric_rational(a)) {
                coef_args[ci++] = expr_copy((Expr*)a);
            } else if (!sqrt_factor && is_sqrt(a)) {
                Expr* g = a->data.function.args[0];
                if (!is_numeric_rational(g) || rational_sign(g) < 0) {
                    for (size_t j = 0; j < ci; j++) expr_free(coef_args[j]);
                    free(coef_args);
                    return false;
                }
                sqrt_factor = a;
            } else {
                for (size_t j = 0; j < ci; j++) expr_free(coef_args[j]);
                free(coef_args);
                return false;
            }
        }
        if (!sqrt_factor) {
            /* Pure rational product. */
            if (ci == 0) { free(coef_args); return false; }
            if (ci == 1) {
                *alpha_inc_out = coef_args[0];
                free(coef_args);
            } else {
                Expr* tn = expr_new_function(expr_new_symbol(SYM_Times),
                                              coef_args, ci);
                *alpha_inc_out = evaluate(tn);
                expr_free(tn);
            }
            return true;
        }
        Expr* coef;
        if (ci == 0) {
            free(coef_args);
            coef = expr_new_integer(1);
        } else if (ci == 1) {
            coef = coef_args[0];
            free(coef_args);
        } else {
            Expr* tn = expr_new_function(expr_new_symbol(SYM_Times),
                                          coef_args, ci);
            coef = evaluate(tn);
            expr_free(tn);
        }
        *beta_inc_out = coef;
        *gamma_out = expr_copy(sqrt_factor->data.function.args[0]);
        return true;
    }
    return false;
}

/* Build a fresh rational expression representing a^2 - b^2 * g, with
 * a, b, g consumed by the call. */
static Expr* compute_discriminant_eval(Expr* a, Expr* b, Expr* g) {
    Expr* a2_args[2] = { a, expr_new_integer(2) };
    Expr* a2_p = expr_new_function(expr_new_symbol(SYM_Power), a2_args, 2);
    Expr* a2 = evaluate(a2_p);
    expr_free(a2_p);

    Expr* b2_args[2] = { b, expr_new_integer(2) };
    Expr* b2_p = expr_new_function(expr_new_symbol(SYM_Power), b2_args, 2);
    Expr* b2 = evaluate(b2_p);
    expr_free(b2_p);

    Expr* bg_args[2] = { b2, g };
    Expr* bg_t = expr_new_function(expr_new_symbol(SYM_Times), bg_args, 2);
    Expr* bg = evaluate(bg_t);
    expr_free(bg_t);

    Expr* neg_args[2] = { expr_new_integer(-1), bg };
    Expr* neg_t = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
    Expr* neg_bg = evaluate(neg_t);
    expr_free(neg_t);

    Expr* D_args[2] = { a2, neg_bg };
    Expr* D_p = expr_new_function(expr_new_symbol(SYM_Plus), D_args, 2);
    Expr* D = evaluate(D_p);
    expr_free(D_p);
    return D;
}

/* If `e` is a Q(Sqrt[gamma]) element with gamma a single nonneg rational
 * and every Plus term consistent with the same gamma, decide
 * nonnegativity by exact sign analysis on alpha vs beta and on the
 * discriminant alpha^2 - beta^2 * gamma. Returns true iff e is
 * provably >= 0. No numeric sampling -- all comparisons run on exact
 * rationals via the evaluator.
 *
 * This is what unlocks the phase-2 branch check on Q values like
 * 11 - 2*Sqrt[29], where 11^2 = 121 > 116 = (2)^2 * 29 proves Q > 0. */
static bool q_sqrt_extension_is_nonneg(const Expr* e) {
    if (!e) return false;

    if (is_numeric_rational(e)) {
        return numeric_is_nonneg(e);
    }

    /* Single-term forms (Sqrt[r] or r*Sqrt[r] or pure-rational Times). */
    if (is_sqrt(e) ||
        (e->type == EXPR_FUNCTION && e->data.function.head &&
         e->data.function.head->type == EXPR_SYMBOL &&
         e->data.function.head->data.symbol.name == SYM_Times)) {
        Expr* a = NULL; Expr* b = NULL; Expr* g = NULL;
        if (!classify_q_sqrt_term(e, &a, &b, &g)) return false;
        bool ok;
        if (a) {
            ok = numeric_is_nonneg(a);
            expr_free(a);
        } else {
            ok = numeric_is_nonneg(b);
            expr_free(b); expr_free(g);
        }
        return ok;
    }

    if (e->type != EXPR_FUNCTION ||
        !e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_Plus) {
        return false;
    }

    /* Plus form: collect alpha parts, beta parts, and a single shared
     * gamma. Any term that doesn't classify cleanly aborts. */
    size_t n = e->data.function.arg_count;
    Expr** alpha_parts = (Expr**)malloc(sizeof(Expr*) * (n + 1));
    Expr** beta_parts  = (Expr**)malloc(sizeof(Expr*) * (n + 1));
    if (!alpha_parts || !beta_parts) {
        free(alpha_parts); free(beta_parts);
        return false;
    }
    size_t ai = 0, bi = 0;
    Expr* common_gamma = NULL;
    bool ok = true;
    for (size_t i = 0; i < n; i++) {
        Expr* a = NULL; Expr* b = NULL; Expr* g = NULL;
        if (!classify_q_sqrt_term(e->data.function.args[i], &a, &b, &g)) {
            ok = false; break;
        }
        if (a) alpha_parts[ai++] = a;
        if (b) {
            if (common_gamma) {
                if (!expr_eq(common_gamma, g)) {
                    ok = false;
                    expr_free(b); expr_free(g);
                    break;
                }
                expr_free(g);
            } else {
                common_gamma = g;
            }
            beta_parts[bi++] = b;
        }
    }
    if (!ok) {
        for (size_t j = 0; j < ai; j++) expr_free(alpha_parts[j]);
        for (size_t j = 0; j < bi; j++) expr_free(beta_parts[j]);
        free(alpha_parts); free(beta_parts);
        if (common_gamma) expr_free(common_gamma);
        return false;
    }

    Expr* alpha;
    if (ai == 0) {
        alpha = expr_new_integer(0);
        free(alpha_parts);
    } else if (ai == 1) {
        alpha = alpha_parts[0];
        free(alpha_parts);
    } else {
        Expr* pn = expr_new_function(expr_new_symbol(SYM_Plus), alpha_parts, ai);
        alpha = evaluate(pn);
        expr_free(pn);
    }
    Expr* beta;
    if (bi == 0) {
        beta = expr_new_integer(0);
        free(beta_parts);
    } else if (bi == 1) {
        beta = beta_parts[0];
        free(beta_parts);
    } else {
        Expr* pn = expr_new_function(expr_new_symbol(SYM_Plus), beta_parts, bi);
        beta = evaluate(pn);
        expr_free(pn);
    }

    if (!is_numeric_rational(alpha) || !is_numeric_rational(beta)) {
        expr_free(alpha); expr_free(beta);
        if (common_gamma) expr_free(common_gamma);
        return false;
    }

    int sa = rational_sign(alpha);
    int sb = rational_sign(beta);

    if (bi == 0 || sb == 0) {
        if (common_gamma) expr_free(common_gamma);
        expr_free(beta);
        bool r = (sa >= 0);
        expr_free(alpha);
        return r;
    }

    if (sa >= 0 && sb >= 0) {
        expr_free(alpha); expr_free(beta); expr_free(common_gamma);
        return true;
    }
    if (sa <= 0 && sb <= 0) {
        bool r = (sa == 0 && sb == 0);
        expr_free(alpha); expr_free(beta); expr_free(common_gamma);
        return r;
    }

    /* Mixed signs: compare alpha^2 to beta^2 * gamma exactly. */
    Expr* D = compute_discriminant_eval(expr_copy(alpha), expr_copy(beta),
                                         expr_copy(common_gamma));
    bool r;
    if (sa > 0 && sb < 0) {
        /* e = alpha - |beta|*Sqrt[gamma]; nonneg iff D = a^2 - b^2 g >= 0. */
        r = is_numeric_rational(D) && numeric_is_nonneg(D);
    } else { /* sa < 0 && sb > 0 */
        /* e = -|alpha| + beta*Sqrt[gamma]; nonneg iff D <= 0. */
        r = is_numeric_rational(D) && rational_sign(D) <= 0;
    }
    expr_free(alpha); expr_free(beta); expr_free(common_gamma); expr_free(D);
    return r;
}

/* Decide whether `e` is provably nonneg under the active context.
 * Numeric literals get a free direct check; Q(Sqrt[k]) elements get
 * exact sign analysis; everything else routes through assume_known_nonneg
 * augmented by a local transitive prover (denest_prov_nonneg) that
 * chains inequality facts. The transitive step is required for case 7's
 * branch check: x > y && y > 0 implies (x+y)/2 >= 0, but
 * assume_known_nonneg alone can't see it. */
static bool denest_is_nonneg(const Expr* e, const AssumeCtx* ctx) {
    if (numeric_is_nonneg(e)) return true;
    if (q_sqrt_extension_is_nonneg(e)) return true;
    return denest_prov_nonneg(ctx, e, 4);
}

/* Compute (P, Q, s) for a specific choice of "outer" sqrt in `plus_node`.
 * `candidate_idx` selects which sqrt-bearing term to use (see
 * split_plus_into_a_plus_b_sqrt_c). Phase 1 single-sqrt radicands hit
 * idx == 0 once. Phase 2 multi-sqrt radicands are explored by
 * denest_compute_pq_s iterating over candidates until one succeeds.
 *
 * Returns false when the chosen split doesn't yield a clean denesting
 * (no closed-form sqrt(A^2 - b^2 C), or branch validity fails). The
 * outer wrapper retries the next candidate on false. */
static bool denest_compute_pq_s_at_candidate(const Expr* plus_node,
                                              const AssumeCtx* ctx,
                                              size_t candidate_idx,
                                              Expr** P_out, Expr** Q_out,
                                              Expr** s_out,
                                              bool* b_is_negative_out) {
    *P_out = NULL; *Q_out = NULL; *s_out = NULL;
    Expr* A = NULL; Expr* b = NULL; Expr* C = NULL;
    if (!split_plus_into_a_plus_b_sqrt_c(plus_node, ctx, candidate_idx,
                                          &A, &b, &C)) return false;

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
    Expr* b2_pow = expr_new_function(expr_new_symbol(SYM_Power), b2_args, 2);
    Expr* b2 = evaluate(b2_pow);
    expr_free(b2_pow);
    Expr* B_args[2] = { b2, expr_copy(C) };
    Expr* B_times = expr_new_function(expr_new_symbol(SYM_Times), B_args, 2);
    Expr* B = evaluate(B_times);
    expr_free(B_times);

    /* D = A^2 - B; Expand to put it in canonical polynomial form so
     * sqrt_if_clean_square's FactorSquareFree path has a chance. */
    Expr* A2_args[2] = { expr_copy(A), expr_new_integer(2) };
    Expr* A2_pow = expr_new_function(expr_new_symbol(SYM_Power), A2_args, 2);
    Expr* A2 = evaluate(A2_pow);
    expr_free(A2_pow);
    Expr* negB_args[2] = { expr_new_integer(-1), B };
    Expr* negB_times = expr_new_function(expr_new_symbol(SYM_Times), negB_args, 2);
    Expr* negB = evaluate(negB_times);
    expr_free(negB_times);
    Expr* D_args[2] = { A2, negB };
    Expr* D_plus = expr_new_function(expr_new_symbol(SYM_Plus), D_args, 2);
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

    /* Compute P = (A + s)/2 and Q = (A - s)/2. Mathilda's evaluator
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
    Expr* Ap = expr_new_function(expr_new_symbol(SYM_Plus), Ap_args, 2);
    Expr* Ap_e = evaluate(Ap);
    expr_free(Ap);
    Expr* halfP_args[2] = { make_rational(1, 2), Ap_e };
    Expr* halfP = expr_new_function(expr_new_symbol(SYM_Times), halfP_args, 2);
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
    Expr* neg_s = expr_new_function(expr_new_symbol(SYM_Times), neg_s_args, 2);
    Expr* neg_s_e = evaluate(neg_s);
    expr_free(neg_s);
    Expr* Am_args[2] = { expr_copy(A), neg_s_e };
    Expr* Am = expr_new_function(expr_new_symbol(SYM_Plus), Am_args, 2);
    Expr* Am_e = evaluate(Am);
    expr_free(Am);
    Expr* halfQ_args[2] = { make_rational(1, 2), Am_e };
    Expr* halfQ = expr_new_function(expr_new_symbol(SYM_Times), halfQ_args, 2);
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

/* Outer iterator for the half-sum denesting identity. Phase 1 inputs
 * (single sqrt-bearing term) hit candidate 0 once with no overhead;
 * phase 2 inputs explore every sqrt-bearing term ranked outer-first by
 * sqrt_nesting_depth so the discriminant lives in Q(Sqrt[k]) for one
 * k, letting sqrt_if_clean_square's recursive case fire. */
static bool denest_compute_pq_s(const Expr* plus_node,
                                 const AssumeCtx* ctx,
                                 Expr** P_out, Expr** Q_out, Expr** s_out,
                                 bool* b_is_negative_out) {
    *P_out = NULL; *Q_out = NULL; *s_out = NULL;
    size_t ncands = count_outer_sqrt_candidates(plus_node, ctx);
    if (ncands == 0) return false;
    for (size_t i = 0; i < ncands; i++) {
        if (denest_compute_pq_s_at_candidate(plus_node, ctx, i,
                                              P_out, Q_out, s_out,
                                              b_is_negative_out)) {
            return true;
        }
        /* On failure, *_out are already cleaned up by the inner
         * function's failure paths -- reset before the next candidate. */
        *P_out = NULL; *Q_out = NULL; *s_out = NULL;
    }
    return false;
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
        e->data.function.head->data.symbol.name != SYM_Power) return false;
    if (e->data.function.arg_count != 2) return false;
    if (!e->data.function.args[0] || !e->data.function.args[1]) return false;
    Expr* exp = e->data.function.args[1];
    if (!is_rational_literal(exp)) return false;
    /* is_rational_literal only guarantees arg_count == 2; the two
     * children may still be NULL during transient build states. */
    if (!exp->data.function.args[0] || !exp->data.function.args[1]) return false;
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
    Expr* pow = expr_new_function(expr_new_symbol(SYM_Power), args, 2);
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
    Expr* sP = expr_new_function(expr_new_symbol(SYM_Power), sP_args, 2);
    Expr* sP_e = evaluate(sP);
    expr_free(sP);
    Expr* sQ_args[2] = { Q, make_rational(1, 2) };
    Expr* sQ = expr_new_function(expr_new_symbol(SYM_Power), sQ_args, 2);
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
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), sum_args, 2);
        R = evaluate(sum);
        expr_free(sum);
    } else {
        Expr* neg_sQ_args[2] = { expr_new_integer(-1), sQ_e };
        Expr* neg_sQ = expr_new_function(expr_new_symbol(SYM_Times),
                                          neg_sQ_args, 2);
        Expr* neg_sQ_e = evaluate(neg_sQ);
        expr_free(neg_sQ);
        Expr* sum_args[2] = { sP_e, neg_sQ_e };
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), sum_args, 2);
        R = evaluate(sum);
        expr_free(sum);
    }

    int64_t k = m < 0 ? -m : m;
    Expr* R_pow = (k == 1) ? R : power_int_eval(R, k);

    Expr* result;
    if (m < 0) {
        /* Divide by s^k. Build R_pow * Power[s, -k] and evaluate. */
        Expr* inv_s_args[2] = { s, expr_new_integer(-k) };
        Expr* inv_s = expr_new_function(expr_new_symbol(SYM_Power),
                                         inv_s_args, 2);
        Expr* inv_s_e = evaluate(inv_s);
        expr_free(inv_s);
        Expr* prod_args[2] = { R_pow, inv_s_e };
        Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
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
                /* OOM: drop the in-flight rewrite and propagate failure
                 * by treating this node as "no change". Returning NULL
                 * (no fired denesting) is correct -- the caller will
                 * fall back to the unmodified subtree. */
                if (!new_args) {
                    expr_free(r);
                    return NULL;
                }
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
            && base->data.function.head->data.symbol.name == SYM_Plus) {
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
 * NULL. When no denesting fires, returns expr_copy(e).
 *
 * Manages the phase-2 depth budget and memoisation cache. Re-entry
 * detection: if g_denest_budget > 0 we're being called from inside an
 * in-flight denesting (e.g. evaluator-driven sub-Simplify) -- share
 * the existing budget and memo rather than resetting them, so the
 * recursion cap is enforced end-to-end. */
Expr* simp_denest_sqrt(const Expr* e, const AssumeCtx* ctx) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    bool top_level = (g_denest_budget == 0);
    DenestMemo memo;
    if (top_level) {
        memo.entries = NULL;
        memo.count = 0;
        memo.capacity = 0;
        g_denest_budget = DENEST_MAX_DEPTH;
        g_denest_memo = &memo;
    }

    Expr* r = simp_denest_sqrt_walk(e, ctx);
    Expr* out = r ? r : expr_copy((Expr*)e);

    if (top_level) {
        denest_memo_clear(&memo);
        g_denest_memo = NULL;
        g_denest_budget = 0;
    }

    if (dbg) simp_debug_log("DenestSqrt", e, out,
                            simp_debug_elapsed_ms(t0));
    return out;
}

