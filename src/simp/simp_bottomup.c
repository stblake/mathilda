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
/* Bottom-up Simplify: memoised recursive descent over subtrees            */
/* ----------------------------------------------------------------------- */

#define SIMP_BOTTOMUP_MAX_DEPTH 64
/* SIMP_MEMO_BUCKETS, SimpMemoEntry, and SimpMemo all live in
 * simp_internal.h since the table is instantiated in simp_builtins.c. */

void simp_memo_init(SimpMemo* m) {
    for (int i = 0; i < SIMP_MEMO_BUCKETS; i++) m->buckets[i] = NULL;
}

void simp_memo_free(SimpMemo* m) {
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
const Expr* simp_memo_get(SimpMemo* m, const Expr* key) {
    uint64_t h = expr_hash(key) % SIMP_MEMO_BUCKETS;
    for (SimpMemoEntry* e = m->buckets[h]; e; e = e->next) {
        if (expr_eq(e->key, key)) return e->value;
    }
    return NULL;
}

/* Stores deep copies of both key and value. */
void simp_memo_put(SimpMemo* m, const Expr* key, const Expr* value) {
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
/* simp_bottomup -- recursive driver                                       */
/* ----------------------------------------------------------------------- */

Expr* simp_bottomup(const Expr* input, const AssumeCtx* ctx,
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

