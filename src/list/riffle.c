/* Riffle — interleave separators into the gaps of a list.
 *
 * Mathematica semantics:
 *
 *   Riffle[list, x]                 x is placed in every gap:
 *                                   Riffle[{1,2,3}, 0] -> {1, 0, 2, 0, 3}
 *   Riffle[list, {x1, ..., xk}]     the xi are consumed in order and cycle
 *                                   back to x1 after xk, filling the gaps
 *                                   left to right:
 *                                   Riffle[{a,b,c,d}, {x,y}] ->
 *                                     {a, x, b, y, c, x, d}
 *
 * THE GAP INVARIANT
 *
 * Separators go only BETWEEN consecutive elements — never before the first and
 * never after the last. A list of n elements therefore has exactly n - 1 gaps,
 * and the output has 2n - 1 slots. Two consequences drive the code below:
 *
 *   - n <= 1 means there are no gaps at all, so the result is the input
 *     unchanged whatever the separator is. This case is checked BEFORE the
 *     2n - 1 output sizing, because with n == 0 that expression underflows
 *     size_t to SIZE_MAX and the allocation would be nonsense.
 *   - separators past the last gap are never indexed, so
 *     Riffle[{a,b,c}, {x,y,z}] -> {a, x, b, y, c} simply never reaches z.
 *
 * An empty separator list has nothing to interleave, so it also passes the
 * list through unchanged; that check doubles as the guard that keeps the
 * cycling index from dividing by zero.
 *
 * The head of the first argument is preserved rather than forced to List, so
 * Riffle[f[a,b], x] gives f[a, x, b]. That is also what makes Riffle[{}, 0]
 * come back as {} with no special case.
 *
 * PERFORMANCE
 *
 * One pass, O(n) element copies, and a single exactly-sized allocation — the
 * output length is known up front from n, so no growable buffer is needed.
 *
 * NOT HANDLED: packed arrays (EXPR_NDARRAY) are a distinct representation from
 * List and are left unevaluated here; see md-2aa. */

#include "list_common.h"
#include "ndarray.h"    /* is_ndarray */
#include "ndstruct.h"   /* ndstruct_delist_repack — packed-argument fallback */
#include "riffle.h"
#include "linkedlist_seq.h" /* LLSeq: LinkedList list or separators */

Expr* builtin_riffle(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* list = res->data.function.args[0];
    Expr* sep = res->data.function.args[1];

    /* Native buffer path first; ndstruct_delist_repack below stays the fallback
     * for every form it declines. See ndstruct.h. */
    if (is_ndarray(list)) {
        Expr* fast = ndstruct_riffle(res);
        if (fast) return fast;
        return ndstruct_delist_repack(res, list);
    }

    /* Anything without arguments to interleave — an atom — stays unevaluated. */
    /* A LinkedList riffles from its slots into a new LinkedList. */
    LLSeq lv;
    if (!llq_open(list, &lv)) return NULL;

    size_t n = lv.n;

    /* Normalise the separator to a borrowed vector: a List (or a LinkedList)
     * cycles through its elements, any other expression is one separator used
     * in every gap. */
    LLSeq sv;                 /* the separators, read on demand */
    bool seq_sep = is_listq(sep) || sep->type == EXPR_LINKEDLIST;
    if (seq_sep) llq_open(sep, &sv);
    else         llq_open(NULL, &sv);
    size_t k = seq_sep ? sv.n : 1;

    /* No gaps (n <= 1), or nothing to put in them (k == 0): copy through.
     * Must precede the 2n - 1 sizing below — see the header comment. */
    if (n <= 1 || k == 0) {
        Expr* r;
        if (llq_linked(&lv)) {
            r = llq_range(&lv, 0, n);          /* a new list: the slots copied */
        } else {
            r = expr_copy(list);
        }
        llq_close(&sv);
        llq_close(&lv);
        return r;
    }

    /* The gap following element i is gap i + 1 in 1-based terms, so it takes
     * separator index i mod k. A LinkedList's elements (and a LinkedList of
     * separators) are copied slot to slot into the new list. */
    Expr* result;
    if (llq_linked(&lv)) {
        LinkedList* ol = ll_new();
        for (size_t i = 0; i < n; i++) {
            llq_push(ol, &lv, i);
            if (i + 1 < n) {
                if (seq_sep) llq_push(ol, &sv, i % k);
                else         ll_push_back(ol, expr_copy(sep));
            }
        }
        result = expr_new_linkedlist(ol);
    } else {
        size_t out_count = 2 * n - 1;
        Expr** out = malloc(out_count * sizeof(Expr*));
        if (!out) { llq_close(&sv); llq_close(&lv); return NULL; }
        for (size_t i = 0; i < n; i++) {
            out[2 * i] = expr_copy(llq_peek(&lv, i));
            if (i + 1 < n) out[2 * i + 1] = expr_copy(seq_sep ? llq_peek(&sv, i % k) : sep);
        }
        result = llq_result(&lv, out, out_count);
        free(out);
    }
    llq_close(&sv);
    llq_close(&lv);
    return result;
}
