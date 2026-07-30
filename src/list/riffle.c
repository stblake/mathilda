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
    if (list->type != EXPR_FUNCTION) return NULL;

    size_t n = list->data.function.arg_count;

    /* Normalise the separator to a borrowed vector: a List cycles through its
     * elements, any other expression is one separator used in every gap. */
    Expr** seps;
    size_t k;
    if (is_listq(sep)) {
        seps = sep->data.function.args;
        k = sep->data.function.arg_count;
    } else {
        seps = &sep;
        k = 1;
    }

    /* No gaps (n <= 1), or nothing to put in them (k == 0): copy through.
     * Must precede the 2n - 1 sizing below — see the header comment. */
    if (n <= 1 || k == 0) return expr_copy(list);

    size_t out_count = 2 * n - 1;
    Expr** out = malloc(out_count * sizeof(Expr*));
    if (!out) return NULL;

    for (size_t i = 0; i < n; i++) {
        out[2 * i] = expr_copy(list->data.function.args[i]);
        /* The gap following element i is gap i + 1 in 1-based terms, so it
         * takes separator index i mod k. */
        if (i + 1 < n) out[2 * i + 1] = expr_copy(seps[i % k]);
    }

    Expr* result = expr_new_function(expr_copy(list->data.function.head), out, out_count);
    free(out);
    return result;
}
