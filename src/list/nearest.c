/* Nearest[list, x] -- the element(s) of `list` closest to the target `x`.
 *
 *   Nearest[{1, 5, 10}, 3]    -> {1, 5}
 *   Nearest[{10, 20, 30}, 100] -> {30}
 *   Nearest[{}, 3]             -> {}
 *
 * Distance is Abs[element - x], composed from the existing internal_subtract
 * and internal_abs rather than a bespoke helper, so a complex element uses its
 * modulus for free: Nearest[{3 + 4 I, 1}, 0] is {1} because the distances are
 * 5 and 1.
 *
 * ALL elements tied at the minimum distance are returned, in their original
 * order. That requirement is what fixes the algorithm. The obvious neighbour,
 * RankedMin (sort.c), selects with a quickselect whose comparator carries an
 * original-index tiebreak (sort.c:932) existing precisely to make ties
 * IMPOSSIBLE, so that exactly one element can win -- the opposite of what is
 * needed here. The shape used instead is MinimalBy's (sort.c:663-716): one pass
 * to find the minimum, a second to collect every distance equal to it. Input
 * order among ties then falls out of the ascending collect pass, so there is no
 * tie logic in this file at all.
 *
 * Nearest diverges from MinimalBy in exactly one respect, and deliberately.
 * Every distance must be a real number, or the whole call stays unevaluated.
 * MinimalBy[{1, a, 3}, Abs[# - 2] &] answers {1, 3}: expr_compare orders
 * symbols after all numbers (sort.c:379-380), so the symbolic element is never
 * minimal and silently vanishes from a result that still looks plausible.
 * Gating on the DISTANCE rather than on the element covers a symbolic element,
 * a symbolic target, and a non-real complex in a single check. A symbolic real
 * such as Pi is rejected too -- Abs[Pi - 3] stays as Abs[-3 + Pi] -- rather
 * than being numericalized the way RankedMin's ranked_numeric_key would.
 *
 * Cost: O(n) distance evaluations, O(n) comparisons, O(n) peak extra memory.
 * The two evaluate passes per element dominate, so this is an interpreter-speed
 * path and not a buffer path. There is no exact-hit short circuit: a later
 * element can tie at distance 0, and dropping it would break the tie contract.
 *
 * Only the two-argument form lives here. The n-nearest, radius, rule,
 * all-pairs, and NearestTo operator forms, and the DistanceFunction option, are
 * separate follow-ups. */

#include "list_common.h"
#include "internal.h"
#include "nearest.h"

/* Abs[e - x], fully evaluated. Caller owns the result.
 *
 * internal_call_impl (internal.c:189-211) consumes the argument array, and when
 * the builtin declines it returns the UNEVALUATED node rather than NULL -- so a
 * symbolic operand comes back as Abs[...], which the caller's numeric gate then
 * rejects. Neither call can yield NULL, so neither needs a null check. The
 * composition and the eval_and_free wrapping follow comparisons.c:313-314,
 * which already builds a difference this way for the Equal zero test. */
static Expr* nearest_distance(Expr* e, Expr* x) {
    Expr* sub_args[2] = { expr_copy(e), expr_copy(x) };
    Expr* diff = eval_and_free(internal_subtract(sub_args, 2));
    Expr* abs_args[1] = { diff };          /* internal_abs takes ownership */
    return eval_and_free(internal_abs(abs_args, 1));
}

/* The numeric gate and the value comparator used below --
 * list_real_number_q / list_numeric_cmp -- live in list_common.{c,h}. They
 * were written here for Nearest and were promoted when FindClusters became a
 * second caller; the block comment in list_common.h records why neither
 * expr_compare nor expr_numeric_sign can serve in their place, and why the
 * mixed exact/inexact path must not go through a subtraction. */

Expr* builtin_nearest(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* list = res->data.function.args[0];
    Expr* x    = res->data.function.args[1];

    /* A PACKED list has already been materialised on the way in, because
     * Nearest is not on pack.c's AWARE list. A VISIBLE NDArray is not, and is
     * not a List either, so it lands here and stays unevaluated rather than
     * being silently truncated. Following RankedMin (sort.c:1014) rather than
     * MinimalBy, which accepts and preserves any head. */
    if (!is_listq(list)) return NULL;

    size_t n = list->data.function.arg_count;
    Expr** elem = list->data.function.args;

    /* Empty in, empty out -- checked before the gate, so Nearest[{}, a] is {}
     * and not unevaluated. Matches MaximalBy (sort.c:684). */
    if (n == 0) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    Expr** dist = malloc(sizeof(Expr*) * n);
    if (!dist) return NULL;

    /* One pass to build the distances. A distance that is not a real number
     * means there is no definite answer: free what we built and decline. */
    for (size_t i = 0; i < n; i++) {
        dist[i] = nearest_distance(elem[i], x);
        if (!list_real_number_q(dist[i])) {
            for (size_t j = 0; j <= i; j++) expr_free(dist[j]);
            free(dist);
            return NULL;
        }
    }

    /* Two passes over list_numeric_cmp: find the minimum, then collect every
     * distance equal to it. An undecidable comparison declines the whole call
     * rather than guessing a tie. */
    bool ok = true;
    size_t best = 0;
    for (size_t i = 1; i < n && ok; i++)
        if (list_numeric_cmp(dist[i], dist[best], &ok) < 0 && ok) best = i;

    Expr** out = ok ? malloc(sizeof(Expr*) * n) : NULL;
    if (!out) {
        for (size_t i = 0; i < n; i++) expr_free(dist[i]);
        free(dist);
        return NULL;
    }

    /* Ascending index order is what preserves input order among ties. */
    size_t nout = 0;
    for (size_t i = 0; i < n && ok; i++)
        if (list_numeric_cmp(dist[i], dist[best], &ok) == 0 && ok)
            out[nout++] = expr_copy(elem[i]);

    for (size_t i = 0; i < n; i++) expr_free(dist[i]);
    free(dist);

    if (!ok) {
        for (size_t i = 0; i < nout; i++) expr_free(out[i]);
        free(out);
        return NULL;
    }

    /* The wrapper is always List. Stated explicitly because the input head is
     * necessarily List here, but the n-nearest follow-up keeps the same rule. */
    Expr* head = expr_new_symbol(SYM_List);
    Expr* result = expr_new_function(head, out, nout);
    if (!result) {                       /* OOM: expr_new_function took nothing */
        for (size_t i = 0; i < nout; i++) expr_free(out[i]);
        expr_free(head);
    }
    free(out);
    return result;
}
