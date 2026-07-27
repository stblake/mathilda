/* ---------------------------------------------------------------------------
 * gather.c — Gather[list], the identity case of GatherBy.
 *
 * Gather partitions a list into sublists of structurally identical elements:
 * every element appears in exactly one sublist, two elements share a sublist
 * iff expr_eq holds between them, sublists appear in order of the first
 * occurrence of their element, and within a sublist elements keep their input
 * order. Unlike Split, grouping is not restricted to adjacent runs:
 *
 *     Gather[{1, 7, 3, 7, 2, 3, 9}]  ->  {{1}, {7, 7}, {3, 3}, {2}, {9}}
 *     Gather[{a, b, a}]              ->  {{a, a}, {b}}
 *
 * The grouping itself is not reimplemented here. assoc_gather_core (assoc.c)
 * is the same hash-indexed O(n) engine that backs GatherBy; passing a NULL key
 * function selects the identity key, so Gather[l] and GatherBy[l, Identity]
 * agree by construction, and the identity path skips the n Identity[x]
 * applications that spelling it as GatherBy[l, Identity] would evaluate.
 * -------------------------------------------------------------------------- */

#include "list_common.h"
#include "gather.h"
#include "assoc.h"

Expr* builtin_gather(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    /* NULL key function == identity: the element is its own group key. */
    return assoc_gather_core(res->data.function.args[0], NULL);
}
