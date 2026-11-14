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
 *
 * Gather[list, test] takes the equivalence from a two-argument test: an
 * element joins the first sublist whose FIRST element e gives test[e, elem] ===
 * True, else it opens a new sublist.  No hash key exists for an arbitrary test,
 * so this is O(n g) test calls for g groups (each element is compared with one
 * representative per group, never with every member).
 *
 *     Gather[{1, 2, 3, 4}, Mod[#1, 2] == Mod[#2, 2] &]  ->  {{1, 3}, {2, 4}}
 * -------------------------------------------------------------------------- */

#include "eval.h"
#include "list_common.h"
#include "gather.h"
#include "assoc.h"

/* Gather[list, test]: see the file comment. */
static Expr* gather_with_test(Expr* res) {
    Expr* list = res->data.function.args[0];
    Expr* test = res->data.function.args[1];
    if (list->type != EXPR_FUNCTION) return builtin_normal_error("Gather", res, 1);
    size_t n = list->data.function.arg_count;
    /* group g holds members[g][0 .. len[g]-1], indices into list */
    size_t* group_of = malloc(sizeof(size_t) * (n ? n : 1));
    size_t* rep = malloc(sizeof(size_t) * (n ? n : 1));     /* first index per group */
    size_t ng = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* elem = list->data.function.args[i];
        size_t g = 0;
        for (; g < ng; g++) {
            Expr* call = expr_new_function(expr_copy(test),
                (Expr*[]){ expr_copy(list->data.function.args[rep[g]]), expr_copy(elem) }, 2);
            Expr* r = eval_callback(call);
            bool same = r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True;
            if (r) expr_free(r);
            expr_free(call);
            if (same) break;
        }
        if (g == ng) rep[ng++] = i;
        group_of[i] = g;
    }
    size_t* len = calloc(ng ? ng : 1, sizeof(size_t));
    for (size_t i = 0; i < n; i++) len[group_of[i]]++;
    Expr** groups = malloc(sizeof(Expr*) * (ng ? ng : 1));
    size_t* fill = calloc(ng ? ng : 1, sizeof(size_t));
    Expr*** members = malloc(sizeof(Expr**) * (ng ? ng : 1));
    for (size_t g = 0; g < ng; g++) members[g] = malloc(sizeof(Expr*) * len[g]);
    for (size_t i = 0; i < n; i++) {
        size_t g = group_of[i];
        members[g][fill[g]++] = expr_copy(list->data.function.args[i]);
    }
    for (size_t g = 0; g < ng; g++) {
        groups[g] = expr_new_function(expr_copy(list->data.function.head), members[g], len[g]);
        free(members[g]);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), groups, ng);
    free(members); free(fill); free(groups); free(len); free(rep); free(group_of);
    return out;
}

Expr* builtin_gather(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count == 2) return gather_with_test(res);
    if (res->data.function.arg_count != 1) return NULL;
    /* NULL key function == identity: the element is its own group key. */
    return assoc_gather_core(res->data.function.args[0], NULL);
}
