#include "list_common.h"
#include "join.h"

/* ------------------- Join ------------------- */

/*
 * join_at_level: recursively join expressions at a given depth.
 * level 1 = concatenate the top-level arguments (basic join).
 * level 2 = descend one level into each list, joining corresponding
 *           sub-elements across the input lists. Ragged arrays are
 *           handled by concatenating successive elements at that level.
 *
 * lists: array of n_lists expressions to join
 * n_lists: number of lists
 * level: remaining depth to descend (>= 1)
 *
 * Returns a new Expr* on success, or NULL if types don't match.
 */
static Expr* join_at_level(Expr** lists, size_t n_lists, int level) {
    if (n_lists == 0) return NULL;

    Expr* first = lists[0];
    if (first->type != EXPR_FUNCTION) return NULL;

    Expr* head = first->data.function.head;

    /* Verify all lists share the same head */
    for (size_t i = 1; i < n_lists; i++) {
        if (lists[i]->type != EXPR_FUNCTION ||
            !expr_eq(lists[i]->data.function.head, head))
            return NULL;
    }

    if (level == 1) {
        /* Base case: concatenate all arguments */
        size_t total = 0;
        for (size_t i = 0; i < n_lists; i++)
            total += lists[i]->data.function.arg_count;

        Expr** new_args = malloc(sizeof(Expr*) * (total > 0 ? total : 1));
        size_t curr = 0;
        for (size_t i = 0; i < n_lists; i++) {
            Expr* li = lists[i];
            for (size_t j = 0; j < li->data.function.arg_count; j++)
                new_args[curr++] = expr_copy(li->data.function.args[j]);
        }
        Expr* result = expr_new_function(expr_copy(head), new_args, total);
        free(new_args);
        return result;
    }

    /* Recursive case: level > 1.
     * Find the maximum number of sub-elements across all lists.
     * For each position k, gather the sub-elements from each list
     * that have a k-th element and recursively join them at level-1.
     * Lists that don't have a k-th element simply contribute nothing. */
    size_t max_len = 0;
    for (size_t i = 0; i < n_lists; i++) {
        if (lists[i]->data.function.arg_count > max_len)
            max_len = lists[i]->data.function.arg_count;
    }

    Expr** result_args = malloc(sizeof(Expr*) * (max_len > 0 ? max_len : 1));
    size_t result_count = 0;

    /* Temporary buffer for gathering sub-elements at position k */
    Expr** subs = malloc(sizeof(Expr*) * n_lists);

    for (size_t k = 0; k < max_len; k++) {
        size_t n_subs = 0;
        for (size_t i = 0; i < n_lists; i++) {
            if (k < lists[i]->data.function.arg_count)
                subs[n_subs++] = lists[i]->data.function.args[k];
        }

        if (n_subs == 0) continue;

        if (n_subs == 1) {
            /* Only one list contributes at this position: copy as-is */
            result_args[result_count++] = expr_copy(subs[0]);
        } else {
            Expr* joined = join_at_level(subs, n_subs, level - 1);
            if (!joined) {
                /* Cleanup on failure */
                for (size_t j = 0; j < result_count; j++)
                    expr_free(result_args[j]);
                free(result_args);
                free(subs);
                return NULL;
            }
            result_args[result_count++] = joined;
        }
    }

    free(subs);
    Expr* result = expr_new_function(expr_copy(head), result_args, result_count);
    free(result_args);
    return result;
}

Expr* builtin_join(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1)
        return NULL;

    size_t n_args = res->data.function.arg_count;
    int level = 1;
    size_t n_lists = n_args;

    /* Check if the last argument is an integer (level specification) */
    Expr* last = res->data.function.args[n_args - 1];
    if (n_args >= 2 && last->type == EXPR_INTEGER) {
        level = (int)last->data.integer;
        if (level < 1) return NULL;
        n_lists = n_args - 1;
    }

    if (n_lists < 1) return NULL;

    Expr* result = join_at_level(res->data.function.args, n_lists, level);
    if (!result) return NULL;

    return result;
}
