/* Tr[list, f, n] -- generalised trace.
 *
 * Walks the diagonal of a rank-n nested list and combines the leaves with
 * the head `f` (default Plus).  For a 2-D matrix this is the usual sum of
 * diagonal entries.
 */

#include "linalg.h"
#include "linkedlist.h"
#include "ndlinalg.h"
#include "eval.h"
#include "sym_names.h"
#include <stdlib.h>

/* The levels of a List or LinkedList (a LinkedList is the List of its
 * elements) nested along the first-element spine. */
static int64_t get_default_trace_depth(Expr* list) {
    int64_t depth = 0;
    const Expr* curr = list;
    bool own = false;
    Expr* held = NULL;                    /* a boxed slot read from a typed list */
    while (ll_is_seq(curr)) {
        depth++;
        int64_t n = curr->type == EXPR_LINKEDLIST ? ll_length(curr->data.llist)
                                                  : (int64_t)curr->data.function.arg_count;
        if (n == 0) break;
        bool o;
        Expr* next = ll_seq_view(curr, 0, &o);
        if (held) ll_elem_done(held, own);
        held = next; own = o;
        curr = next;
    }
    if (held) ll_elem_done(held, own);
    return depth;
}

/* Element {index, index, ..., index} (n levels) of a nest of Lists /
 * LinkedLists, as a new reference; NULL when a level is not a container or
 * is too short. */
static Expr* extract_diagonal_element(Expr* list, int64_t n, size_t index) {
    Expr* curr = expr_copy(list);
    for (int64_t level = 0; level < n; level++) {
        if (!ll_is_seq(curr)) { expr_free(curr); return NULL; }
        int64_t len = curr->type == EXPR_LINKEDLIST ? ll_length(curr->data.llist)
                                                    : (int64_t)curr->data.function.arg_count;
        if ((int64_t)index >= len) { expr_free(curr); return NULL; }
        Expr* next = curr->type == EXPR_LINKEDLIST ? ll_get(curr->data.llist, (int64_t)index)
                                                   : expr_copy(curr->data.function.args[index]);
        expr_free(curr);
        curr = next;
    }
    return curr;
}

Expr* builtin_tr(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (linalg_call_has_ndarray(res)) return ndla_tr(res);
    size_t count = res->data.function.arg_count;
    if (count < 1 || count > 3) return NULL;

    Expr* list = res->data.function.args[0];

    if (!ll_is_seq(list)) return expr_copy(list);

    Expr* f = NULL;
    bool free_f = false;
    if (count >= 2) {
        f = res->data.function.args[1];
    } else {
        f = expr_new_symbol(SYM_Plus);
        free_f = true;
    }

    int64_t n = 0;
    if (count == 3) {
        Expr* n_expr = res->data.function.args[2];
        if (n_expr->type == EXPR_INTEGER) {
            n = n_expr->data.integer;
            if (n < 0) n = 0;
        } else {
            if (free_f) expr_free(f);
            return NULL;
        }
    } else {
        n = get_default_trace_depth(list);
        if (n == 0) n = 1;
    }

    if (n == 0) {
        if (free_f) expr_free(f);
        return expr_copy(list);
    }

    size_t cap = 16;
    size_t num_elems = 0;
    Expr** elements = malloc(sizeof(Expr*) * cap);

    size_t i = 0;
    while (true) {
        Expr* elem = extract_diagonal_element(list, n, i);
        if (!elem) {
            break;
        }
        if (num_elems == cap) {
            cap *= 2;
            elements = realloc(elements, sizeof(Expr*) * cap);
        }
        elements[num_elems++] = elem;
        i++;
    }

    Expr* result;
    if (num_elems == 0) {
        result = eval_and_free(expr_new_function(expr_copy(f), NULL, 0));
    } else {
        result = eval_and_free(expr_new_function(expr_copy(f), elements, num_elems));
    }

    free(elements);
    if (free_f) expr_free(f);

    return result;
}
