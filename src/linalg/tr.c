/* Tr[list, f, n] -- generalised trace.
 *
 * Walks the diagonal of a rank-n nested list and combines the leaves with
 * the head `f` (default Plus).  For a 2-D matrix this is the usual sum of
 * diagonal entries.
 */

#include "linalg.h"
#include "eval.h"
#include "sym_names.h"
#include <stdlib.h>

static int64_t get_default_trace_depth(Expr* list) {
    int64_t depth = 0;
    Expr* curr = list;
    while (curr->type == EXPR_FUNCTION && curr->data.function.head->type == EXPR_SYMBOL && curr->data.function.head->data.symbol == SYM_List) {
        depth++;
        if (curr->data.function.arg_count == 0) break;
        curr = curr->data.function.args[0];
    }
    return depth;
}

static Expr* extract_diagonal_element(Expr* list, int64_t n, size_t index) {
    Expr* curr = list;
    for (int64_t level = 0; level < n; level++) {
        if (curr->type != EXPR_FUNCTION || curr->data.function.head->type != EXPR_SYMBOL || curr->data.function.head->data.symbol != SYM_List) {
            return NULL; // Not a list at this level
        }
        if (index >= curr->data.function.arg_count) {
            return NULL; // Index out of bounds
        }
        curr = curr->data.function.args[index];
    }
    return expr_copy(curr);
}

Expr* builtin_tr(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t count = res->data.function.arg_count;
    if (count < 1 || count > 3) return NULL;

    Expr* list = res->data.function.args[0];

    if (list->type != EXPR_FUNCTION || list->data.function.head->type != EXPR_SYMBOL || list->data.function.head->data.symbol != SYM_List) {
        return expr_copy(list);
    }

    Expr* f = NULL;
    bool free_f = false;
    if (count >= 2) {
        f = res->data.function.args[1];
    } else {
        f = expr_new_symbol("Plus");
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
