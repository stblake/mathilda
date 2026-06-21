#include "list_common.h"
#include "flatten.h"

static void flatten_rec(Expr* e, const char* h, int64_t level, Expr*** results, size_t* count, size_t* cap) {
    if (level != 0 && head_is(e, intern_symbol(h))) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            flatten_rec(e->data.function.args[i], h, level == -1 ? -1 : level - 1, results, count, cap);
        }
    } else {
        if (*count == *cap) {
            *cap *= 2;
            *results = realloc(*results, sizeof(Expr*) * (*cap));
        }
        (*results)[(*count)++] = expr_copy(e);
    }
}

Expr* builtin_flatten(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 3) return NULL;

    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION) return expr_copy(list);

    int64_t n = -1; // -1 means infinity
    if (res->data.function.arg_count >= 2) {
        if (res->data.function.args[1]->type != EXPR_INTEGER) return NULL;
        n = res->data.function.args[1]->data.integer;
    }

    const char* h = SYM_List;
    if (res->data.function.arg_count == 3) {
        if (res->data.function.args[2]->type != EXPR_SYMBOL) return NULL;
        h = res->data.function.args[2]->data.symbol;
    }

    size_t cap = 16;
    size_t count = 0;
    Expr** results = malloc(sizeof(Expr*) * cap);

    // Initial call: we flatten children of the head if they also have head h.
    for (size_t i = 0; i < list->data.function.arg_count; i++) {
        flatten_rec(list->data.function.args[i], h, n, &results, &count, &cap);
    }

    Expr* result = expr_new_function(expr_copy(list->data.function.head), results, count);
    free(results);
    return result;
}
