#include "list_common.h"
#include "rotate.h"

static Expr* rotate_rec(Expr* expr, Expr* n_spec, size_t level_idx) {
    if (expr->type != EXPR_FUNCTION) return expr_copy(expr);

    int64_t n = 0;
    if (n_spec->type == EXPR_INTEGER) {
        if (level_idx == 0) n = n_spec->data.integer;
    } else if (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->data.symbol.name == SYM_List) {
        if (level_idx < n_spec->data.function.arg_count) {
            Expr* sub_n = n_spec->data.function.args[level_idx];
            if (sub_n->type == EXPR_INTEGER) n = sub_n->data.integer;
        }
    }

    size_t len = expr->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * len);

    if (len > 0) {
        int64_t offset = n % (int64_t)len;
        if (offset < 0) offset += (int64_t)len;

        for (size_t i = 0; i < len; i++) {
            size_t old_idx = (i + (size_t)offset) % len;
            new_args[i] = rotate_rec(expr->data.function.args[old_idx], n_spec, level_idx + 1);
        }
    }

    Expr* result = expr_new_function(expr_copy(expr->data.function.head), new_args, len);
    free(new_args);
    return result;
}

Expr* builtin_rotateleft(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* n_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    
    Expr* default_n = NULL;
    if (!n_spec) {
        default_n = expr_new_integer(1);
        n_spec = default_n;
    }

    Expr* ret = rotate_rec(expr, n_spec, 0);
    if (default_n) expr_free(default_n);
    return ret;
}

Expr* builtin_rotateright(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* n_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    
    Expr* neg_n_spec = NULL;
    if (!n_spec) {
        neg_n_spec = expr_new_integer(-1);
    } else if (n_spec->type == EXPR_INTEGER) {
        neg_n_spec = expr_new_integer(-n_spec->data.integer);
    } else if (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->data.symbol.name == SYM_List) {
        Expr** neg_args = malloc(sizeof(Expr*) * n_spec->data.function.arg_count);
        for (size_t i = 0; i < n_spec->data.function.arg_count; i++) {
            if (n_spec->data.function.args[i]->type == EXPR_INTEGER) {
                neg_args[i] = expr_new_integer(-n_spec->data.function.args[i]->data.integer);
            } else {
                neg_args[i] = expr_copy(n_spec->data.function.args[i]);
            }
        }
        neg_n_spec = expr_new_function(expr_new_symbol(SYM_List), neg_args, n_spec->data.function.arg_count);
        free(neg_args);
    } else {
        return NULL;
    }

    Expr* ret = rotate_rec(expr, neg_n_spec, 0);
    expr_free(neg_n_spec);
    return ret;
}
