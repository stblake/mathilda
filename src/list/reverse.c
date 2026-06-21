#include "list_common.h"
#include "reverse.h"

static bool should_reverse_at_level(Expr* level_spec, size_t current_level) {
    if (!level_spec) return current_level == 1;
    if (level_spec->type == EXPR_INTEGER) return (size_t)level_spec->data.integer == current_level;
    if (level_spec->type == EXPR_FUNCTION && level_spec->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < level_spec->data.function.arg_count; i++) {
            if (level_spec->data.function.args[i]->type == EXPR_INTEGER && 
                (size_t)level_spec->data.function.args[i]->data.integer == current_level) return true;
        }
    }
    return false;
}

static Expr* reverse_rec(Expr* expr, Expr* level_spec, size_t current_level) {
    if (expr->type != EXPR_FUNCTION) return expr_copy(expr);

    size_t len = expr->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * len);
    bool do_rev = should_reverse_at_level(level_spec, current_level);

    for (size_t i = 0; i < len; i++) {
        size_t src_idx = do_rev ? (len - 1 - i) : i;
        new_args[i] = reverse_rec(expr->data.function.args[src_idx], level_spec, current_level + 1);
    }

    Expr* result = expr_new_function(expr_copy(expr->data.function.head), new_args, len);
    free(new_args);
    return result;
}

Expr* builtin_reverse(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* level_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    return reverse_rec(expr, level_spec, 1);
}
