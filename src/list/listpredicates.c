#include "list_common.h"
#include "listpredicates.h"

Expr* builtin_listq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (is_listq(arg)) {
        return expr_new_symbol(SYM_True);
    }
    return expr_new_symbol(SYM_False);
}

Expr* builtin_vectorq(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 1 && res->data.function.arg_count != 2)) return NULL;
    Expr* arg = res->data.function.args[0];
    if (!is_listq(arg)) return expr_new_symbol(SYM_False);

    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    for (size_t i = 0; i < arg->data.function.arg_count; i++) {
        Expr* elem = arg->data.function.args[i];
        if (test == NULL) {
            if (is_listq(elem)) return expr_new_symbol(SYM_False);
        } else {
            Expr* call_args[1] = { expr_copy(elem) };
            Expr* call = expr_new_function(expr_copy(test), call_args, 1);
            Expr* eval_res = evaluate(call);
            bool is_true = (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol == SYM_True);
            expr_free(eval_res);
            expr_free(call);
            if (!is_true) return expr_new_symbol(SYM_False);
        }
    }
    return expr_new_symbol(SYM_True);
}

Expr* builtin_matrixq(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 1 && res->data.function.arg_count != 2)) return NULL;
    Expr* arg = res->data.function.args[0];
    if (!is_listq(arg)) return expr_new_symbol(SYM_False);

    if (arg->data.function.arg_count == 0) return expr_new_symbol(SYM_False);

    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    size_t col_count = 0;
    bool first_row = true;

    for (size_t i = 0; i < arg->data.function.arg_count; i++) {
        Expr* row = arg->data.function.args[i];
        if (!is_listq(row)) return expr_new_symbol(SYM_False);

        if (first_row) {
            col_count = row->data.function.arg_count;
            first_row = false;
        } else {
            if (row->data.function.arg_count != col_count) return expr_new_symbol(SYM_False);
        }

        for (size_t j = 0; j < row->data.function.arg_count; j++) {
            Expr* elem = row->data.function.args[j];
            if (test == NULL) {
                if (is_listq(elem)) return expr_new_symbol(SYM_False);
            } else {
                Expr* call_args[1] = { expr_copy(elem) };
                Expr* call = expr_new_function(expr_copy(test), call_args, 1);
                Expr* eval_res = evaluate(call);
                bool is_true = (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol == SYM_True);
                expr_free(eval_res);
                expr_free(call);
                if (!is_true) return expr_new_symbol(SYM_False);
            }
        }
    }

    return expr_new_symbol(SYM_True);
}
