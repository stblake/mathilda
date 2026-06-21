#include "list_common.h"
#include "split.h"

Expr* builtin_split(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* list = res->data.function.args[0];
    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    
    if (list->type != EXPR_FUNCTION) return expr_copy(list);
    
    size_t count = list->data.function.arg_count;
    if (count == 0) return expr_copy(list);
    
    Expr** result_runs = malloc(sizeof(Expr*) * count);
    size_t num_runs = 0;
    
    size_t run_start = 0;
    for (size_t i = 1; i <= count; i++) {
        bool split_here = (i == count);
        if (!split_here) {
            Expr* prev = list->data.function.args[i-1];
            Expr* curr = list->data.function.args[i];
            bool identical = false;
            if (test == NULL) {
                identical = expr_eq(prev, curr);
            } else {
                Expr* call_args[2] = { expr_copy(prev), expr_copy(curr) };
                Expr* call = expr_new_function(expr_copy(test), call_args, 2);
                Expr* eval_res = evaluate(call);
                if (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol == SYM_True) {
                    identical = true;
                }
                expr_free(eval_res);
                expr_free(call);
            }
            if (!identical) split_here = true;
        }
        
        if (split_here) {
            size_t run_len = i - run_start;
            Expr** run_args = malloc(sizeof(Expr*) * run_len);
            for (size_t j = 0; j < run_len; j++) {
                run_args[j] = expr_copy(list->data.function.args[run_start + j]);
            }
            result_runs[num_runs++] = expr_new_function(expr_copy(list->data.function.head), run_args, run_len);
            free(run_args);
            run_start = i;
        }
    }
    
    Expr* result = expr_new_function(expr_copy(list->data.function.head), result_runs, num_runs);
    free(result_runs);
    return result;
}
