#include "eval.h"
#include "list_common.h"
#include "split.h"
#include "linkedlist_seq.h"

Expr* builtin_split(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* list = res->data.function.args[0];
    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    
    /* A LinkedList splits into a LinkedList of LinkedList runs, read from its
     * slots (a snapshot: the test may change the list). */
    LLView lv;
    if (!ll_view(list, &lv)) return builtin_normal_error("Split", res, 1);
    
    size_t count = lv.n;
    if (count == 0) {
        Expr* r = ll_view_linked(&lv) ? ll_view_result(&lv, NULL, 0) : expr_copy(list);
        ll_view_done(&lv);
        return r;
    }
    
    Expr** result_runs = malloc(sizeof(Expr*) * count);
    size_t num_runs = 0;
    
    size_t run_start = 0;
    for (size_t i = 1; i <= count; i++) {
        bool split_here = (i == count);
        if (!split_here) {
            Expr* prev = lv.a[i-1];
            Expr* curr = lv.a[i];
            bool identical = false;
            if (test == NULL) {
                identical = expr_eq(prev, curr);
            } else {
                Expr* call_args[2] = { expr_copy(prev), expr_copy(curr) };
                Expr* call = expr_new_function(expr_copy(test), call_args, 2);
                Expr* eval_res = eval_callback(call);
                if (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol.name == SYM_True) {
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
                run_args[j] = expr_copy(lv.a[run_start + j]);
            }
            result_runs[num_runs++] = ll_view_result(&lv, run_args, run_len);
            free(run_args);
            run_start = i;
        }
    }
    
    Expr* result = ll_view_result(&lv, result_runs, num_runs);
    free(result_runs);
    ll_view_done(&lv);
    return result;
}
