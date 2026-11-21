#include "eval.h"
#include "list_common.h"
#include "split.h"
#include "linkedlist_seq.h"
#include <math.h>

/* Adjacent slots i-1, i of a typed buffer equal as expr_eq compares the
 * machine numbers they hold (NaN equals NaN; a complex pair both parts). */
static bool split_raw_eq(LLKind k, const void* buf, size_t i) {
    if (k == LL_I64) return ((const int64_t*)buf)[i - 1] == ((const int64_t*)buf)[i];
    const double* d = buf;
    size_t w = (k == LL_CPLX) ? 2 : 1;
    for (size_t c = 0; c < w; c++) {
        double x = d[(i - 1) * w + c], y = d[i * w + c];
        if (!((isnan(x) && isnan(y)) || x == y)) return false;
    }
    return true;
}

Expr* builtin_split(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* list = res->data.function.args[0];
    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    
    /* A LinkedList splits into a LinkedList of LinkedList runs: its elements
     * read from a snapshot one at a time (the test may change the list), or a
     * typed list's slots compared raw, and each run copied span-wise. */
    LLSeq q;
    if (!llq_open(list, &q)) return builtin_normal_error("Split", res, 1);
    
    size_t count = q.n;
    if (count == 0) {
        Expr* r = llq_linked(&q) ? llq_emit(&q, NULL, 0) : expr_copy(list);
        llq_close(&q);
        return r;
    }
    LLKind kind = llq_kind(&q);
    void* raw = NULL;
    if (!test && llq_linked(&q) && kind != LL_BOXED) {
        raw = malloc(ll_slot_size(q.snap) * count);
        ll_gather(q.snap, 0, (int64_t)count, raw);
    }
    
    Expr** result_runs = malloc(sizeof(Expr*) * count);
    size_t num_runs = 0;
    
    size_t run_start = 0;
    for (size_t i = 1; i <= count; i++) {
        bool split_here = (i == count);
        if (!split_here) {
            bool identical = false;
            if (raw) {
                identical = split_raw_eq(kind, raw, i);
            } else {
                Expr* prev = expr_copy(llq_peek(&q, i - 1));
                Expr* curr = expr_copy(llq_peek(&q, i));
                if (test == NULL) {
                    identical = expr_eq(prev, curr);
                    expr_free(prev); expr_free(curr);
                } else {
                    Expr* call_args[2] = { prev, curr };
                    Expr* call = expr_new_function(expr_copy(test), call_args, 2);
                    Expr* eval_res = eval_callback(call);
                    if (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol.name == SYM_True) {
                        identical = true;
                    }
                    expr_free(eval_res);
                    expr_free(call);
                }
            }
            if (!identical) split_here = true;
        }
        
        if (split_here) {
            result_runs[num_runs++] = llq_range(&q, run_start, i);
            run_start = i;
        }
    }
    
    Expr* result = llq_result(&q, result_runs, num_runs);
    free(result_runs);
    free(raw);
    llq_close(&q);
    return result;
}
