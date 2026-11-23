/* quartiles.c -- Quartiles[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "arithmetic.h"
#include "sym_names.h"
#include "ndreduce.h"
#include "sort.h"
#include "linkedlist.h" /* LLReader */
#include "message.h"
#include "print.h"     /* expr_to_string */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

Expr* builtin_quartiles(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* data = res->data.function.args[0];
    Expr* param_expr = NULL;
    if (res->data.function.arg_count == 2) {
        param_expr = res->data.function.args[1];
    }

    /* NDArray fast path: sort the flat buffer, default quartile method (ndreduce.c). */
    if (ndred_call_has_ndarray(res)) return ndred_quartiles(res);

    /* LinkedList: the packed kernel on the chunks, or columnwise (stats_ll.c). */
    { Expr* r; if (stats_ll_arm(res, STATS_QUARTILES, &r)) return r; }

    if (!ll_is_seq(data)) return expr_copy(res);
    if (ll_seq_length(data) == 0) return expr_copy(res);

    if (data->type != EXPR_LINKEDLIST &&
        data->data.function.args[0]->type == EXPR_FUNCTION &&
        data->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
        data->data.function.args[0]->data.function.head->data.symbol.name == SYM_List) {

        /* Ragged / mixed depth: Transpose stays unevaluated, and mapping over
         * it re-entered Quartiles with the same data until $RecursionLimit. */
        if (!stats_is_matrix(data)) {
            stats_msg_rectn("Quartiles", res);
            return expr_copy(res);
        }
        Expr* t_expr = expr_new_function(expr_new_symbol(SYM_Transpose), (Expr*[]){expr_copy(data)}, 1);
        Expr* t_eval = evaluate(t_expr);
        expr_free(t_expr);

        if (t_eval->type != EXPR_FUNCTION || t_eval->data.function.head->type != EXPR_SYMBOL ||
            t_eval->data.function.head->data.symbol.name != SYM_List) {
            expr_free(t_eval);
            return expr_copy(res);
        }

        size_t cols = t_eval->data.function.arg_count;
        Expr** q_args = malloc(sizeof(Expr*) * cols);
        for (size_t i = 0; i < cols; i++) {
            Expr* call_args[2];
            call_args[0] = expr_copy(t_eval->data.function.args[i]);
            size_t call_cnt = 1;
            if (param_expr) {
                call_args[1] = expr_copy(param_expr);
                call_cnt = 2;
            }
            Expr* q_call = expr_new_function(expr_new_symbol(SYM_Quartiles), call_args, call_cnt);
            q_args[i] = evaluate(q_call);
            expr_free(q_call);
        }
        Expr* res_list = expr_new_function(expr_new_symbol(SYM_List), q_args, cols);
        free(q_args);
        expr_free(t_eval);
        return res_list;
    }

    LLReader e;
    ll_rd_open(data, &e);
    size_t n = e.n;

    bool all_real = true;
    for (size_t i = 0; i < n; i++) {
        if (!stats_is_real_numeric(ll_rd_at(&e, i))) {
            all_real = false;
            break;
        }
    }

    if (!all_real) {
        ll_rd_close(&e);
        char* str = expr_to_string(res);
        mth_message("Quartiles", "rectn", "Rectangular array of real numbers is expected at position 1 in %s.", str);
        free(str);
        return expr_copy(res);
    }

    Expr *a = NULL, *b = NULL, *c = NULL, *d = NULL;
    if (param_expr) {
        if (param_expr->type == EXPR_FUNCTION && param_expr->data.function.arg_count == 2) {
            Expr* row1 = param_expr->data.function.args[0];
            Expr* row2 = param_expr->data.function.args[1];
            if (row1->type == EXPR_FUNCTION && row1->data.function.arg_count == 2 &&
                row2->type == EXPR_FUNCTION && row2->data.function.arg_count == 2) {
                a = expr_copy(row1->data.function.args[0]);
                b = expr_copy(row1->data.function.args[1]);
                c = expr_copy(row2->data.function.args[0]);
                d = expr_copy(row2->data.function.args[1]);
            }
        }
    }

    if (!a || !b || !c || !d) {
        if (a) expr_free(a);
        if (b) expr_free(b);
        if (c) expr_free(c);
        if (d) expr_free(d);
        a = make_rational(1, 2);
        b = expr_new_integer(0);
        c = expr_new_integer(0);
        d = expr_new_integer(1);
    }

    /* Numeric order (stats_sort_order: by value, not Sort's canonical order),
     * as a permutation over the elements where they lie. */
    size_t* perm = stats_sort_order(&e);
    if (!perm) {
        ll_rd_close(&e);
        expr_free(a); expr_free(b); expr_free(c); expr_free(d);
        return expr_copy(res);
    }
    StatsOrder ord = { &e, perm };
    Expr* q_vals[3];
    q_vals[0] = make_rational(1, 4);
    q_vals[1] = make_rational(1, 2);
    q_vals[2] = make_rational(3, 4);

    Expr* results[3];
    /* Shared engine (stats_common.c) -- extracted from the loop that lived here;
     * see stats_quantile_point's header comment for the one semantic difference
     * (integer h), which Quartiles' c=0 default cannot observe. */
    for (int k = 0; k < 3; k++) {
        results[k] = stats_quantile_point(&ord, n, q_vals[k], a, b, c, d);
    }

    expr_free(q_vals[0]);
    expr_free(q_vals[1]);
    expr_free(q_vals[2]);

    free(perm);
    ll_rd_close(&e);
    expr_free(a);
    expr_free(b);
    expr_free(c);
    expr_free(d);

    return expr_new_function(expr_new_symbol(SYM_List), results, 3);
}
