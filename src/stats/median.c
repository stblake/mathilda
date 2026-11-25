/* median.c -- Median[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "sym_names.h"
#include "assoc.h"
#include "ndreduce.h"
#include "sort.h"
#include "linkedlist.h" /* LLReader */
#include "message.h"
#include "print.h"     /* expr_to_string */
#include <stdio.h>
#include <stdlib.h>
#include "dist.h"    /* ml_dist_property: distribution moments */
#include "common.h"  /* head_is */

Expr* builtin_median(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* Median[assoc] is the median of the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* Median of a distribution (NormalDistribution, UniformDistribution, a 1-D
     * SmoothKernel): its closed form, see src/ml/dist.h. */
    if (data->type == EXPR_FUNCTION && !head_is(data, SYM_List)) {
        Expr* r = ml_dist_property(ML_PROP_MEDIAN, data, NULL);
        if (r) return r;
    }

    /* NDArray fast path: quickselect on the flat buffer (see ndreduce.c). */
    if (ndred_call_has_ndarray(res)) return ndred_median(res);

    /* LinkedList: the packed kernel on the chunks, or columnwise (stats_ll.c). */
    { Expr* r; if (stats_ll_arm(res, STATS_MEDIAN, &r)) return r; }

    // Check if it's a vector or tensor. If it's empty or not a list, return NULL.
    if (!ll_is_seq(data)) return expr_copy(res);

    LLReader e;
    ll_rd_open(data, &e);
    size_t n = e.n;
    if (n == 0) { ll_rd_close(&e); return expr_copy(res); }

    // Check if it's a matrix/tensor by checking if the first element is a List.
    if (data->type != EXPR_LINKEDLIST &&
        data->data.function.args[0]->type == EXPR_FUNCTION &&
        data->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
        data->data.function.args[0]->data.function.head->data.symbol.name == SYM_List) {
        // Apply columnwise via Transpose and Map
        ll_rd_close(&e);
        Expr* r = stats_apply_columnwise("Median", data);
        if (!r) stats_msg_rectn("Median", res);    /* ragged / mixed depth */
        return r ? r : expr_copy(res);
    }

    // Now it's a 1D vector. Check if all elements are real numbers.
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
        mth_message("Median", "rectn", "Rectangular array of real numbers is expected at position 1 in %s.", str);
        free(str);
        return expr_copy(res);
    }

    /* Numeric order (stats_sort_order: by value, not Sort's canonical order),
     * as a permutation over the elements where they lie: no Sort call and no
     * element array. */
    size_t* perm = stats_sort_order(&e);
    if (!perm) { ll_rd_close(&e); return expr_copy(res); }
    StatsOrder ord = { &e, perm };

    Expr* result = NULL;
    if (n % 2 == 1) {
        result = expr_copy(stats_order_at(&ord, n / 2));
    } else {
        Expr* m1 = expr_copy(stats_order_at(&ord, n / 2 - 1));
        Expr* m2 = expr_copy(stats_order_at(&ord, n / 2));
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){m1, m2}, 2);
        Expr* sum_eval = evaluate(sum);
        expr_free(sum);

        Expr* div = expr_new_function(expr_new_symbol(SYM_Divide), (Expr*[]){sum_eval, expr_new_integer(2)}, 2);
        result = evaluate(div);
        expr_free(div);
    }

    free(perm);
    ll_rd_close(&e);
    return result;
}
