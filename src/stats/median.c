/* median.c -- Median[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "sym_names.h"
#include "assoc.h"
#include "ndreduce.h"
#include "pack.h"      /* pack_eval_plain — the internal Sort can return a buffer */
#include "print.h"     /* expr_to_string */
#include <stdio.h>
#include <stdlib.h>

Expr* builtin_median(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* Median[assoc] is the median of the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* NDArray fast path: quickselect on the flat buffer (see ndreduce.c). */
    if (ndred_call_has_ndarray(res)) return ndred_median(res);

    // Check if it's a vector or tensor. If it's empty or not a list, return NULL.
    if (data->type != EXPR_FUNCTION || data->data.function.head->type != EXPR_SYMBOL || data->data.function.head->data.symbol.name != SYM_List) {
        return expr_copy(res);
    }

    size_t n = data->data.function.arg_count;
    if (n == 0) return expr_copy(res);

    // Check if it's a matrix/tensor by checking if the first element is a List.
    if (data->data.function.args[0]->type == EXPR_FUNCTION &&
        data->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
        data->data.function.args[0]->data.function.head->data.symbol.name == SYM_List) {
        // Apply columnwise via Transpose and Map
        return stats_apply_columnwise("Median", data);
    }

    // Now it's a 1D vector. Check if all elements are real numbers.
    bool all_real = true;
    for (size_t i = 0; i < n; i++) {
        Expr* elem = data->data.function.args[i];
        if (!stats_is_real_numeric(elem)) {
            all_real = false;
            break;
        }
    }

    if (!all_real) {
        char* str = expr_to_string(res);
        printf("Median::rectn: Rectangular array of real numbers is expected at position 1 in %s.\n", str);
        free(str);
        return expr_copy(res);
    }

    // Sort the list
    /* pack_eval_plain, not evaluate: Sort now returns a PACKED list for a large
     * machine-number input, and the args walk below reads .data.function.args.
     * The gate only covers a builtin's arguments, never what its own internal
     * evaluate() hands back. */
    Expr* sort_expr = expr_new_function(expr_new_symbol(SYM_Sort), (Expr*[]){expr_copy(data)}, 1);
    Expr* sorted = pack_eval_plain(sort_expr);
    expr_free(sort_expr);

    if (sorted->type != EXPR_FUNCTION || sorted->data.function.arg_count != n) {
        expr_free(sorted);
        return expr_copy(res);
    }

    Expr* result = NULL;
    if (n % 2 == 1) {
        result = expr_copy(sorted->data.function.args[n / 2]);
    } else {
        Expr* m1 = sorted->data.function.args[n / 2 - 1];
        Expr* m2 = sorted->data.function.args[n / 2];
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(m1), expr_copy(m2)}, 2);
        Expr* sum_eval = evaluate(sum);
        expr_free(sum);

        Expr* div = expr_new_function(expr_new_symbol(SYM_Divide), (Expr*[]){sum_eval, expr_new_integer(2)}, 2);
        result = evaluate(div);
        expr_free(div);
    }

    expr_free(sorted);
    return result;
}
