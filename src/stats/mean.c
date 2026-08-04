/* mean.c -- Mean[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "arithmetic.h"
#include "sym_names.h"
#include "assoc.h"
#include "ndreduce.h"

Expr* builtin_mean(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* Mean[assoc] is the mean of the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* NDArray fast path (see ndreduce.c). */
    if (ndred_call_has_ndarray(res)) return ndred_mean(res);

    // Check if it's a matrix
    Expr* matrixq_args[1] = { expr_copy(data) };
    Expr* matrixq_call = expr_new_function(expr_new_symbol(SYM_MatrixQ), matrixq_args, 1);
    Expr* is_matrix = evaluate(matrixq_call);
    expr_free(matrixq_call);
    if (is_matrix->type == EXPR_SYMBOL && is_matrix->data.symbol.name == SYM_True) {
        expr_free(is_matrix);
        return stats_apply_columnwise("Mean", data);
    }
    expr_free(is_matrix);

    // Check if it's a vector (list)
    Expr* listq_args[1] = { expr_copy(data) };
    Expr* listq_call = expr_new_function(expr_new_symbol(SYM_ListQ), listq_args, 1);
    Expr* is_list = evaluate(listq_call);
    expr_free(listq_call);
    if (is_list->type == EXPR_SYMBOL && is_list->data.symbol.name == SYM_False) {
        expr_free(is_list);
        return NULL;
    }
    expr_free(is_list);

    size_t n = data->data.function.arg_count;
    if (n == 0) return NULL;

    // Check if all elements are purely real numeric
    bool all_numeric = true;
    bool has_real = false;

    for (size_t i = 0; i < n; i++) {
        Expr* elem = data->data.function.args[i];
        if (elem->type == EXPR_REAL) {
            has_real = true;
        } else if (elem->type == EXPR_INTEGER) {
            // Keep as rational n/1
        } else if (is_rational(elem, NULL, NULL)) {
            // Keep as rational
        } else {
            all_numeric = false;
            break;
        }
    }

    if (all_numeric) {
        if (has_real) {
            double sum_val = 0;
            for (size_t i = 0; i < n; i++) {
                double v = 0.0;
                stats_is_numeric(data->data.function.args[i], &v, NULL);
                sum_val += v;
            }
            return expr_new_real(sum_val / (double)n);
        } else {
            // Exact rational arithmetic
            int64_t sum_n = 0;
            int64_t sum_d = 1;
            for (size_t i = 0; i < n; i++) {
                Expr* elem = data->data.function.args[i];
                int64_t cur_n, cur_d;
                if (elem->type == EXPR_INTEGER) {
                    cur_n = elem->data.integer;
                    cur_d = 1;
                } else {
                    is_rational(elem, &cur_n, &cur_d);
                }
                // sum = sum_n/sum_d + cur_n/cur_d = (sum_n*cur_d + cur_n*sum_d) / (sum_d*cur_d)
                int64_t new_n = sum_n * cur_d + cur_n * sum_d;
                int64_t new_d = sum_d * cur_d;
                int64_t common = gcd(new_n < 0 ? -new_n : new_n, new_d);
                sum_n = new_n / common;
                sum_d = new_d / common;
            }
            // mean = sum / n = sum_n / (sum_d * n)
            return make_rational(sum_n, sum_d * (int64_t)n);
        }
    }

    // Fallback to symbolic: Mean = (1/n) * Plus @@ data
    Expr* sum_call = expr_new_function(expr_new_symbol(SYM_Apply), (Expr*[]){expr_new_symbol(SYM_Plus), expr_copy(data)}, 2);
    Expr* sum = evaluate(sum_call);
    expr_free(sum_call);

    Expr* n_inv = make_rational(1, (int64_t)n);
    Expr* mean_args[2] = { n_inv, sum };
    Expr* mean_call = expr_new_function(expr_new_symbol(SYM_Times), mean_args, 2);
    Expr* result = evaluate(mean_call);
    expr_free(mean_call);
    return result;
}
