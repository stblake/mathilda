/* mean.c -- Mean[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "arithmetic.h"
#include "checked_int.h"  /* ci_mul_i64 etc: the int64 fast path must decline, not wrap */
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
            /* Exact rational arithmetic in int64, WITH OVERFLOW DETECTION.
             *
             * Every product below can overflow, and unchecked it did: a list of
             * 10^4 rationals sharing a denominator of 10^4 (exactly what
             * `data - Mean[data]` produces, so exactly what CentralMoment and
             * therefore Skewness/Kurtosis feed back in) drove `sum_d * cur_d`
             * past INT64_MAX and returned a silently wrong mean — 0.0177 where
             * the answer was 84807.4. The failure was erratic rather than
             * monotone in n, because whether it overflows depends on how far the
             * running gcd reduces, which is why a small-denominator list of the
             * same length was unaffected and hid the bug.
             *
             * On overflow we do NOT return: control falls through to the
             * symbolic path below, which computes `Plus @@ data / n` through the
             * evaluator and is exact for any magnitude via GMP. The int64 block
             * is a fast path, and a fast path that cannot represent the answer
             * must decline, not approximate. */
            int64_t sum_n = 0;
            int64_t sum_d = 1;
            bool overflow = false;
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
                int64_t a, b, new_n, new_d;
                if (ci_mul_i64(sum_n, cur_d, &a) ||
                    ci_mul_i64(cur_n, sum_d, &b) ||
                    ci_add_i64(a, b, &new_n) ||
                    ci_mul_i64(sum_d, cur_d, &new_d)) { overflow = true; break; }
                int64_t common = gcd(new_n < 0 ? -new_n : new_n, new_d);
                sum_n = new_n / common;
                sum_d = new_d / common;
            }
            // mean = sum / n = sum_n / (sum_d * n)
            int64_t den;
            if (!overflow && !ci_mul_i64(sum_d, (int64_t)n, &den)) {
                return make_rational(sum_n, den);
            }
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
