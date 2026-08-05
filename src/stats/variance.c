/* variance.c -- Variance[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "arithmetic.h"
#include "checked_int.h"  /* ci_mul_i64 etc: the int64 fast path must decline, not wrap */
#include "sym_names.h"
#include "assoc.h"
#include "ndreduce.h"
#include <stdlib.h>

Expr* builtin_variance(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* Variance[assoc] is the variance of the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* NDArray fast path (see ndreduce.c). */
    if (ndred_call_has_ndarray(res)) return ndred_variance(res);

    // Check if it's a matrix
    Expr* matrixq_args[1] = { expr_copy(data) };
    Expr* matrixq_call = expr_new_function(expr_new_symbol(SYM_MatrixQ), matrixq_args, 1);
    Expr* is_matrix = evaluate(matrixq_call);
    expr_free(matrixq_call);
    if (is_matrix->type == EXPR_SYMBOL && is_matrix->data.symbol.name == SYM_True) {
        expr_free(is_matrix);
        return stats_apply_columnwise("Variance", data);
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
    if (n <= 1) return NULL;

    // Check if all elements are purely real numeric
    bool all_numeric = true;
    bool has_real = false;
    for (size_t i = 0; i < n; i++) {
        Expr* elem = data->data.function.args[i];
        if (elem->type == EXPR_REAL) {
            has_real = true;
        } else if (elem->type == EXPR_INTEGER || is_rational(elem, NULL, NULL)) {
            // Numeric
        } else {
            all_numeric = false;
            break;
        }
    }

    if (all_numeric) {
        if (has_real) {
            // Welford's algorithm
            double m = 0;
            double s = 0;
            for (size_t i = 0; i < n; i++) {
                double x = 0.0;
                stats_is_numeric(data->data.function.args[i], &x, NULL);
                double old_m = m;
                m += (x - m) / (double)(i + 1);
                s += (x - m) * (x - old_m);
            }
            return expr_new_real(s / (double)(n - 1));
        } else {
            /* Exact calculation for Variance, WITH OVERFLOW DETECTION.
             * Var = (Sum[x^2] - n*Mean[x]^2) / (n-1)
             * Using Sum[(n*x_i - Sum[x_j])^2] / (n^2 * (n-1))
             *
             * This squares an already-scaled numerator: term_n is ~n*value, so
             * term_n^2 is ~(n*value)^2, and the running accumulator is then
             * multiplied by that square. Unchecked, a 10^4-element integer list
             * overflowed and returned 0.0176 where the answer was 84815.8 —
             * and, being an overflowed sum of squares, sometimes returned a
             * NEGATIVE variance, which is impossible by construction. The
             * erratic-in-n behaviour (n=4096 wrong, n=4097 right) is the
             * signature: whether it overflows depends on how far each gcd
             * reduces, not on n alone.
             *
             * On overflow we fall through to the symbolic path below, which is
             * exact for any magnitude via GMP. A fast path that cannot
             * represent the answer must decline, not approximate. */
            int64_t sum_n = 0;
            int64_t sum_d = 1;
            bool overflow = false;
            for (size_t i = 0; i < n && !overflow; i++) {
                Expr* elem = data->data.function.args[i];
                int64_t cur_n, cur_d;
                if (elem->type == EXPR_INTEGER) { cur_n = elem->data.integer; cur_d = 1; }
                else is_rational(elem, &cur_n, &cur_d);
                int64_t a, b, new_n, new_d;
                if (ci_mul_i64(sum_n, cur_d, &a) ||
                    ci_mul_i64(cur_n, sum_d, &b) ||
                    ci_add_i64(a, b, &new_n) ||
                    ci_mul_i64(sum_d, cur_d, &new_d)) { overflow = true; break; }
                int64_t common = gcd(new_n < 0 ? -new_n : new_n, new_d);
                sum_n = new_n / common;
                sum_d = new_d / common;
            }
            // Now sum_n/sum_d is the sum of elements
            int64_t sq_sum_n = 0;
            int64_t sq_sum_d = 1;
            for (size_t i = 0; i < n && !overflow; i++) {
                Expr* elem = data->data.function.args[i];
                int64_t cur_n, cur_d;
                if (elem->type == EXPR_INTEGER) { cur_n = elem->data.integer; cur_d = 1; }
                else is_rational(elem, &cur_n, &cur_d);

                // (x - mean)^2 = (cur_n/cur_d - sum_n/(n*sum_d))^2
                // = ( (cur_n * n * sum_d - sum_n * cur_d) / (n * sum_d * cur_d) )^2
                int64_t t1, t2, t3, term_n, term_d;
                if (ci_mul_i64(cur_n, (int64_t)n, &t1) ||
                    ci_mul_i64(t1, sum_d, &t2) ||
                    ci_mul_i64(sum_n, cur_d, &t3) ||
                    ci_sub_i64(t2, t3, &term_n) ||
                    ci_mul_i64((int64_t)n, sum_d, &t1) ||
                    ci_mul_i64(t1, cur_d, &term_d)) { overflow = true; break; }
                int64_t common = gcd(term_n < 0 ? -term_n : term_n, term_d);
                term_n /= common; term_d /= common;

                int64_t term_sq_n, term_sq_d;
                if (ci_mul_i64(term_n, term_n, &term_sq_n) ||
                    ci_mul_i64(term_d, term_d, &term_sq_d)) { overflow = true; break; }

                int64_t a, b, new_sq_sum_n, new_sq_sum_d;
                if (ci_mul_i64(sq_sum_n, term_sq_d, &a) ||
                    ci_mul_i64(term_sq_n, sq_sum_d, &b) ||
                    ci_add_i64(a, b, &new_sq_sum_n) ||
                    ci_mul_i64(sq_sum_d, term_sq_d, &new_sq_sum_d)) {
                    overflow = true; break;
                }
                common = gcd(new_sq_sum_n < 0 ? -new_sq_sum_n : new_sq_sum_n, new_sq_sum_d);
                sq_sum_n = new_sq_sum_n / common;
                sq_sum_d = new_sq_sum_d / common;
            }
            // Variance = sq_sum / (n-1)
            int64_t den;
            if (!overflow && !ci_mul_i64(sq_sum_d, (int64_t)n - 1, &den)) {
                return make_rational(sq_sum_n, den);
            }
        }
    }

    // Fallback to symbolic: Compute Mean
    Expr* mean_args[1] = { expr_copy(data) };
    Expr* mean_call = expr_new_function(expr_new_symbol(SYM_Mean), mean_args, 1);
    Expr* mu = evaluate(mean_call);
    expr_free(mean_call);

    // Sum[(x - mu) * Conjugate[x - mu]]
    Expr** diffs = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* x = data->data.function.args[i];
        Expr* sub_args[2] = { expr_copy(x), expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), expr_copy(mu)}, 2) };
        Expr* diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), sub_args, 2));

        Expr* conj_args[1] = { expr_copy(diff) };
        Expr* conj_call = expr_new_function(expr_new_symbol(SYM_Conjugate), conj_args, 1);
        Expr* conj_diff = evaluate(conj_call);
        expr_free(conj_call);

        Expr* prod_args[2] = { diff, conj_diff };
        diffs[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), prod_args, 2));
    }
    expr_free(mu);

    Expr* sum_diffs = expr_new_function(expr_new_symbol(SYM_Plus), diffs, n);
    Expr* sum = evaluate(sum_diffs);
    expr_free(sum_diffs);
    free(diffs);

    Expr* n_minus_1_inv = make_rational(1, (int64_t)n - 1);
    Expr* var_call = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ n_minus_1_inv, sum }, 2);
    Expr* result = evaluate(var_call);
    expr_free(var_call);
    return result;
}
