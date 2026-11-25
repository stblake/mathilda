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
#include "numeric.h"   /* numeric_machine_subnormal */
#include "linkedlist.h" /* LLReader */
#include <stdlib.h>
#include "dist.h"    /* ml_dist_property: distribution moments */
#include "common.h"  /* head_is */

/* Mean of the n elements v[0..n): the vector core, shared by a List's
 * arguments and a LinkedList's elements (LLReader). */
Expr* stats_mean_vec(LLReader* rd) {
    size_t n = rd->n;
    if (n == 0) return NULL;

    // Check if all elements are purely real numeric
    bool all_numeric = true;
    bool has_real = false;

    for (size_t i = 0; i < n; i++) {
        Expr* elem = ll_rd_at(rd, i);
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
                stats_is_numeric(ll_rd_at(rd, i), &v, NULL);
                sum_val += v;
            }
            /* A sum that overflowed (Mean[{1.*^308, 1.*^308}] was Overflow[])
             * or a mean below DBL_MIN cannot be answered here: decline to the
             * symbolic Plus path below, which keeps the value at machine
             * precision in MPFR (numeric-tower review C3). */
            double mean_val = sum_val / (double)n;
            if (stats_machine_ok(mean_val))
                return expr_new_real(mean_val);
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
                Expr* elem = ll_rd_at(rd, i);
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
    Expr* sum = eval_and_free(ll_rd_call(expr_new_symbol(SYM_Plus), rd));

    Expr* n_inv = make_rational(1, (int64_t)n);
    Expr* mean_args[2] = { n_inv, sum };
    Expr* mean_call = expr_new_function(expr_new_symbol(SYM_Times), mean_args, 2);
    Expr* result = evaluate(mean_call);
    expr_free(mean_call);
    return result;
}

Expr* builtin_mean(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* Mean[assoc] is the mean of the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* Mean of a distribution (NormalDistribution, UniformDistribution, a 1-D
     * SmoothKernel): its closed form, see src/ml/dist.h. */
    if (data->type == EXPR_FUNCTION && !head_is(data, SYM_List)) {
        Expr* r = ml_dist_property(ML_PROP_MEAN, data, NULL);
        if (r) return r;
    }

    /* NDArray fast path (see ndreduce.c). */
    if (ndred_call_has_ndarray(res)) return ndred_mean(res);

    /* LinkedList: the packed kernel on the chunks, or columnwise (stats_ll.c). */
    { Expr* r; if (stats_ll_arm(res, STATS_MEAN, &r)) return r; }

    // Check if it's a matrix
    if (data->type != EXPR_LINKEDLIST) {
        Expr* matrixq_args[1] = { expr_copy(data) };
        Expr* matrixq_call = expr_new_function(expr_new_symbol(SYM_MatrixQ), matrixq_args, 1);
        Expr* is_matrix = evaluate(matrixq_call);
        expr_free(matrixq_call);
        bool mat = is_matrix->type == EXPR_SYMBOL && is_matrix->data.symbol.name == SYM_True;
        expr_free(is_matrix);
        if (mat) return stats_apply_columnwise("Mean", data);
    }

    /* A vector: a List's arguments or a LinkedList's elements. */
    LLReader e;
    if (!ll_rd_open(data, &e)) return NULL;
    Expr* r = stats_mean_vec(&e);
    ll_rd_close(&e);
    return r;
}
