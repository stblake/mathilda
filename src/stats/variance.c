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
#include "linkedlist.h" /* LLReader */
#include <stdlib.h>
#include "dist.h"    /* ml_dist_property: distribution moments */
#include "common.h"  /* head_is */

Expr* stats_variance_vec(LLReader* rd) {
    size_t n = rd->n;
    if (n <= 1) return NULL;

    // Check if all elements are purely real numeric
    bool all_numeric = true;
    bool has_real = false;
    for (size_t i = 0; i < n; i++) {
        Expr* elem = ll_rd_at(rd, i);
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
            /* The machine kernel's own two-pass reduction, not a second
             * algorithm: a boxed list too short to pack must give the packed,
             * NDArray and Compile[] answer bit for bit (Welford's one-pass
             * update, used here before, differed in the last ulp on glibc). */
            double* xs = malloc(sizeof(double) * n);
            if (!xs) return NULL;
            for (size_t i = 0; i < n; i++) {
                xs[i] = 0.0;
                stats_is_numeric(ll_rd_at(rd, i), &xs[i], NULL);
            }
            double v = ndred_variance_f64(xs, n);
            /* Squares that underflowed to 0. (distinct data of size
             * 1.*^-200) are no zero variance: decline, as for an overflow. */
            bool underflow = false;
            if (v == 0.0)
                for (size_t i = 1; i < n && !underflow; i++) underflow = xs[i] != xs[0];
            free(xs);
            /* Overflowed squares decline to the symbolic path (C3). */
            if (stats_machine_ok(v) && !underflow) return expr_new_real(v);
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
                Expr* elem = ll_rd_at(rd, i);
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
                Expr* elem = ll_rd_at(rd, i);
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
    Expr* mu = stats_mean_vec(rd);
    if (!mu) return NULL;

    // Sum[(x - mu) * Conjugate[x - mu]]
    Expr** diffs = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* x = ll_rd_at(rd, i);
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

Expr* builtin_variance(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* Variance[assoc] is the variance of the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* Variance of a distribution (NormalDistribution, UniformDistribution, a 1-D
     * SmoothKernel): its closed form, see src/ml/dist.h. */
    if (data->type == EXPR_FUNCTION && !head_is(data, SYM_List)) {
        Expr* r = ml_dist_property(ML_PROP_VARIANCE, data, NULL);
        if (r) return r;
    }

    /* NDArray fast path (see ndreduce.c). */
    if (ndred_call_has_ndarray(res)) return ndred_variance(res);

    /* LinkedList: the packed kernel on the chunks, or columnwise (stats_ll.c). */
    { Expr* r; if (stats_ll_arm(res, STATS_VARIANCE, &r)) return r; }

    if (stats_reject_short("Variance", data)) return NULL;

    // Check if it's a matrix
    if (data->type != EXPR_LINKEDLIST) {
        Expr* matrixq_args[1] = { expr_copy(data) };
        Expr* matrixq_call = expr_new_function(expr_new_symbol(SYM_MatrixQ), matrixq_args, 1);
        Expr* is_matrix = evaluate(matrixq_call);
        expr_free(matrixq_call);
        bool mat = is_matrix->type == EXPR_SYMBOL && is_matrix->data.symbol.name == SYM_True;
        expr_free(is_matrix);
        if (mat) return stats_apply_columnwise("Variance", data);
        if (stats_reject_ragged("Variance", res, data)) return NULL;
    }

    /* A vector: a List's arguments or a LinkedList's elements. */
    LLReader e;
    if (!ll_rd_open(data, &e)) return NULL;
    Expr* r = stats_variance_vec(&e);
    ll_rd_close(&e);
    return r;
}
