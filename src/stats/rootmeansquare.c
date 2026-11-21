/* rootmeansquare.c -- RootMeanSquare[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "arithmetic.h"
#include "sym_names.h"
#include "ndreduce.h"
#include "linkedlist.h" /* LLReader */
#include <stdlib.h>
#include <math.h>

static Expr* rms_vec(LLReader* rd) {
    size_t n = rd->n;
    if (n == 0) return NULL;

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

    if (all_numeric && has_real) {
        double sum_sq = 0;
        for (size_t i = 0; i < n; i++) {
            double v = 0.0;
            stats_is_numeric(ll_rd_at(rd, i), &v, NULL);
            sum_sq += v * v;
        }
        double r = sqrt(sum_sq / (double)n);
        /* Squares that overflowed, or underflowed to 0 from nonzero data
         * (RootMeanSquare[{1.*^-200, 1.*^-200}] was 0.), decline to the
         * symbolic path below (numeric-tower review C3). */
        bool any_nonzero = false;
        for (size_t i = 0; i < n && !any_nonzero; i++) {
            double v = 0.0;
            stats_is_numeric(ll_rd_at(rd, i), &v, NULL);
            any_nonzero = v != 0.0;
        }
        if (stats_machine_ok(r) && (r != 0.0 || !any_nonzero)) return expr_new_real(r);
    }

    // Fallback to symbolic (also used for exact rationals/integers to ensure correct square root distribution)
    Expr** sq_args = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        Expr* power_args[2] = { expr_copy(ll_rd_at(rd, i)), expr_new_integer(2) };
        Expr* power_call = expr_new_function(expr_new_symbol(SYM_Power), power_args, 2);
        sq_args[i] = evaluate(power_call);
        expr_free(power_call);
    }

    Expr* sum_call = expr_new_function(expr_new_symbol(SYM_Plus), sq_args, n);
    Expr* sum = evaluate(sum_call);
    expr_free(sum_call);
    free(sq_args);

    bool sum_is_numeric = (sum->type == EXPR_INTEGER || sum->type == EXPR_REAL || sum->type == EXPR_BIGINT || is_rational(sum, NULL, NULL));

    int64_t root_n = (int64_t)round(sqrt((double)n));
    if (!sum_is_numeric && root_n * root_n == (int64_t)n) {
        Expr* sqrt_sum_args[2] = { sum, make_rational(1, 2) };
        Expr* sqrt_sum_call = expr_new_function(expr_new_symbol(SYM_Power), sqrt_sum_args, 2);
        Expr* sqrt_sum = evaluate(sqrt_sum_call);
        expr_free(sqrt_sum_call);

        Expr* n_inv = make_rational(1, root_n);
        Expr* times_args[2] = { n_inv, sqrt_sum };
        Expr* times_call = expr_new_function(expr_new_symbol(SYM_Times), times_args, 2);
        Expr* result = evaluate(times_call);
        expr_free(times_call);
        return result;
    } else {
        Expr* n_inv = make_rational(1, (int64_t)n);
        Expr* mean_sq_args[2] = { n_inv, sum };
        Expr* mean_sq_call = expr_new_function(expr_new_symbol(SYM_Times), mean_sq_args, 2);
        Expr* mean_sq = evaluate(mean_sq_call);
        expr_free(mean_sq_call);

        int64_t num, den;
        if (is_rational(mean_sq, &num, &den)) {
            int64_t root_den = (int64_t)round(sqrt((double)den));
            if (root_den * root_den == den && root_den > 1) {
                Expr* sqrt_num_args[2] = { expr_new_integer(num), make_rational(1, 2) };
                Expr* sqrt_num_call = expr_new_function(expr_new_symbol(SYM_Power), sqrt_num_args, 2);
                Expr* sqrt_num = evaluate(sqrt_num_call);
                expr_free(sqrt_num_call);

                Expr* times_args[2] = { make_rational(1, root_den), sqrt_num };
                Expr* times_call = expr_new_function(expr_new_symbol(SYM_Times), times_args, 2);
                Expr* result = evaluate(times_call);
                expr_free(times_call);
                expr_free(mean_sq);
                return result;
            }
        }

        Expr* sqrt_args[2] = { mean_sq, make_rational(1, 2) };
        Expr* sqrt_call = expr_new_function(expr_new_symbol(SYM_Power), sqrt_args, 2);
        Expr* result = evaluate(sqrt_call);
        expr_free(sqrt_call);
        return result;
    }
}

Expr* builtin_rootmeansquare(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* NDArray fast path (see ndreduce.c). */
    if (ndred_call_has_ndarray(res)) return ndred_rms(res);

    /* LinkedList: the packed kernel on the chunks, or columnwise (stats_ll.c). */
    { Expr* r; if (stats_ll_arm(res, STATS_RMS, &r)) return r; }

    if (data->type != EXPR_LINKEDLIST) {
        Expr* matrixq_args[1] = { expr_copy(data) };
        Expr* matrixq_call = expr_new_function(expr_new_symbol(SYM_MatrixQ), matrixq_args, 1);
        Expr* is_matrix = evaluate(matrixq_call);
        expr_free(matrixq_call);
        bool mat = is_matrix->type == EXPR_SYMBOL && is_matrix->data.symbol.name == SYM_True;
        expr_free(is_matrix);
        if (mat) return stats_apply_columnwise("RootMeanSquare", data);
    }

    /* A vector: a List's arguments or a LinkedList's elements. */
    LLReader e;
    if (!ll_rd_open(data, &e)) return NULL;
    Expr* r = rms_vec(&e);
    ll_rd_close(&e);
    return r;
}
