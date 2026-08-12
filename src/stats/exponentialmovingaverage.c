/* exponentialmovingaverage.c -- ExponentialMovingAverage[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "sym_names.h"
#include "ndarray.h"
#include "ndreduce.h"
#include <stdlib.h>

/*
 * ExponentialMovingAverage[list, alpha]
 *   gives the exponential moving average of list with smoothing constant alpha.
 *
 * The recurrence is y[1] = x[1], y[i+1] = y[i] + alpha * (x[i+1] - y[i]).
 * Output has the same length as the input list. Stays unevaluated when the
 * first argument is not a List, or when the list is empty.
 *
 * Two evaluation strategies:
 *   1. Fast path (machine precision): if at least one element of the list or
 *      alpha itself is a real number, we evaluate the recurrence in C using
 *      doubles. This is O(n) with no Expr allocations beyond the output and
 *      keeps numeric stability comparable to Mathematica's N[...] form.
 *   2. Symbolic path: build the recurrence as Plus / Times nodes and let the
 *      evaluator handle exact rational arithmetic, bignum promotion, and
 *      symbolic simplification. Works for arbitrary alpha (including symbolic).
 */
Expr* builtin_exponential_moving_average(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* data = res->data.function.args[0];
    Expr* alpha = res->data.function.args[1];

    /* NDArray fast path: EMA recurrence on the buffer (see ndreduce.c). */
    if (is_ndarray(data)) return ndred_ema(res);

    if (data->type != EXPR_FUNCTION ||
        data->data.function.head->type != EXPR_SYMBOL ||
        data->data.function.head->data.symbol.name != SYM_List) {
        return NULL;
    }

    size_t n = data->data.function.arg_count;
    if (n == 0) return NULL;

    /* Decide whether the fast double-precision path applies: at least one
     * element of list or alpha must be EXPR_REAL, and every list element plus
     * alpha must be a real-valued numeric (Integer / Real / Rational, no
     * Complex, no symbolic). Bignums fall through to the symbolic path so we
     * don't lose precision in cases where the exact value is wanted. */
    bool any_real = (alpha->type == EXPR_REAL);
    if (!any_real) {
        for (size_t i = 0; i < n; i++) {
            if (data->data.function.args[i]->type == EXPR_REAL) {
                any_real = true;
                break;
            }
        }
    }
    bool fast_ok = any_real;
    if (fast_ok) {
        bool cplx = false;
        if (!stats_is_numeric(alpha, NULL, &cplx) || cplx) fast_ok = false;
    }
    if (fast_ok) {
        for (size_t i = 0; i < n && fast_ok; i++) {
            bool cplx = false;
            if (!stats_is_numeric(data->data.function.args[i], NULL, &cplx) || cplx) {
                fast_ok = false;
            }
        }
    }

    if (fast_ok) {
        double a = 0.0;
        stats_is_numeric(alpha, &a, NULL);
        Expr** out = malloc(sizeof(Expr*) * n);
        if (!out) return NULL;
        double y = 0.0;
        stats_is_numeric(data->data.function.args[0], &y, NULL);
        out[0] = expr_new_real(y);
        for (size_t i = 1; i < n; i++) {
            double x = 0.0;
            stats_is_numeric(data->data.function.args[i], &x, NULL);
            y = y + a * (x - y);
            out[i] = expr_new_real(y);
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, n);
        free(out);
        return result;
    }

    /* Symbolic / exact-rational / bignum path. */
    Expr** out = malloc(sizeof(Expr*) * n);
    if (!out) return NULL;
    out[0] = expr_copy(data->data.function.args[0]);
    for (size_t i = 1; i < n; i++) {
        Expr* x_i = expr_copy(data->data.function.args[i]);
        Expr* neg_y_prev = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_new_integer(-1), expr_copy(out[i-1]) }, 2));
        Expr* diff = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Plus),
            (Expr*[]){ x_i, neg_y_prev }, 2));
        Expr* alpha_diff = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_copy(alpha), diff }, 2));
        out[i] = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy(out[i-1]), alpha_diff }, 2));
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, n);
    free(out);
    return result;
}
