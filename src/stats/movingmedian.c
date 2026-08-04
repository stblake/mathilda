/* movingmedian.c -- MovingMedian[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "sym_names.h"
#include "ndarray.h"
#include "ndreduce.h"
#include "print.h"     /* expr_to_string */
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

/*
 * MovingMedian[list, r]
 *   gives the moving median of list, computed using spans of r elements.
 *
 * Output length is Length[list] - r + 1 when 1 <= r <= Length[list]; the
 * function stays unevaluated otherwise. Operates on real-valued vectors and
 * matrices; for matrices each row-window is reduced via Median, which yields
 * a column-wise median vector. Non-numeric data triggers the MovingMedian::arg1
 * message and the expression remains unevaluated.
 */
Expr* builtin_moving_median(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* data = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* NDArray fast path: sliding-window median on the buffer (see ndreduce.c). */
    if (is_ndarray(data)) return ndred_moving_median(res);

    if (data->type != EXPR_FUNCTION ||
        data->data.function.head->type != EXPR_SYMBOL ||
        data->data.function.head->data.symbol.name != SYM_List) {
        return NULL;
    }

    size_t n = data->data.function.arg_count;
    size_t r = 0;
    if (spec->type == EXPR_INTEGER) {
        if (spec->data.integer < 1) return NULL;
        if ((uint64_t)spec->data.integer > (uint64_t)n) return NULL;
        r = (size_t)spec->data.integer;
    } else if (spec->type == EXPR_BIGINT) {
        if (mpz_sgn(spec->data.bigint) <= 0) return NULL;
        if (!mpz_fits_ulong_p(spec->data.bigint)) return NULL;
        unsigned long rr = mpz_get_ui(spec->data.bigint);
        if ((size_t)rr > n) return NULL;
        r = (size_t)rr;
    } else {
        return NULL;
    }

    if (n == 0) return NULL;

    /* Decide vector vs matrix based on whether the first element is a List. */
    bool matrix_mode = (data->data.function.args[0]->type == EXPR_FUNCTION &&
                        data->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
                        data->data.function.args[0]->data.function.head->data.symbol.name == SYM_List);

    /* Validate that every leaf is a real-valued numeric. Matrices must be rectangular. */
    bool ok = true;
    if (matrix_mode) {
        size_t cols = data->data.function.args[0]->data.function.arg_count;
        for (size_t i = 0; i < n && ok; i++) {
            Expr* row = data->data.function.args[i];
            if (row->type != EXPR_FUNCTION ||
                row->data.function.head->type != EXPR_SYMBOL ||
                row->data.function.head->data.symbol.name != SYM_List ||
                row->data.function.arg_count != cols) {
                ok = false;
                break;
            }
            for (size_t j = 0; j < cols && ok; j++) {
                if (!stats_is_real_numeric(row->data.function.args[j])) ok = false;
            }
        }
    } else {
        for (size_t i = 0; i < n && ok; i++) {
            if (!stats_is_real_numeric(data->data.function.args[i])) ok = false;
        }
    }

    if (!ok) {
        char* str = expr_to_string(data);
        printf("MovingMedian::arg1: The first argument %s must be a vector or matrix of real values.\n", str);
        free(str);
        return expr_copy(res);
    }

    size_t out_n = n - r + 1;
    Expr** out = malloc(sizeof(Expr*) * out_n);
    if (!out) return NULL;

    for (size_t i = 0; i < out_n; i++) {
        Expr** sub = malloc(sizeof(Expr*) * r);
        for (size_t k = 0; k < r; k++) {
            sub[k] = expr_copy(data->data.function.args[i + k]);
        }
        Expr* sublist = expr_new_function(expr_new_symbol(SYM_List), sub, r);
        free(sub);
        out[i] = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Median), (Expr*[]){ sublist }, 1));
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, out_n);
    free(out);
    return result;
}
