/* movingaverage.c -- MovingAverage[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "eval.h"
#include "sym_names.h"
#include "ndarray.h"
#include "ndreduce.h"
#include <stdlib.h>
#include <gmp.h>

/*
 * MovingAverage[list, r]                — simple moving average over runs of r elements.
 * MovingAverage[list, {w1, ..., wr}]    — weighted moving average with weights w_k
 *                                         (output uses w_k / Sum[w_k] as effective weights).
 *
 * Output length is Length[list] - r + 1 when 1 <= r <= Length[list]; the function
 * stays unevaluated otherwise.  The unweighted form delegates to Mean so it inherits
 * Mean's exact rational / bigint / real / symbolic handling.  The weighted form builds
 * Plus[Times[w_k/wsum, x_{i+k}], ...] and lets the evaluator simplify.
 */
Expr* builtin_moving_average(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* data = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* NDArray fast path: sliding-window mean on the buffer (see ndreduce.c). */
    if (is_ndarray(data)) return ndred_moving_average(res);

    if (data->type != EXPR_FUNCTION ||
        data->data.function.head->type != EXPR_SYMBOL ||
        data->data.function.head->data.symbol.name != SYM_List) {
        return NULL;
    }

    size_t n = data->data.function.arg_count;
    size_t r = 0;
    Expr** weights = NULL; /* borrowed; non-NULL means weighted form */

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
    } else if (spec->type == EXPR_FUNCTION &&
               spec->data.function.head->type == EXPR_SYMBOL &&
               spec->data.function.head->data.symbol.name == SYM_List) {
        r = spec->data.function.arg_count;
        if (r == 0 || r > n) return NULL;
        weights = spec->data.function.args;
    } else {
        return NULL;
    }

    size_t out_n = n - r + 1;
    Expr** out = malloc(sizeof(Expr*) * out_n);
    if (!out) return NULL;

    if (weights == NULL) {
        for (size_t i = 0; i < out_n; i++) {
            Expr** sub = malloc(sizeof(Expr*) * r);
            for (size_t k = 0; k < r; k++) {
                sub[k] = expr_copy(data->data.function.args[i + k]);
            }
            Expr* sublist = expr_new_function(expr_new_symbol(SYM_List), sub, r);
            free(sub);
            out[i] = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Mean), (Expr*[]){ sublist }, 1));
        }
    } else {
        Expr** w_copies = malloc(sizeof(Expr*) * r);
        for (size_t k = 0; k < r; k++) w_copies[k] = expr_copy(weights[k]);
        Expr* wsum = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Plus), w_copies, r));
        free(w_copies);

        Expr* wsum_inv = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Power),
            (Expr*[]){ wsum, expr_new_integer(-1) }, 2));

        Expr** coefs = malloc(sizeof(Expr*) * r);
        for (size_t k = 0; k < r; k++) {
            coefs[k] = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_copy(weights[k]), expr_copy(wsum_inv) }, 2));
        }
        expr_free(wsum_inv);

        for (size_t i = 0; i < out_n; i++) {
            Expr** terms = malloc(sizeof(Expr*) * r);
            for (size_t k = 0; k < r; k++) {
                terms[k] = eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_copy(coefs[k]),
                               expr_copy(data->data.function.args[i + k]) }, 2));
            }
            out[i] = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Plus), terms, r));
            free(terms);
        }

        for (size_t k = 0; k < r; k++) expr_free(coefs[k]);
        free(coefs);
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, out_n);
    free(out);
    return result;
}
