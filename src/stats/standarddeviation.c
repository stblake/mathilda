/* standarddeviation.c -- StandardDeviation[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "arithmetic.h"
#include "sym_names.h"
#include "assoc.h"
#include "ndreduce.h"
#include "linkedlist.h" /* LLReader */
#include <math.h>
#include <stdlib.h>

Expr* builtin_standard_deviation(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* StandardDeviation[assoc] uses the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* NDArray fast path (see ndreduce.c). */
    if (ndred_call_has_ndarray(res)) return ndred_std(res);

    /* LinkedList: the packed kernel on the chunks, or columnwise (stats_ll.c). */
    { Expr* r; if (stats_ll_arm(res, STATS_STD, &r)) return r; }

    if (stats_reject_short("StandardDeviation", data)) return NULL;

    // Check if it's a matrix
    if (data->type != EXPR_LINKEDLIST) {
        Expr* matrixq_args[1] = { expr_copy(data) };
        Expr* matrixq_call = expr_new_function(expr_new_symbol(SYM_MatrixQ), matrixq_args, 1);
        Expr* is_matrix = evaluate(matrixq_call);
        expr_free(matrixq_call);
        bool mat = is_matrix->type == EXPR_SYMBOL && is_matrix->data.symbol.name == SYM_True;
        expr_free(is_matrix);
        if (mat) return stats_apply_columnwise("StandardDeviation", data);
        if (stats_reject_ragged("StandardDeviation", res, data)) return NULL;
    }

    // Optimization for real numeric data: a List's arguments or a LinkedList's elements
    LLReader e;
    bool seq = false;
    if (data->type == EXPR_LINKEDLIST) {
        seq = ll_rd_open(data, &e);
    } else if (data->type == EXPR_FUNCTION) {
        ll_rd_wrap(&e, data->data.function.args, data->data.function.arg_count);
        seq = true;
    }
    if (seq) {
        bool all_real = true;
        size_t n = e.n;
        for (size_t i = 0; i < n; i++) {
            bool complex;
            if (!stats_is_numeric(ll_rd_at(&e, i), NULL, &complex) || complex) {
                all_real = false;
                break;
            }
        }
        /* Machine data (a Real among them): the buffer kernel's own
         * reduction, so the answer is the packed / NDArray one bit for bit,
         * and extreme magnitudes are rescaled rather than squared to 0. or
         * Infinity (StandardDeviation[{1.*^-200, 2.*^-200, 3.*^-200}]). */
        bool any_real = false;
        for (size_t i = 0; all_real && i < n; i++)
            if (ll_rd_at(&e, i)->type == EXPR_REAL) any_real = true;
        if (all_real && any_real && n > 1) {
            double* xs = malloc(sizeof(double) * n);
            if (xs) {
                for (size_t i = 0; i < n; i++) {
                    xs[i] = 0.0;
                    stats_is_numeric(ll_rd_at(&e, i), &xs[i], NULL);
                }
                double v = ndred_std_f64(xs, n);
                free(xs);
                if (stats_machine_ok(v)) { ll_rd_close(&e); return expr_new_real(v); }
            }
        }
        ll_rd_close(&e);
        if (all_real && n > 1) {
            Expr* var_call = expr_new_function(expr_new_symbol(SYM_Variance), (Expr*[]){ expr_copy(data) }, 1);
            Expr* var_e = evaluate(var_call);
            expr_free(var_call);
            if (var_e->type == EXPR_REAL) {
                double val = sqrt(var_e->data.real);
                expr_free(var_e);
                return expr_new_real(val);
            }
            expr_free(var_e);
        }
    }

    // StandardDeviation[data] = Variance[data]^(1/2)
    Expr* var_call = expr_new_function(expr_new_symbol(SYM_Variance), (Expr*[]){ expr_copy(data) }, 1);
    Expr* var = evaluate(var_call);
    expr_free(var_call);

    if (!var) return NULL;
    /* Variance declined (unevaluated): so does StandardDeviation. */
    if (var->type == EXPR_FUNCTION && var->data.function.head->type == EXPR_SYMBOL &&
        var->data.function.head->data.symbol.name == SYM_Variance) {
        expr_free(var);
        return NULL;
    }

    Expr* sqrt_call = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ var, make_rational(1, 2) }, 2);
    Expr* result = evaluate(sqrt_call);
    expr_free(sqrt_call);
    return result;
}
