/* standarddeviation.c -- StandardDeviation[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "arithmetic.h"
#include "sym_names.h"
#include "assoc.h"
#include "ndreduce.h"
#include <math.h>

Expr* builtin_standard_deviation(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* StandardDeviation[assoc] uses the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* NDArray fast path (see ndreduce.c). */
    if (ndred_call_has_ndarray(res)) return ndred_std(res);

    // Check if it's a matrix
    Expr* matrixq_args[1] = { expr_copy(data) };
    Expr* matrixq_call = expr_new_function(expr_new_symbol(SYM_MatrixQ), matrixq_args, 1);
    Expr* is_matrix = evaluate(matrixq_call);
    expr_free(matrixq_call);
    if (is_matrix->type == EXPR_SYMBOL && is_matrix->data.symbol.name == SYM_True) {
        expr_free(is_matrix);
        return stats_apply_columnwise("StandardDeviation", data);
    }
    expr_free(is_matrix);

    // Optimization for real numeric data
    if (data->type == EXPR_FUNCTION) {
        bool all_real = true;
        size_t n = data->data.function.arg_count;
        for (size_t i = 0; i < n; i++) {
            bool complex;
            if (!stats_is_numeric(data->data.function.args[i], NULL, &complex) || complex) {
                all_real = false;
                break;
            }
        }
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

    Expr* sqrt_call = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ var, make_rational(1, 2) }, 2);
    Expr* result = evaluate(sqrt_call);
    expr_free(sqrt_call);
    return result;
}
