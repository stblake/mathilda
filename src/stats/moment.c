/* moment.c -- Moment[] (raw / power moment).
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout.
 *
 * Moment[data, r]              — the r-th raw (power) moment,
 *                                mu_r = (1/n) Sum[x_i^r]. For a matrix / array the
 *                                reduction is columnwise over the first axis
 *                                (equivalently ArrayReduce[Moment[#,r]&, x, 1]).
 * Moment[data, {r1, ..., rm}]  — the multivariate mixed raw moment,
 *                                (1/n) Sum_i Product_j x[[i,j]]^r_j,
 *                                summing the first axis and taking a product over
 *                                the second (its length must equal Length[{r1,...}]).
 *
 * The raw moment is CentralMoment without the mean subtraction. Because there is
 * no mean to subtract, Mean[data^r] threads correctly for a vector, a matrix
 * (columnwise), AND a higher-rank array in a single expression — Power is
 * Listable so data^r threads elementwise at every rank, and the outer Mean
 * collapses the first axis by n. So (unlike CentralMoment, whose data - Mean[data]
 * would thread row-wise) the scalar-order case needs no separate columnwise
 * routine. Numeric real vectors take a tight C loop (and packed / NDArray inputs a
 * machine-buffer fast path via ndred_moment); every other case — exact, symbolic,
 * matrix/array, multivariate — is built as an expression and handed to the
 * evaluator, which already knows how to be exact or symbolic. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "arithmetic.h"   /* is_rational */
#include "sym_names.h"
#include "assoc.h"
#include "ndreduce.h"
#include <stdlib.h>

/* Integer power by squaring on a double (exponent >= 0). */
static double mom_ipow(double b, int64_t e) {
    double r = 1.0;
    while (e > 0) {
        if (e & 1) r *= b;
        b *= b;
        e >>= 1;
    }
    return r;
}

/* True when e is a List[...] expression (a vector, matrix row, or subarray). */
static bool mom_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* Exact / symbolic base case: Mean[data^r]. Power is Listable, so data^r threads
 * over the elements at every rank (a vector's scalars, a matrix's / array's
 * entries), and the outer Mean divides the first axis by n — giving the vector's
 * scalar moment, the matrix's columnwise vector of moments, or the array's
 * columnwise array. Reuses Mean's exact-rational and symbolic paths (and the
 * evaluator's MPFR path for arbitrary-precision reals). */
static Expr* mom_symbolic(Expr* data, Expr* order) {
    Expr* powered = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_copy(data), expr_copy(order) }, 2));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Mean),
               (Expr*[]){ powered }, 1));
}

/* Multivariate order {r1, ..., rm} over an array of depth >= 2. For each
 * first-axis slice sub, term = Times @@ sub^rvec: the elementwise Power threads
 * rvec over the slice's second axis (raising block j to r_j) and Times @@ takes
 * the product over that axis, leaving a rank k-2 block. The result is Mean over
 * the slices (dividing the first axis by n). Requires the slice's second-axis
 * length to equal Length[rvec]. Mirrors CentralMoment's multivariate case with
 * the centering removed. */
static Expr* mom_multivariate(Expr* data, Expr* rvec) {
    Expr* first = data->data.function.args[0];
    if (!mom_is_list(first)) return NULL;                 /* need depth >= 2 */
    size_t m = rvec->data.function.arg_count;
    if (first->data.function.arg_count != m) return NULL; /* second axis must match order length */

    size_t n = data->data.function.arg_count;
    Expr** terms = malloc(sizeof(Expr*) * n);
    if (!terms) return NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* sub = data->data.function.args[i];             /* borrowed, rank k-1 */
        Expr* powered = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ expr_copy(sub), expr_copy(rvec) }, 2));  /* threads r over axis 2 */
        terms[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Apply),
                       (Expr*[]){ expr_new_symbol(SYM_Times), powered }, 2)); /* Times @@ powered */
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), terms, n);
    free(terms);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Mean),
               (Expr*[]){ list }, 1));
}

Expr* builtin_moment(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* data = res->data.function.args[0];
    Expr* order = res->data.function.args[1];

    /* Moment[assoc, r] works on the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* NDArray / packed fast path (scalar integer order); see ndreduce.c. A
     * list-valued or non-integer order, an int64 buffer, or a complex buffer
     * degrades back to a nested List and re-enters here. */
    if (ndred_call_has_ndarray(res)) return ndred_moment(res);

    if (!mom_is_list(data)) return NULL;
    size_t n = data->data.function.arg_count;
    if (n == 0) return NULL;

    /* Multivariate order {r1, ..., rm}. */
    if (mom_is_list(order)) return mom_multivariate(data, order);

    /* Scalar order over a flat vector (depth 1). Fast real path: every element is
     * a real number, at least one is inexact (approximate input -> approximate
     * output), and the order is a non-negative machine integer. A matrix / array
     * (data[[1]] is a List) skips this and takes the symbolic columnwise path. */
    if (!mom_is_list(data->data.function.args[0]) &&
        order->type == EXPR_INTEGER && order->data.integer >= 0) {
        bool all_numeric = true, has_real = false;
        for (size_t i = 0; i < n; i++) {
            Expr* e = data->data.function.args[i];
            if (e->type == EXPR_REAL) has_real = true;
            else if (e->type == EXPR_INTEGER || is_rational(e, NULL, NULL)) { /* numeric */ }
            else { all_numeric = false; break; }
        }
        if (all_numeric && has_real) {
            int64_t r = order->data.integer;
            double acc = 0.0;
            for (size_t i = 0; i < n; i++) {
                double v = 0.0;
                stats_is_numeric(data->data.function.args[i], &v, NULL);
                acc += mom_ipow(v, r);
            }
            return expr_new_real(acc / (double)n);
        }
    }

    /* Exact / symbolic vector, matrix, or array; or symbolic / non-integer order. */
    return mom_symbolic(data, order);
}
