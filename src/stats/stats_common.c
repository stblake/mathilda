/* stats_common.c -- helpers shared across the src/stats/ per-builtin files.
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout. */

#include "stats_common.h"
#include "arithmetic.h"   /* is_rational, make_rational */
#include "complex.h"      /* is_complex */
#include "eval.h"         /* evaluate, eval_and_free */
#include "sym_names.h"    /* SYM_* */
#include "assoc.h"        /* is_association, assoc_apply_over_values */
#include "ndreduce.h"     /* ndred_call_has_ndarray, ndred_skewness, ndred_kurtosis */
#include <math.h>         /* floor -- stats_quantile_point index fallback */

bool stats_is_numeric(Expr* e, double* val, bool* out_complex) {
    if (e->type == EXPR_INTEGER) {
        if (val) *val = (double)e->data.integer;
        if (out_complex) *out_complex = false;
        return true;
    }
    if (e->type == EXPR_REAL) {
        if (val) *val = e->data.real;
        if (out_complex) *out_complex = false;
        return true;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        if (val) *val = (double)n / (double)d;
        if (out_complex) *out_complex = false;
        return true;
    }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        if (out_complex) *out_complex = true;
        return true;
    }
    return false;
}

Expr* stats_apply_columnwise(const char* func_name, Expr* matrix) {
    // Result is Map[func_name, Transpose[matrix]]
    Expr* transpose_args[1] = { expr_copy(matrix) };
    Expr* transpose_call = expr_new_function(expr_new_symbol(SYM_Transpose), transpose_args, 1);
    Expr* transposed = evaluate(transpose_call);
    expr_free(transpose_call);

    if (transposed->type != EXPR_FUNCTION) {
        expr_free(transposed);
        return NULL;
    }

    Expr* map_args[2] = { expr_new_symbol(func_name), transposed };
    Expr* map_call = expr_new_function(expr_new_symbol(SYM_Map), map_args, 2);
    Expr* result = evaluate(map_call);
    expr_free(map_call);
    return result;
}

/* True when e is a function call with the given symbol head. */
static bool stats_head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

Expr* stats_standardized_moment(Expr* res, int p) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* data = res->data.function.args[0];

    /* Skewness[assoc] / Kurtosis[assoc] work on the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* NDArray / packed fast path (see ndreduce.c). */
    if (ndred_call_has_ndarray(res))
        return (p == 3) ? ndred_skewness(res) : ndred_kurtosis(res);

    /* m_p = CentralMoment[data, p] and m_2 = CentralMoment[data, 2]. Both reduce
     * for any list/vector/matrix/array (numeric, exact, or symbolic); if either
     * stays a CentralMoment[...] the data was not reducible, so leave the caller's
     * head unevaluated too. */
    Expr* mp = eval_and_free(expr_new_function(expr_new_symbol(SYM_CentralMoment),
                   (Expr*[]){ expr_copy(data), expr_new_integer(p) }, 2));
    Expr* m2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_CentralMoment),
                   (Expr*[]){ expr_copy(data), expr_new_integer(2) }, 2));
    if (stats_head_is(mp, SYM_CentralMoment) || stats_head_is(m2, SYM_CentralMoment)) {
        expr_free(mp); expr_free(m2);
        return NULL;
    }

    /* result = m_p / m_2^(p/2). Power and Divide thread, so a columnwise vector
     * of moments yields a columnwise vector of standardized moments. */
    Expr* denom = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ m2, make_rational(p, 2) }, 2));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Divide),
               (Expr*[]){ mp, denom }, 2));
}

bool stats_is_real_numeric(Expr* e) {
    Expr* numq = expr_new_function(expr_new_symbol(SYM_NumericQ), (Expr*[]){expr_copy(e)}, 1);
    Expr* numq_eval = evaluate(numq);
    expr_free(numq);
    if (numq_eval->type != EXPR_SYMBOL || numq_eval->data.symbol.name != SYM_True) {
        expr_free(numq_eval);
        return false;
    }
    expr_free(numq_eval);

    Expr* freeq = expr_new_function(expr_new_symbol(SYM_FreeQ), (Expr*[]){expr_copy(e), expr_new_symbol(SYM_I)}, 2);
    Expr* freeq_eval = evaluate(freeq);
    expr_free(freeq);
    if (freeq_eval->type != EXPR_SYMBOL || freeq_eval->data.symbol.name != SYM_True) {
        expr_free(freeq_eval);
        return false;
    }
    expr_free(freeq_eval);

    return true;
}

Expr* stats_quantile_point(Expr** sorted_args, size_t n, Expr* q,
                           Expr* a, Expr* b, Expr* c, Expr* d) {
    Expr* n_expr = expr_new_integer((int64_t)n);
    Expr* n_plus_b = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){n_expr, expr_copy(b)}, 2));
    Expr* times_q = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){n_plus_b, expr_copy(q)}, 2));
    Expr* h = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(a), times_q}, 2));

    double h_val = 0;
    if (!stats_is_numeric(h, &h_val, NULL)) {
        expr_free(h);
        return expr_new_symbol(SYM_Indeterminate);
    }

    if (h_val <= 1.0) {
        expr_free(h);
        return expr_copy(sorted_args[0]);
    }
    if (h_val >= (double)n) {
        expr_free(h);
        return expr_copy(sorted_args[n - 1]);
    }

    Expr* j_expr = eval_and_free(expr_new_function(expr_new_symbol(SYM_Floor), (Expr*[]){expr_copy(h)}, 1));
    int64_t j_idx = 0;
    if (j_expr->type == EXPR_INTEGER) j_idx = j_expr->data.integer;
    else j_idx = (int64_t)floor(h_val);
    expr_free(j_expr);

    if (j_idx < 1) j_idx = 1;
    if (j_idx >= (int64_t)n) j_idx = n - 1;

    Expr* j_expr2 = expr_new_integer(j_idx);
    Expr* neg_j = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), j_expr2}, 2));
    Expr* g = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(h), neg_j}, 2));
    expr_free(h);

    /* Integer h: Wolfram's Floor/Ceiling neighbors coincide -- the answer is
     * sorted_args[h-1] (== [j_idx-1], since g == 0 means h == j) for ANY c,d. */
    bool g_is_zero = (g->type == EXPR_INTEGER && g->data.integer == 0) ||
                     (g->type == EXPR_REAL && g->data.real == 0.0);
    if (g_is_zero) {
        expr_free(g);
        return expr_copy(sorted_args[j_idx - 1]);
    }

    Expr* d_times_g = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_copy(d), expr_copy(g)}, 2));
    Expr* g_weight = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(c), d_times_g}, 2));
    expr_free(g);

    Expr* neg_Aj1 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), expr_copy(sorted_args[j_idx-1])}, 2));
    Expr* diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(sorted_args[j_idx]), neg_Aj1}, 2));

    Expr* weight_diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){g_weight, diff}, 2));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(sorted_args[j_idx-1]), weight_diff}, 2));
}
